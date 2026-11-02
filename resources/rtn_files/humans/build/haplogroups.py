"""One sequence per mitoLEAF haplogroup, for build_humans.py.

Each sequence is rCRS with the haplogroup's motif applied:
- substitutions (73G), deletions per base (8281-) and insertions (573.1C) as written;
- lowercase marks (315.1c, 16193c: a base that may be missing) as present;
- a mark that allows two or more states (an IUPAC code, 16519Y) gets the state most of the
  haplogroup's genomes carry. A haplogroup without genomes takes the genomes of its nearest
  ancestral clade (the ancestor and all its descendants). A tie goes to the parent's clade,
  counting only the tied states, and further up while it stays tied; rCRS only if the whole
  tree is tied.
The genomes are the GenBank genomes and re-called people kept in the list, with the haplogroup
scan_humans.py gave them. A genome's ambiguous bases are not counted.
"""
import re
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np

IUPAC = {"R": "AG", "Y": "CT", "M": "AC", "K": "GT", "S": "CG", "W": "AT",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
TOKEN = re.compile(r"(\d+)(?:\.(\d+))?([A-Za-z-])")
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
REF = None


def rcrs_pos(p):
    """Position on rCRS without its N at 3107 -> rCRS numbering."""
    return p if p <= 3106 else p + 1


def _init(ref):
    global REF
    REF = ref


def record_tokens(line):
    """A genome's differences from rCRS in mitoLEAF's notation and the stretch it covers, from its
    alignment (as scan_humans.py: deletions per base and insertions placed 3')."""
    f = line.split("\t")
    flag = int(f[1])
    if flag & 4 or flag & 256:
        return f[0], [], None
    start, seq = int(f[3]), f[9]
    r, q, toks = start - 1, 0, []
    for n, op in CIGAR.findall(f[5]):
        n = int(n)
        if op == "M":
            a = np.frombuffer(REF[r:r + n].encode(), dtype=np.uint8)
            b = np.frombuffer(seq[q:q + n].encode(), dtype=np.uint8)
            for i in np.nonzero(a != b)[0]:
                toks.append(f"{rcrs_pos(r + i + 1)}{seq[q + i]}")
            r, q = r + n, q + n
        elif op == "I":
            p, s = r, seq[q:q + n]
            while p < len(REF) and REF[p] == s[0]:
                s, p = s[1:] + s[0], p + 1
            toks += [f"{rcrs_pos(p)}.{k + 1}{c}" for k, c in enumerate(s)]
            q += n
        elif op == "D":
            p, e = r + 1, r + n
            while e < len(REF) and REF[e] == REF[p - 1]:
                p, e = p + 1, e + 1
            toks += [f"{rcrs_pos(x)}-" for x in range(p, e + 1)]
            r += n
        elif op == "S":
            q += n
    return f[0], toks, (rcrs_pos(start), rcrs_pos(r))


def genome_states(sams, names, ref):
    """{genome: ({rCRS position: base or "-"}, [covered stretches], {positions with an ambiguous base})}
    for the named genomes, from their alignments (SAM files)."""
    out = {}
    for sam in sams:
        with open(sam) as fh, Pool(8, initializer=_init, initargs=(ref,)) as pool:
            lines = (l for l in fh if l[0] != "@" and l.split("\t", 1)[0] in names)
            for name, toks, iv in pool.imap(record_tokens, lines, chunksize=200):
                bases, cover, amb = out.setdefault(name, ({}, [], set()))
                if iv:
                    cover.append(iv)
                for t in toks:
                    m = TOKEN.fullmatch(t)
                    if m[2]:
                        continue
                    pos, s = int(m[1]), m[3]
                    if s == "-" or s in "ACGT":
                        bases[pos] = s
                    else:
                        amb.add(pos)
    return out


def parse(motif):
    """The motif as (substitutions {pos: base}, deletions {pos}, insertions {pos: {k: base}},
    two-state marks {pos: IUPAC code})."""
    subs, dels, ins, two = {}, set(), defaultdict(dict), {}
    for t in motif.split():
        m = TOKEN.fullmatch(t)
        if not m:
            raise ValueError(f"unknown motif token {t}")
        pos, k, s = int(m[1]), m[2], m[3]
        if s == "-":
            dels.add(pos)
        elif k:
            ins[pos][int(k)] = s.upper()          # lowercase insertions (315.1c) are written as present
        elif s.upper() in "ACGT":
            subs[pos] = s.upper()                 # lowercase substitutions (16193c) too
        else:
            two[pos] = s
    return subs, dels, ins, two


class Tree:
    """mitoLEAF's tree (tree.json): parents, children and the genomes of each haplogroup."""

    def __init__(self, tree_json, genome_hg):
        self.parent, self.children = {}, defaultdict(list)

        def walk(n, p):
            self.parent[n["name"]] = p
            if p:
                self.children[p].append(n["name"])
            for c in n.get("children") or []:
                walk(c, n["name"])
        walk(tree_json, None)
        self.members = defaultdict(list)
        for g, h in genome_hg.items():
            self.members[h].append(g)
        self._clade = {}

    def clade(self, h):
        if h not in self._clade:
            out, stack = [], [h]
            while stack:
                x = stack.pop()
                out += self.members[x]
                stack += self.children[x]
            self._clade[h] = out
        return self._clade[h]


_COUNTS = {}


def decide(h, pos, tree, states, rcrs):
    """(state, where it was decided, counts) for a two-state mark at pos in haplogroup h."""
    def count(node, whole_clade):
        key = (node, whole_clade, pos)
        if key not in _COUNTS:
            c = Counter()
            for g in (tree.clade(node) if whole_clade else tree.members[node]):
                bases, cover, amb = states[g]
                if pos not in amb and any(a <= pos <= b for a, b in cover):
                    c[bases.get(pos, rcrs[pos - 1])] += 1
            _COUNTS[key] = c
        return _COUNTS[key]
    src, c = h, count(h, False)
    while not c and tree.parent[src]:
        src = tree.parent[src]
        c = count(src, True)
    counts = c
    if not c:
        return rcrs[pos - 1], "no genomes", counts
    top = max(c.values())
    tied = sorted(s for s, n in c.items() if n == top)
    where = "own" if src == h else f"clade {src}"
    while len(tied) > 1 and tree.parent[src]:
        src = tree.parent[src]
        cc = count(src, True)
        sub = {s: cc.get(s, 0) for s in tied}
        top = max(sub.values())
        tied = sorted(s for s, n in sub.items() if n == top)
        where = f"tie, clade {src}"
    if len(tied) > 1:
        return (rcrs[pos - 1] if rcrs[pos - 1] in tied else tied[0]), "tie at the root", counts
    return tied[0], where, counts


def build(motifs, tree, states, rcrs):
    """{haplogroup: sequence (without the N at 3107)}, {haplogroup: {rCRS position: base, "" if deleted}}
    for the positions the motif changes, and one row per two-state mark (for haplogroup_marks.tsv)."""
    seqs, changed, rows = {}, {}, []
    for h, motif in motifs.items():
        subs, dels, ins, two = parse(motif)
        for pos, code in sorted(two.items()):
            state, where, counts = decide(h, pos, tree, states, rcrs)
            rows.append([h, f"{pos}{code}", where, " ".join(f"{s}:{n}" for s, n in counts.most_common()), state])
            if state == "-":
                dels.add(pos)
            else:
                subs[pos] = state
        out = []
        for pos in range(1, len(rcrs) + 1):
            if pos != 3107 and pos not in dels:
                out.append(subs.get(pos, rcrs[pos - 1]))
            out.extend(ins[pos][k] for k in sorted(ins[pos]))
        seqs[h] = "".join(out)
        changed[h] = {**{p: b for p, b in subs.items()}, **{p: "" for p in dels}}
    return seqs, changed, rows
