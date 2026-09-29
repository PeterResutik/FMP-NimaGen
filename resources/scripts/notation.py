"""How a molecule differs from rCRS, in EMPOP-style labels (the general rule).

Used outside the C-stretches 16180-16193 and 300-315 (which have their own
frames) and for rebuilding a caller's major calls into one spelling:

- fewest changes;
- on an equal count, substitutions beat gaps, and gaps that lengthen or
  shorten a run beat other gaps;
- a gap is one block, placed as far 3' as the sequence allows, also across
  repeat copies;
- in a tandem repeat whole copies are added or removed at the 3'-most copy
  (AC in 514-524 gives A523- C524-); with a unit of 3+ bases this holds even
  if it costs one extra substitution;
- the rCRS N at 3107 is left out, so a base the reads put there is an
  insertion placed by the same rules (N3107T becomes -3109.1T).

A window must lie within 1-16569; nothing is shifted across the origin.
"""
import re

SUB, SLIP, GAP, OPEN = 100, 101, 102, 5  # every change costs ~100, so fewest changes wins first
INF = float("inf")


def reference_window(reference, start, end):
    """[(position, base)] of rCRS start..end (1-based), without the N at 3107."""
    return [(p, reference[p - 1]) for p in range(start, end + 1) if reference[p - 1] != "N"]


def align(molecule, ref):
    """Columns (position or None, ref base or '-', molecule base or '-') of the
    cheapest alignment, gaps moved as far 3' as the sequence allows."""
    n, m = len(ref), len(molecule)
    rb = [b for _, b in ref]

    def deletion_cost(i):
        b = rb[i - 1]
        return SLIP if (i >= 2 and rb[i - 2] == b) or (i < n and rb[i] == b) else GAP

    def insertion_cost(i, j):
        b = molecule[j - 1]
        near = {rb[i - 1] if i >= 1 else None, rb[i] if i < n else None,
                molecule[j - 2] if j >= 2 else None, molecule[j] if j < m else None}
        return SLIP if b in near else GAP

    M = [[INF] * (m + 1) for _ in range(n + 1)]  # last column a match or substitution
    X = [[INF] * (m + 1) for _ in range(n + 1)]  # last column a deletion of ref[i-1]
    Y = [[INF] * (m + 1) for _ in range(n + 1)]  # last column an insertion of molecule[j-1]
    M[0][0] = 0
    for i in range(n + 1):
        for j in range(m + 1):
            if i and j:
                M[i][j] = min(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]) + (SUB if rb[i - 1] != molecule[j - 1] else 0)
            if i:
                X[i][j] = min(M[i - 1][j] + OPEN, X[i - 1][j], Y[i - 1][j] + OPEN) + deletion_cost(i)
            if j:
                Y[i][j] = min(M[i][j - 1] + OPEN, Y[i][j - 1], X[i][j - 1] + OPEN) + insertion_cost(i, j)
    cols, i, j = [], n, m
    state = min((X[n][m], 0, "X"), (Y[n][m], 1, "Y"), (M[n][m], 2, "M"))[2]  # gaps first on ties: 3' placement
    while i or j:
        if state == "M":
            cols.append((ref[i - 1][0], rb[i - 1], molecule[j - 1]))
            prev = M[i][j] - (SUB if rb[i - 1] != molecule[j - 1] else 0)
            i -= 1
            j -= 1
            state = "X" if X[i][j] == prev else "Y" if Y[i][j] == prev else "M"
        elif state == "X":
            cols.append((ref[i - 1][0], rb[i - 1], "-"))
            prev = X[i][j] - deletion_cost(i)
            i -= 1
            state = "X" if X[i][j] == prev else "M" if M[i][j] + OPEN == prev else "Y"
        else:
            cols.append((None, "-", molecule[j - 1]))
            prev = Y[i][j] - insertion_cost(i, j)
            j -= 1
            state = "Y" if Y[i][j] == prev else "M" if M[i][j] + OPEN == prev else "X"
    cols.reverse()
    return _shift_gaps_3prime(cols)


def _shift_gaps_3prime(cols):
    """Move every gap block one column to the right while the sequence allows it."""
    cols = [list(c) for c in cols]
    changed = True
    while changed:
        changed = False
        k = 0
        while k < len(cols):
            if cols[k][2] == "-" or cols[k][1] == "-":
                kind = 2 if cols[k][2] == "-" else 1  # 2: deletion (molecule gap), 1: insertion (ref gap)
                other = 3 - kind
                e = k
                while e + 1 < len(cols) and cols[e + 1][kind] == "-" and cols[e + 1][other] != "-":
                    e += 1
                c = e + 1
                match_next = c < len(cols) and cols[c][1] != "-" and cols[c][2] != "-" and cols[c][1] == cols[c][2]
                start = None
                if match_next and cols[k][other] == cols[c][other]:
                    start = k  # the whole block can move
                elif match_next:
                    # a mixed block (e.g. inserted TTC): its last run moves on its own
                    s2 = e
                    while s2 - 1 >= k and cols[s2 - 1][other] == cols[e][other]:
                        s2 -= 1
                    if s2 > k and cols[s2][other] == cols[c][other]:
                        start = s2
                if start is not None:
                    if kind == 2:  # deletion: ref bases stay, the molecule base moves left
                        cols[start][2], cols[c][2] = cols[c][2], "-"
                    else:  # insertion: molecule bases stay, the ref base (and its position) moves left
                        cols[start][0], cols[start][1], cols[c][0], cols[c][1] = cols[c][0], cols[c][1], None, "-"
                    changed = True
                k = e + 1
            else:
                k += 1
    return [tuple(c) for c in cols]


def labels_from_columns(cols):
    """Labels of an alignment: A73G, C8281-, -315.1C (inserted after the last ref position)."""
    labels, last_pos, k_ins = set(), None, 0
    for pos, r, q in cols:
        if pos is not None:
            last_pos, k_ins = pos, 0
            if q == "-":
                labels.add(f"{r}{pos}-")
            elif q != r:
                labels.add(f"{r}{pos}{q}")
        else:
            k_ins += 1
            labels.add(f"-{last_pos}.{k_ins}{q}")
    return labels


def _tandem_repeats(ref):
    """(start index, unit length, copies) of tandem repeats with a unit of 2+ bases."""
    seq = "".join(b for _, b in ref)
    out = []
    for unit in range(2, 13):
        i = 0
        while i + 2 * unit <= len(seq):
            if seq[i:i + unit] == seq[i + unit:i + 2 * unit] and len(set(seq[i:i + unit])) > 1:
                n = 2
                while seq[i + n * unit:i + (n + 1) * unit] == seq[i:i + unit]:
                    n += 1
                out.append((i, unit, n))
                i += n * unit
            else:
                i += 1
    return out


def _repeat_readings(ref, molecule):
    """Readings with whole repeat copies removed or added at the 3'-most copy,
    the rest written as substitutions (at most one; none for a 2-base unit)."""
    seq = "".join(b for _, b in ref)
    readings = []
    for i, unit, n in _tandem_repeats(ref):
        max_subs = 0 if unit == 2 else 1
        diff = len(seq) - len(molecule)
        if diff > 0 and diff % unit == 0 and diff // unit < n:
            c0, c1 = i + (n - diff // unit) * unit, i + n * unit - 1  # the last copies
            while c1 + 1 < len(seq) and seq[c0] == seq[c1 + 1]:
                c0 += 1
                c1 += 1  # as far 3' as the sequence allows
            rest = [j for j in range(len(seq)) if not c0 <= j <= c1]
            subs = {f"{seq[j]}{ref[j][0]}{molecule[t]}" for t, j in enumerate(rest) if seq[j] != molecule[t]}
            if len(subs) <= max_subs:
                readings.append({f"{seq[j]}{ref[j][0]}-" for j in range(c0, c1 + 1)} | subs)
        elif diff < 0 and (-diff) % unit == 0:
            ins = seq[i:i + unit] * (-diff // unit)
            after = i + n * unit - 1
            while after + 1 < len(seq) and ins[0] == seq[after + 1]:
                ins = ins[1:] + seq[after + 1]
                after += 1  # as far 3' as the sequence allows
            extended = seq[:after + 1] + ins + seq[after + 1:]
            positions = [ref[j][0] for j in range(after + 1)] + [None] * len(ins) + [ref[j][0] for j in range(after + 1, len(seq))]
            subs, ok = set(), True
            for t, (e, q) in enumerate(zip(extended, molecule)):
                if e != q:
                    if positions[t] is None:
                        ok = False
                        break
                    subs.add(f"{e}{positions[t]}{q}")
            if ok and len(subs) <= max_subs:
                readings.append({f"-{ref[after][0]}.{t + 1}{b}" for t, b in enumerate(ins)} | subs)
    return readings


def label_position(label):
    """Sort key: position, then insertion index (-315.1C after C315T)."""
    m = re.match(r"^-?[ACGTN]?(\d+)(?:\.(\d+))?", label)
    return (int(m.group(1)), int(m.group(2) or 0))


def describe(molecule, reference, start, end):
    """Labels for a molecule covering rCRS start..end (without a base for the N at 3107)."""
    ref = reference_window(reference, start, end)
    labels = labels_from_columns(align(molecule, ref))
    for reading in _repeat_readings(ref, molecule):
        if len(reading) <= len(labels) + 1 and reading != labels:
            labels = reading
            break
    return sorted(labels, key=label_position)
