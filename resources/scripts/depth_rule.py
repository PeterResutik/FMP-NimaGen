"""The depth rule: the share a point heteroplasmy must be above at its read depth, from a noise model of control samples
(resources/noise_model/).

The model holds, for each control, how many wrong bases reach each share of the reads, every 0.05% from 0.1 to 2% (with
at least 2 reads). Their mean over the controls, J(s), is turned into a density of shares: log-log interpolation
between the counted shares, the slope of the first three below 0.1%, falling linearly to 0 after the last share with
any; spread evenly over the possible wrong bases (3 per position, outside the C-stretches and the controls' own
heteroplasmies). For a sample with n reads per amplicon, a position read by c amplicons has c * n reads; the wrong
bases' reads are binomial and their number Poisson, so S = the expected wrong bases with at least
max(2, ceil(t * c * n)) reads, and the threshold is the lowest share t (to 0.01%) at which
1 - exp(-S) <= the risk (1 in --depth_risk samples). A row passes when its share is above t.
"""
import math
from functools import lru_cache

import numpy as np

SHARES = np.round(np.arange(0.001, 0.02 + 1e-9, 0.0005), 5)       # 0.1% ... 2%, every 0.05% (as fractions)
SHARE_COLUMNS = [f"{100 * s:.2f}" for s in SHARES]
C_STRETCHES = ((57, 60), (300, 315), (16180, 16193))               # not counted in the model
GENOME = 16569


def read_models(paths):
    """(counts, masked, names) of the controls in one or more model files, pooled: counts is an array of controls x
    SHARES, masked the positions left out in any of them."""
    counts, masked, names = [], set(), []
    for path in paths:
        header = None
        for line in open(path):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                if header[:3] != ["control", "median_reads_per_amplicon", "masked"] or header[3:] != SHARE_COLUMNS:
                    raise ValueError(f"{path}: not a noise model (expected the columns control, "
                                     f"median_reads_per_amplicon, masked, {SHARE_COLUMNS[0]} ... {SHARE_COLUMNS[-1]})")
                continue
            row = dict(zip(header, fields))
            names.append(row["control"])
            masked |= {int(p) for p in row["masked"].split(",") if p.strip()}
            counts.append([float(row[c]) for c in SHARE_COLUMNS])
    if not counts:
        raise ValueError(f"no controls in {', '.join(map(str, paths))}")
    return np.array(counts), masked, names


def library_segments(library_path):
    """{amplicon: [(start, end), ...]} from the FDSTOOLS library's [genome_position] block (mtNG_100 has two parts,
    across the origin)."""
    segments, in_block = {}, False
    for line in open(library_path):
        line = line.strip()
        if line.startswith("[genome_position]"):
            in_block = True
            continue
        if in_block and line.startswith("["):
            break
        if in_block and "=" in line:
            name, values = line.split("=", 1)
            numbers = [int(v) for v in values.split(",")[1:]]
            segments[name.strip()] = list(zip(numbers[0::2], numbers[1::2]))
    return segments


def counted(position, masked):
    """Whether the model counts wrong bases at the position."""
    return position not in masked and not any(a <= position <= b for a, b in C_STRETCHES)


def coverage_groups(segments, reference, masked):
    """{c: possible wrong bases at positions read by c amplicons}: 3 per counted position with an rCRS base."""
    groups = {}
    for p in range(1, GENOME + 1):
        if reference[p - 1] not in "ACGT" or not counted(p, masked):
            continue
        c = sum(1 for parts in segments.values() for a, b in parts if a <= p <= b)
        if c:
            groups[c] = groups.get(c, 0) + 3
    return groups


def density(J):
    """J(s) at SHARES (mean wrong bases per control at share s or more) -> shares sg and the expected wrong bases per
    sample at each, ws."""
    pos = J > 0
    if pos.sum() < 3:
        raise ValueError("the noise model needs wrong bases at 3 or more of its shares")
    s_pos, e_pos = SHARES[pos], J[pos]
    later = (SHARES > s_pos.max()) & (J <= 0)
    last = SHARES[later].min() if later.any() else 0.05
    slope = np.polyfit(np.log(s_pos[:3]), np.log(e_pos[:3]), 1)[0]

    def interpolated(s):
        out = np.zeros_like(s)
        lo = s < s_pos.min()
        out[lo] = e_pos[0] * (s[lo] / s_pos.min()) ** slope
        mid = (s >= s_pos.min()) & (s <= s_pos.max())
        out[mid] = np.exp(np.interp(np.log(s[mid]), np.log(s_pos), np.log(e_pos)))
        hi = (s > s_pos.max()) & (s < last)
        out[hi] = e_pos[-1] * (last - s[hi]) / (last - s_pos.max())
        return out

    sg = np.geomspace(1e-4, last, 300)
    return sg, np.clip(-np.gradient(interpolated(sg)), 0, None)


_LOG_FACTORIAL = np.zeros(1)


def _log_factorials(n):
    """log(i!) for i = 0 ... n (at least)."""
    global _LOG_FACTORIAL
    if len(_LOG_FACTORIAL) <= n:
        size = max(n + 1, 2 * len(_LOG_FACTORIAL))
        _LOG_FACTORIAL = np.array([math.lgamma(i + 1) for i in range(size)])
    return _LOG_FACTORIAL


def binomial_tail(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p), for each p of an array (all below 1). Terms beyond the largest p's mean plus
    12 standard deviations are left out (below 1e-30)."""
    p = np.asarray(p, float)
    if k <= 0:
        return np.ones_like(p)
    if k > n:
        return np.zeros_like(p)
    mean = n * p.max()
    last = min(n, int(math.ceil(mean + 12 * math.sqrt(mean) + 20)))
    if k > last:
        return np.zeros_like(p)
    lf = _log_factorials(n)
    j = np.arange(k, last + 1)
    log_choose = lf[n] - lf[j] - lf[n - j]
    log_pmf = log_choose[:, None] + j[:, None] * np.log(p)[None, :] + (n - j)[:, None] * np.log1p(-p)[None, :]
    return np.exp(log_pmf).sum(axis=0)


class DepthRule:
    """Thresholds from a pooled noise model for one amplicon panel (the FDSTOOLS library) at one risk."""

    def __init__(self, model_paths, library_path, reference, risk):
        counts, masked, self.controls = read_models(model_paths)
        self.sg, self.ws = density(counts.mean(axis=0))
        self.groups = coverage_groups(library_segments(library_path), reference, masked)
        self.possible = sum(self.groups.values())
        self.risk = risk
        self.threshold = lru_cache(maxsize=None)(self._threshold)

    def expected(self, n, t):
        """S: expected wrong bases reaching the share t (%) in a sample with n reads per amplicon."""
        S = 0.0
        for c, count in self.groups.items():
            reads = c * n
            k = max(2, math.ceil(t * reads / 100 - 1e-9))
            S += count / self.possible * (self.ws * binomial_tail(k, reads, self.sg)).sum()
        return S

    def _threshold(self, n):
        """The share (%) a point heteroplasmy must be above at n reads per amplicon."""
        n = max(int(n), 2)
        lo, hi = 0.0, 100.0
        while hi - lo > 0.01:
            mid = (lo + hi) / 2
            if 1 - math.exp(-self.expected(n, mid)) <= self.risk:
                hi = mid
            else:
                lo = mid
        return round(hi, 2)
