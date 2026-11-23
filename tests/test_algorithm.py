"""Independent exhaustive oracle: encode via explicit Kronecker generator."""
from collections import Counter
from math import comb
from pathlib import Path
import csv
import json
import random
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from polar_spectrum import algorithm1, WorkLimitError, get_profile, checks_for_profile, get_rate_matching
from polar_spectrum.native import library_path
from polar_spectrum.algorithm import Budget, Family, Solver
from polar_spectrum.profiles import selection_matrices

try:
    library_path()
    BACKENDS = ("python", "cpp")
except FileNotFoundError:
    BACKENDS = ("python",)


def encode(u, n):
    # G[i,j]=1 iff the binary expansion of j is a subset of i.
    return sum((sum((u >> i & 1) for i in range(n) if (i & j) == j) % 2) << j for j in range(n))


def brute(n, checks, D, p=(), s=()):
    kept = sum(1 << j for j in range(n) if j not in p and j not in s)
    smask = sum(1 << j for j in s)
    distinct = set()
    for u in range(1 << n):
        if any((u & row).bit_count() % 2 for row in checks):
            continue
        x = encode(u, n)
        if not x & smask:
            distinct.add(x & kept)
    return dict(sorted(Counter(x.bit_count() for x in distinct if x.bit_count() <= D).items()))


def small_cases():
    rng = random.Random(20260907)
    for n in (1, 2, 4, 8):
        for _ in range(16):
            # Noncanonical, redundant checks exercise row reduction too.
            rows = [rng.randrange(1 << n) for _ in range(rng.randrange(n + 2))]
            labels = [rng.randrange(3) for _ in range(n)]
            p = tuple(j for j in range(n) if labels[j] == 1)
            s = tuple(j for j in range(n) if labels[j] == 2)
            yield n, rows, p, s
    # Causal dynamic checks with genuine lower-half Gen expansion.
    for _ in range(12):
        n = 16
        rows = [(1 << i) | rng.randrange(1 << i) for i in range(n) if rng.random() < .6]
        yield n, rows, (0, 2), (14, 15)


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("policy", ("auto", "enumerate", "dp"))
def test_independent_exhaustive(backend, policy):
    for n, rows, p, s in small_cases():
        for D in sorted({0, min(3, n), n}):
            expected = brute(n, rows, D, p, s)
            actual = algorithm1(n, rows, D, punctured=p, shortened=s, backend=backend, policy=policy)
            assert actual["spectrum"] == expected, (n, rows, D, p, s, backend, policy)
            if D == n:
                assert sum(expected.values()) == 1 << actual["transmitted_dimension"]


@pytest.mark.parametrize("backend", BACKENDS)
def test_exact_large_integers(backend):
    # 128 choose 24 is beyond both 2^53 and 2^64.
    result = algorithm1(128, (), 24, backend=backend)
    assert result["spectrum"] == {w: comb(128, w) for w in range(25)}


@pytest.mark.parametrize("backend", BACKENDS)
def test_projection_and_shortening(backend):
    assert algorithm1(8, (), 8, punctured=range(8), backend=backend)["spectrum"] == {0: 1}
    r = algorithm1(8, (), 8, punctured=(0, 1), shortened=(6, 7), backend=backend)
    assert r["spectrum"] == {w: comb(4, w) for w in range(5)}
    assert (r["original_dimension"], r["shortened_dimension"], r["transmitted_dimension"], r["kernel_dimension"]) == (8, 6, 4, 2)


@pytest.mark.parametrize("backend", BACKENDS)
def test_work_limit(backend):
    with pytest.raises(WorkLimitError):
        algorithm1(64, checks_for_profile(get_profile()), 8, backend=backend, max_work=1)


@pytest.mark.parametrize("kwargs", [
    {"n": 3}, {"D": -1}, {"checks": [256]}, {"punctured": [0, 0]},
    {"punctured": [0], "shortened": [0]}, {"policy": "invalid"}, {"max_work": 0},
    {"D": True}, {"punctured": [8]}, {"backend": "unknown"},
])
def test_validation(kwargs):
    args = dict(n=8, checks=(), D=2)
    args.update(kwargs)
    with pytest.raises(ValueError):
        algorithm1(**args)


def test_pac_convolution_definition():
    cfg = dict(n=8, information=(3, 5, 6, 7), taps=(1, 0, 1, 1, 0, 1, 1))
    H = checks_for_profile(cfg)
    words = set()
    for bits in range(16):
        v = [0] * 8
        for k, i in enumerate(cfg["information"]):
            v[i] = bits >> k & 1
        u = sum((sum(cfg["taps"][lag] * v[i-lag] for lag in range(min(i+1, 7))) % 2) << i for i in range(8))
        assert not any((u & h).bit_count() % 2 for h in H)
        words.add(encode(u, 8))
    for backend in BACKENDS:
        assert algorithm1(8, H, 8, backend=backend)["spectrum"] == dict(Counter(x.bit_count() for x in words))


@pytest.mark.parametrize("n", (64, 128, 256))
@pytest.mark.parametrize("kind", ("polar", "pac"))
@pytest.mark.parametrize("mode", ("none", "puncture", "shorten", "mixed"))
def test_fixed_profiles(n, kind, mode):
    H = checks_for_profile(get_profile(n, kind))
    matching = get_rate_matching(n, mode)
    reference = algorithm1(n, H, 8, **matching)
    for backend in BACKENDS:
        assert algorithm1(n, H, 8, backend=backend, **matching)["spectrum"] == reference["spectrum"]


def test_paper_affine_example():
    # Paper's N=8 example: a=(0,1,0,0), s=(1,1,0,0), b3=b0.
    for policy in ("auto", "enumerate", "dp"):
        family = Family(2, 3, (1,) * 8, [(0b0111, 0)], 4, Budget(10000))
        assert family.histogram(policy) == {2: 2, 4: 4}


def test_weighted_wef_pair():
    rng = random.Random(13)
    for n in (2, 4, 8):
        for _ in range(12):
            q = tuple(rng.randrange(5) for _ in range(n))
            length = rng.randrange(n)
            prefix = tuple(rng.randrange(2) for _ in range(length))
            solver = Solver(n, {}, 5, "auto", Budget(100000))
            pair = solver.wef_pair(q, prefix)
            for bit in (0, 1):
                p = prefix + (bit,)
                expected = Counter()
                for tail in range(1 << (n-len(p))):
                    u = sum(b << j for j, b in enumerate(p)) | (tail << len(p))
                    x = encode(u, n)
                    weight = sum(q[j] for j in range(n) if x >> j & 1)
                    if weight <= 5:
                        expected[weight] += 1
                assert pair[bit] == expected


def test_fixture_provenance_and_matrices():
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / "tools"))
    from freeze_profiles import ga
    fixtures = json.loads((root / "profiles/ga_2db.json").read_text())
    for n in (64, 128, 256):
        info, _ = ga(n, n//2)
        for kind in ("polar", "pac"):
            cfg = get_profile(n, kind)
            frozen = fixtures[f"{kind}_{n}_{n//2}"]
            assert list(cfg["information"]) == info == frozen["information"]
            assert list(cfg["taps"]) == frozen["taps"]
        for mode in ("none", "puncture", "shorten", "mixed"):
            ET, EP, ES = selection_matrices(n, mode)
            assert sorted(row.index(1) for row in ET + EP + ES) == list(range(n))
            for name, matrix in (("punctured", EP), ("shortened", ES)):
                with (root / f"profiles/matrices/{n}_{mode}_{name}.csv").open(newline="") as stream:
                    assert [[int(x) for x in row] for row in csv.reader(stream)] == matrix
