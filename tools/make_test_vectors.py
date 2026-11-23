"""Freeze small independent exhaustive-oracle cases for MATLAB checks."""
import json
from math import comb
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
from test_algorithm import brute, small_cases

cases = []
for n, rows, p, s in small_cases():
    for D in sorted({0, min(3, n), n}):
        expected = brute(n, rows, D, p, s)
        cases.append(dict(n=n, rows=rows, D=D, P=list(p), S=list(s),
                          weights=list(expected), counts=[str(v) for v in expected.values()]))
cases.append(dict(n=128, rows=[], D=24, P=[], S=[], weights=list(range(25)),
                  counts=[str(comb(128, w)) for w in range(25)]))
(ROOT / "tests/oracle_vectors.json").write_text(json.dumps(cases, indent=2) + "\n", encoding="utf-8")
print(f"Saved {len(cases)} independent cases.")
