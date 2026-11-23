"""Automatic minimum layer: independent truth, cumulative budget, and examples."""
import json
from pathlib import Path
import subprocess
import sys

import pytest

from test_algorithm import BACKENDS, brute, small_cases
from polar_spectrum import algorithm1, minimum_distance, checks_for_profile, get_profile, get_rate_matching, WorkLimitError


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("policy", ("auto", "enumerate", "dp"))
def test_minimum_exhaustive(backend, policy):
    for n, rows, p, s in small_cases():
        truth = brute(n, rows, n, p, s)
        dmin = min((w for w in truth if w > 0), default=None)
        result = minimum_distance(n, rows, punctured=p, shortened=s, backend=backend, policy=policy)
        assert result["dmin"] == dmin
        assert result["A_dmin"] == (truth[dmin] if dmin is not None else 0)
        assert result["spectrum"] == ({0: 1, dmin: truth[dmin]} if dmin is not None else {0: 1})
        assert result["cutoff"] is None
        if dmin is not None:
            assert result["search_cutoffs"][-1] == dmin
            assert all(d <= dmin for d in result["search_cutoffs"])
        else:
            assert result["search_cutoffs"] == []


@pytest.mark.parametrize("backend", BACKENDS)
def test_multiple_trials_and_shared_budget(backend):
    # Only u=0000 and u=0111 are allowed; nonzero output weight is 3.
    result = minimum_distance(4, (3, 5, 8), backend=backend)
    assert (result["dmin"], result["A_dmin"]) == (3, 1)
    assert result["search_cutoffs"] == [1, 2, 3]
    with pytest.raises(WorkLimitError):
        minimum_distance(4, (3, 5, 8), backend=backend,
                         max_work=result["stats"]["work_units"] - 1)


@pytest.mark.parametrize("backend", BACKENDS)
def test_distance_above_default_cutoff(backend):
    result = minimum_distance(16, tuple(1 << i for i in range(15)), backend=backend)
    assert (result["dmin"], result["A_dmin"], result["search_cutoffs"]) == (16, 1, [16])


def test_all_fixed_profiles():
    for n in (64, 128, 256):
        for kind in ("polar", "pac"):
            H = checks_for_profile(get_profile(n, kind))
            for mode in ("none", "puncture", "shorten", "mixed"):
                matching = get_rate_matching(n, mode)
                truth = algorithm1(n, H, 8, **matching)["spectrum"]
                dmin = min(w for w in truth if w)
                for backend in BACKENDS:
                    r = minimum_distance(n, H, **matching, backend=backend)
                    assert (r["dmin"], r["A_dmin"]) == (dmin, truth[dmin])


def test_example_cli(tmp_path):
    root = Path(__file__).resolve().parents[1]
    command = [sys.executable, str(root / "examples/python/run_spectrum.py"),
               "--kind", "pac", "--mode", "mixed", "--dmin-only"]
    output = tmp_path / "minimum.json"
    completed = subprocess.run(command + ["--json", str(output)], check=True, capture_output=True, text=True)
    assert "dmin = 6" in completed.stdout and "A_dmin = 3" in completed.stdout
    result = json.loads(output.read_text())
    assert result["dmin"] == 6 and result["A_dmin"] == "3"
    assert set(result["spectrum"]) == {"0", "6"}
    assert subprocess.run(command + ["-D", "8"], capture_output=True).returncode != 0
