"""Reproduce Python/C++ exact agreement and timings for the 24 fixed cases."""
import argparse
import json
from pathlib import Path
import platform
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from polar_spectrum import algorithm1, checks_for_profile, get_profile, get_rate_matching


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeat", type=int, default=7)
    parser.add_argument("-D", "--cutoff", type=int, default=8)
    parser.add_argument("--output", type=Path, default=ROOT / "results/python_cpp_validation.json")
    args = parser.parse_args()
    if args.repeat < 1:
        parser.error("--repeat must be positive")
    report = dict(platform=platform.platform(), processor=platform.processor(),
                  python=sys.version, repeat=args.repeat, cutoff=args.cutoff, cases=[])
    for n in (64, 128, 256):
        for kind in ("polar", "pac"):
            H = checks_for_profile(get_profile(n, kind))
            for mode in ("none", "puncture", "shorten", "mixed"):
                matching = get_rate_matching(n, mode)
                reference = None
                for backend in ("python", "cpp"):
                    algorithm1(n, H, args.cutoff, backend=backend, **matching)  # warm-up
                    runs = [algorithm1(n, H, args.cutoff, backend=backend, **matching) for _ in range(args.repeat)]
                    r = runs[-1]
                    assert all(x["spectrum"] == r["spectrum"] for x in runs)
                    if reference is not None:
                        assert r["spectrum"] == reference
                    reference = r["spectrum"]
                    r.update(kind=kind, mode=mode,
                             median_elapsed_seconds=statistics.median(x["elapsed_seconds"] for x in runs),
                             minimum_elapsed_seconds=min(x["elapsed_seconds"] for x in runs))
                    r["spectrum"] = {str(w): str(a) for w, a in r["spectrum"].items()}
                    report["cases"].append(r)
    matlab_file = ROOT / "results/matlab_validation.json"
    if matlab_file.exists():
        matlab = json.loads(matlab_file.read_text(encoding="utf-8"))
        compared = 0
        for m in matlab["cases"]:
            if m["cutoff"] != args.cutoff:
                continue
            reference = next(r for r in report["cases"] if (r["n"],r["kind"],r["mode"],r["backend"]) ==
                             (m["n"],m["kind"],m["mode"],"python"))
            weights = m["weights"] if isinstance(m["weights"], list) else [m["weights"]]
            counts = m["counts"] if isinstance(m["counts"], list) else [m["counts"]]
            assert dict(zip(map(str, weights), counts)) == reference["spectrum"]
            for field in ("original_dimension","shortened_dimension","transmitted_dimension","kernel_dimension"):
                assert m[field] == reference[field]
            compared += 1
        report["matlab_crosschecked_cases"] = compared
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: 24 cases x 2 backends, {args.repeat} repetitions; MATLAB crosschecks: {report.get('matlab_crosschecked_cases', 0)}.")
    for backend in ("python", "cpp"):
        times = [r["median_elapsed_seconds"] for r in report["cases"] if r["backend"] == backend]
        print(f"{backend}: median-case range {min(times):.6f} .. {max(times):.6f} seconds")


if __name__ == "__main__":
    main()
