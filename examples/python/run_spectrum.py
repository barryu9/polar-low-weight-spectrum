"""Run the paper's Algorithm 1 with a frozen GA profile and rate-matching pattern."""
import argparse
import json
from pathlib import Path
import statistics
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
from polar_spectrum import algorithm1, checks_for_profile, get_profile, get_rate_matching


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, choices=(64, 128, 256), default=64)
    parser.add_argument("--kind", choices=("polar", "pac"), default="polar")
    parser.add_argument("--backend", choices=("python", "cpp"), default="python")
    parser.add_argument("--mode", choices=("none", "puncture", "shorten", "mixed"), default="none")
    target = parser.add_mutually_exclusive_group()
    target.add_argument("-D", "--cutoff", type=int, default=8)
    target.add_argument("--dmin-only", action="store_true",
                        help="Automatically find dmin and A_dmin; stop at the first positive weight.")
    parser.add_argument("--policy", choices=("auto", "enumerate", "dp"), default="auto")
    parser.add_argument("--max-work", type=int, default=10_000_000)
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.repeat < 1:
        parser.error("--repeat must be positive")
    cfg = get_profile(args.n, args.kind)
    H = checks_for_profile(cfg)
    matching = get_rate_matching(args.n, args.mode)
    cutoff = None if args.dmin_only else args.cutoff
    runs = [algorithm1(args.n, H, cutoff, backend=args.backend, policy=args.policy,
                       max_work=args.max_work, **matching) for _ in range(args.repeat)]
    result = runs[-1]
    result.update(kind=args.kind, mode=args.mode,
                  median_elapsed_seconds=statistics.median(r["elapsed_seconds"] for r in runs))
    print(f"{args.kind.upper()} mother ({args.n},{args.n//2}), GA Eb/N0=2 dB, {args.mode}, {args.backend}")
    if args.dmin_only:
        print(f"Transmitted dimension={result['transmitted_dimension']}; automatic dmin; exact={result['exact']}")
        print(f"dmin = {result['dmin'] if result['dmin'] is not None else 'undefined (zero-dimensional code)'}")
        print(f"A_dmin = {result['A_dmin']}")
    else:
        print(f"Transmitted dimension={result['transmitted_dimension']}; D={args.cutoff}; exact={result['exact']}")
        print("weight   multiplicity")
        for weight, count in result["spectrum"].items():
            print(f"{weight:6d}   {count}")
    print(f"Algorithm time: {result['median_elapsed_seconds']:.6f} s (median of {args.repeat})")
    if "binding_elapsed_seconds" in result:
        print(f"Last native call including binding: {result['binding_elapsed_seconds']:.6f} s")
    if args.json:
        # Decimal strings preserve exact counts in downstream JSON consumers.
        result["spectrum"] = {str(k): str(v) for k, v in result["spectrum"].items()}
        if args.dmin_only:
            result["A_dmin"] = str(result["A_dmin"])
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
