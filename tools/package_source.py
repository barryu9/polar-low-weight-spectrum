"""Create a portable source ZIP, excluding local builds, caches, and git metadata."""
import argparse
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT.parent / "polar-low-weight-spectrum-v0.1.1.zip")
    args = parser.parse_args()
    destination = args.output.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    excluded = {".git", "build", "dist", "__pycache__", ".pytest_cache", ".venv"}
    with zipfile.ZipFile(destination, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(ROOT.rglob("*")):
            relative = path.relative_to(ROOT)
            if any(p in excluded or p.startswith("build-") or p.endswith(".egg-info") for p in relative.parts):
                continue
            if not path.is_file() or path.resolve() == destination or path.suffix in (".pyc", ".log", ".asv"):
                continue
            archive.write(path, Path(ROOT.name) / relative)
    print(destination)


if __name__ == "__main__":
    main()
