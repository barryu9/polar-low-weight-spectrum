"""Build the optional standard-library-only C++ kernel without CMake."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default=os.environ.get("CXX"))
    parser.add_argument("--output", type=Path, default=ROOT / "build")
    args = parser.parse_args()
    compiler = args.compiler or shutil.which("g++") or shutil.which("clang++")
    if not compiler and sys.platform == "win32":
        candidates = [Path(os.environ.get("MW_MINGW64_LOC", "C:/MinGW64")) / "bin/g++.exe"]
        compiler = next((str(p) for p in candidates if p.exists()), None)
    if not compiler:
        parser.error("Install g++/clang++, use --compiler, or build with CMake and MSVC.")
    args.output.mkdir(parents=True, exist_ok=True)
    name = "polar_spectrum.dll" if sys.platform == "win32" else "libpolar_spectrum.dylib" if sys.platform == "darwin" else "libpolar_spectrum.so"
    command = [compiler, "-O3", "-DNDEBUG", "-std=c++17", "-shared",
               str(ROOT / "cpp/src/kernel.cpp"), "-o", str(args.output / name)]
    command += ["-static-libgcc", "-static-libstdc++"] if sys.platform == "win32" else ["-fPIC"]
    subprocess.run(command, check=True)
    print(args.output / name)


if __name__ == "__main__":
    main()
