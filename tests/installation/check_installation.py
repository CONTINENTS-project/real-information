"""Build/install checks. Run on a workstation or an approved compute allocation.

Creates a fresh output directory and retains all output for inspection. Builds
run sequentially with at most two compiler processes and two OpenMP threads.
"""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys


def run(command, *, cwd, env):
    print("Running:", " ".join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), cwd=cwd, env=env, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, required=True,
                        help="fresh output directory; existing directories are refused")
    args = parser.parse_args()
    output = args.work_dir.resolve()
    source = Path(__file__).resolve().parents[2]
    fixtures = source / "tests"
    consumer = fixtures / "installation"
    cmake = shutil.which("cmake")
    ctest = shutil.which("ctest")
    if not cmake or not ctest:
        parser.error("cmake and ctest must be available")
    # Fail before any build if this interpreter cannot run the Python consumer.
    import numpy  # noqa: F401

    output.mkdir(parents=True, exist_ok=False)
    environment = os.environ.copy()
    environment.pop("REAL_INFO_LIB_PATH", None)
    environment["OMP_NUM_THREADS"] = "2"
    environment["PYTHONDONTWRITEBYTECODE"] = "1"

    for label, build_type, libdir, testing in (
        ("release-lib", "Release", "lib", "OFF"),
        ("debug-lib64", "Debug", "lib64", "ON"),
    ):
        case = output / label
        case.mkdir()
        build = case / "build"
        prefix = case / "prefix"
        configure = [cmake, "-S", source, "-B", build,
                     f"-DCMAKE_BUILD_TYPE={build_type}", f"-DBUILD_TESTING={testing}",
                     f"-DCMAKE_INSTALL_PREFIX={prefix}", f"-DCMAKE_INSTALL_LIBDIR={libdir}"]
        if testing == "ON":
            configure.append(f"-DPython3_EXECUTABLE={sys.executable}")
        run(configure, cwd=case, env=environment)
        run([cmake, "--build", build, "--target", "bit_info", "--parallel", "2"],
            cwd=case, env=environment)
        run([cmake, "--install", build], cwd=case, env=environment)
        if testing == "ON":
            run([ctest, "--output-on-failure", "-R", "^py_test_library_discovery$"],
                cwd=build, env=environment)

        relative_library = next(line.split('"')[1] for line in
                                (prefix / "real_info.py").read_text().splitlines()
                                if line.startswith("_LIBRARY_RELATIVE_PATH ="))
        for relocated in (False, True):
            if relocated:
                moved = case / "relocated prefix"
                prefix.rename(moved)
                prefix = moved
            python_env = dict(environment, PYTHONPATH=str(prefix))
            run([sys.executable, "-B", consumer / "python_consumer.py",
                 prefix, relative_library, fixtures / "fixed_point.yaml"],
                cwd=case, env=python_env)

        # Exercise a distinct explicit library override, then a missing override.
        override = case / "override.so"
        shutil.copyfile(prefix / relative_library, override)
        run([sys.executable, "-B", consumer / "python_consumer.py",
             prefix, relative_library, fixtures / "fixed_point.yaml"], cwd=case,
            env=dict(python_env, REAL_INFO_LIB_PATH=str(override)))
        missing = case / "missing.so"
        failure = subprocess.run([sys.executable, "-B", "-c", "import real_info"],
                                 cwd=case, env=dict(python_env, REAL_INFO_LIB_PATH=str(missing)),
                                 text=True, capture_output=True)
        if failure.returncode == 0 or "FileNotFoundError" not in failure.stderr \
                or str(missing) not in failure.stderr:
            raise RuntimeError(f"Missing override did not fail as expected: {failure.stderr}")

        # The old prefix no longer exists, so source/prefix leaks cannot hide here.
        package_dir = prefix / libdir / "cmake" / "real_information"
        for name in ("real_informationTargets.cmake", "real_informationConfig.cmake"):
            text = (package_dir / name).read_text()
            if str(source) in text or str(case / "prefix") in text:
                raise RuntimeError(f"Installed package contains a fixed source/prefix path: {name}")
        consumer_build = case / "consumer-build"
        run([cmake, "-S", consumer, "-B", consumer_build,
             f"-DCMAKE_BUILD_TYPE={build_type}", f"-DCMAKE_PREFIX_PATH={prefix}"],
            cwd=case, env=environment)
        run([cmake, "--build", consumer_build, "--parallel", "2"], cwd=case, env=environment)
        run([ctest, "--output-on-failure"], cwd=consumer_build, env=environment)

    print(f"Installation checks passed; outputs retained in {output}", flush=True)


if __name__ == "__main__":
    main()
