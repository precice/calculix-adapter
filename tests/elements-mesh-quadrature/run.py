import os
from itertools import product
from math import dist, sqrt
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
TOLERANCE = 1e-12


def check_points(workdir):
    a = 1 / sqrt(3)
    references = {
        "Tets": [(11, 1, 0.125)],
        "Hexes": [(x - 3, y + 1, z + 1) for x, y, z in product((-a, a), repeat=3)],
    }
    success = True
    for mesh, expected in references.items():
        points = [
            tuple(map(float, line.split()))
            for line in (workdir / f"{mesh}.txt").read_text().splitlines()
        ]
        error = max(
            min((dist(ref, point) for point in points), default=float("inf"))
            for ref in expected
        )
        print(
            f"{mesh}: received {len(points)}, expected {len(expected)}. "
            f"Maximum distance from a reference to its nearest received point: {error:.3g}."
        )
        if len(points) != len(expected) or not error <= TOLERANCE:
            print(f"Received coordinates: {points}")
            success = False
    return success


def run(solver):
    with tempfile.TemporaryDirectory() as directory:
        workdir = Path(directory)
        for name in ("body.inp", "config.yml", "precice-config.xml"):
            shutil.copy(HERE / name, workdir)
        flags = subprocess.check_output(
            ["pkg-config", "--cflags", "--libs", "libprecice"], text=True
        )
        subprocess.run(
            ["c++", "-std=c++17", str(HERE / "probe.cpp"), *shlex.split(flags), "-o", "probe"],
            cwd=workdir,
            check=True,
        )
        env = os.environ | {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
        with (
            (workdir / "solid.log").open("w") as solid_log,
            (workdir / "probe.log").open("w") as probe_log,
        ):
            commands = [
                ([str(solver), "-i", "body", "-precice-participant", "Solid"], solid_log),
                ([str(workdir / "probe")], probe_log),
            ]
            processes = [
                subprocess.Popen(command, cwd=workdir, env=env, stdout=log, stderr=subprocess.STDOUT)
                for command, log in commands
            ]
            try:
                statuses = [process.wait(timeout=30) for process in processes]
            except subprocess.TimeoutExpired:
                statuses = ["timed out"]
            finally:
                for process in processes:
                    if process.poll() is None:
                        process.kill()
                    process.wait()
        if statuses != [0, 0]:
            for name in ("solid.log", "probe.log"):
                print((workdir / name).read_text(), file=sys.stderr)
            raise SystemExit(f"Coupling failed ({statuses})")
        return check_points(workdir)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(f"Usage: {sys.argv[0]} /path/to/ccx_preCICE")
    raise SystemExit(0 if run(Path(sys.argv[1]).resolve(strict=True)) else 1)
