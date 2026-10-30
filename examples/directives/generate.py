"""Generate one runnable IVSolver case for every wire-forming directive."""

from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
SOLVER = ROOT / "build" / "ivsolver"
OUT = Path(__file__).resolve().parent

CASES = {
    "linear_polygon_l": ("Three-turn square spiral", [
        "l 13", "-0.04 -0.04 0", "0.04 -0.04 0", "0.04 0.04 0", "-0.04 0.04 0",
        "-0.05 -0.05 0", "0.05 -0.05 0", "0.05 0.05 0", "-0.05 0.05 0",
        "-0.06 -0.06 0", "0.06 -0.06 0", "0.06 0.06 0", "-0.06 0.06 0", "-0.06 -0.06 0",
    ]),
    "origin_spiral_s": ("Origin-centred planar spiral", [
        "s 3", "0.06 0 0", "0.07 0 0", "0.08 0 0",
    ]),
    "offset_arc_a": ("Offset-axis rounded coil", [
        "a 0.04 0 1 3", "0.10 0 0", "0.11 0 0", "0.12 0 0",
    ]),
    "parameter_sag_g": ("Suspended lead with selected sag parameter", [
        "g 2 0.08", "-0.08 0 0", "0.08 0 0",
    ]),
    "length_sag_c": ("Three-sided sagging triangle", [
        "c 4", "0 1.154700538 0 0", "-1 -0.577350269 0 2.2",
        "1 -0.577350269 0 2.2", "0 1.154700538 0 2.2",
    ]),
    "elliptic_loop_e": ("Partial elliptical pickup", [
        "e 0 0 0 0 0 1 1 0 0 0.10 0.06 10 350",
    ]),
    "parabolic_curve_p": ("Back-to-back parabolic loop", [
        "p 0 0 0 0 0 1 1 0 0 3 0 0 -0.10 0.10",
        "p 0 0 0 0 0 1 1 0 0 -3 0 0.06 0.10 -0.10",
    ]),
    "cylindrical_bend_b": ("Three-turn cylindrical pickup", [
        "s 4", "0.04 0 0", "0.05 0 0", "0.06 0 0", "0.07 0 0",
        "b 0.10",
    ]),
}


def main():
    if not SOLVER.is_file():
        raise SystemExit("Build ivsolver first with make")
    for name, (title, commands) in CASES.items():
        folder = OUT / name
        folder.mkdir(exist_ok=True)
        directive_count = sum(line[0] in "lsagcepb" and " " in line for line in commands)
        (folder / "wire_def").write_text(f"1\n0.001 {directive_count}\n" + "\n".join(commands) + "\n")
        prefix = folder / "result"
        subprocess.run([str(SOLVER), str(folder / "wire_def"), str(prefix)], check=True, capture_output=True, text=True)
        matrix_file = folder / "result.lmat"
        matrix_file.write_text("\n".join(line.rstrip() for line in matrix_file.read_text().splitlines()) + "\n")
        points = [[float(value) for value in line.split()] for line in (folder / "result_0.trj").read_text().splitlines()]
        x, y, z = zip(*points)
        fig = plt.figure(figsize=(8, 4), layout="constrained")
        ax = fig.add_subplot(121)
        ax3 = fig.add_subplot(122, projection="3d")
        ax.plot(x, y, color="#2563eb", linewidth=2)
        ax.scatter([x[0], x[-1]], [y[0], y[-1]], color=["#16a34a", "#dc2626"], s=25)
        ax.set(aspect="equal", xlabel="x (m)", ylabel="y (m)", title="Top view")
        ax3.plot(x, y, z, color="#2563eb", linewidth=2)
        ax3.set(xlabel="x (m)", ylabel="y (m)", zlabel="z (m)", title="3D view")
        fig.suptitle(title)
        fig.savefig(folder / "winding.png", dpi=150)
        plt.close(fig)
        print(name, (folder / "result.lmat").read_text().strip().replace("\n", " | "))


if __name__ == "__main__":
    main()
