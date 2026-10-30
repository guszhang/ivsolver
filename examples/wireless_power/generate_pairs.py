"""Rebuild the hexagon and bowl examples with the plain ivsolver binary."""

from pathlib import Path
import math
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
SOLVER = ROOT / "build" / "ivsolver"
OUT = Path(__file__).resolve().parent


def write_wires(path, wires, radius=0.001):
    lines = [str(len(wires))]
    for points in wires:
        lines.extend((f"{radius} 1", f"l {len(points)}"))
        lines.extend(f"{x:.9g} {y:.9g} {z:.9g}" for x, y, z in points)
    path.write_text("\n".join(lines) + "\n")


def solve(folder, wires, name="result"):
    definition = folder / ("wire_def" if name == "result" else f"{name}.def")
    write_wires(definition, wires)
    prefix = folder / name
    subprocess.run([str(SOLVER), str(definition), str(prefix)], check=True, capture_output=True, text=True)
    matrix_file = prefix.with_suffix(".lmat")
    matrix_file.write_text("\n".join(line.rstrip() for line in matrix_file.read_text().splitlines()) + "\n")
    rows = [[float(value) for value in line.split()] for line in matrix_file.read_text().splitlines()]
    return rows


def hex_spiral(z, dx=0):
    points = []
    for i in range(6 * 3 + 1):
        angle = 2 * math.pi * i / 6
        radius = 0.075 + 0.005 * i / 6
        points.append((dx + radius * math.cos(angle), radius * math.sin(angle), z))
    return points


def circular_spiral(z_fn, dx=0, r0=0.07, pitch=0.02, turns=4):
    count = 240
    points = []
    for i in range(count + 1):
        theta = turns * 2 * math.pi * i / count
        radius = r0 + pitch * theta / (2 * math.pi)
        points.append((dx + radius * math.cos(theta), radius * math.sin(theta), z_fn(radius)))
    return points


def draw(folder, wires, title):
    fig = plt.figure(figsize=(9, 4.5), layout="constrained")
    ax = fig.add_subplot(121)
    ax3 = fig.add_subplot(122, projection="3d")
    for points, color, label in zip(wires, ("#2563eb", "#dc2626"), ("transmitter", "receiver")):
        x, y, z = zip(*points)
        ax.plot(x, y, color=color, linewidth=1.7, label=label)
        ax3.plot(x, y, z, color=color, linewidth=1.5)
    ax.set(aspect="equal", xlabel="x (m)", ylabel="y (m)", title="Top view")
    ax.legend(fontsize=8)
    ax3.set(xlabel="x (m)", ylabel="y (m)", zlabel="z (m)", title="3D view")
    fig.suptitle(title)
    fig.savefig(folder / "winding.png", dpi=160)
    plt.close(fig)


def main():
    if not SOLVER.is_file():
        raise SystemExit("Build the plain solver first: make")
    hex_dir = OUT / "hexagon_pair"
    hex_dir.mkdir(exist_ok=True)
    hex_wires = [hex_spiral(0), hex_spiral(0.04, 0.03)]
    hex_result = solve(hex_dir, hex_wires)
    draw(hex_dir, hex_wires, "Hexagonal winding pair, 30 mm offset and 40 mm gap")
    print("hexagon", hex_result)

    bowl_dir = OUT / "bowl_misalignment"
    bowl_dir.mkdir(exist_ok=True)
    offsets = (0, 0.02, 0.04, 0.06)
    bowl_tx = circular_spiral(lambda r: 3.125 * (r - 0.07) ** 2)
    flat_tx = circular_spiral(lambda r: 0)
    comparisons = []
    for offset in offsets:
        receiver = circular_spiral(lambda r: 0.06, offset, 0.06, 0.015, 3)
        bowl_rows = solve(bowl_dir, [bowl_tx, receiver], f"bowl_{int(offset * 1000)}mm")
        flat_rows = solve(bowl_dir, [flat_tx, receiver], f"flat_{int(offset * 1000)}mm")
        comparisons.append((offset, bowl_rows[0][1], flat_rows[0][1]))
        for variant in ("bowl", "flat"):
            name = f"{variant}_{int(offset * 1000)}mm"
            for suffix in (".def", ".lmat", "_0.trj", "_1.trj"):
                (bowl_dir / f"{name}{suffix}").unlink()
    receiver = circular_spiral(lambda r: 0.06, 0.04, 0.06, 0.015, 3)
    solve(bowl_dir, [bowl_tx, receiver])
    draw(bowl_dir, [bowl_tx, receiver], "Bowl transmitter and receiver, 40 mm offset")
    lines = ["offset_m,bowl_mutual_h,flat_mutual_h"]
    lines += [f"{offset},{bowl},{flat}" for offset, bowl, flat in comparisons]
    (bowl_dir / "misalignment.csv").write_text("\n".join(lines) + "\n")
    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    ax.plot([x[0] * 1000 for x in comparisons], [abs(x[1]) * 1e6 for x in comparisons], "o-", label="Bowl transmitter")
    ax.plot([x[0] * 1000 for x in comparisons], [abs(x[2]) * 1e6 for x in comparisons], "s-", label="Flat reference")
    ax.set(xlabel="Receiver lateral offset (mm)", ylabel="|Mutual inductance| (µH)", title="Misalignment sweep, 60 mm receiver height")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.savefig(bowl_dir / "misalignment.png", dpi=160)
    plt.close(fig)
    print("bowl sweep", comparisons)


if __name__ == "__main__":
    main()
