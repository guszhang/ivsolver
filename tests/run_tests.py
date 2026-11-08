#!/usr/bin/env python3
"""IVSolver regression and physics checks (standard library only).

Usage: python3 tests/run_tests.py [path/to/ivsolver]
"""

import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "build" / "ivsolver").resolve()
MU0 = 4e-7 * math.pi

failures = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}{(' - ' + detail) if detail else ''}")
    if not ok:
        failures.append(name)


def run_json(wire_def, mode="-l", obs_def=None, extra=()):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "wire_def").write_text(wire_def)
        args = [str(BIN), str(tmp / "wire_def"), "-", mode]
        if obs_def is not None:
            (tmp / "obs_def").write_text(obs_def)
            args.append(str(tmp / "obs_def"))
        args += ["--json", *extra]
        done = subprocess.run(args, capture_output=True, text=True, check=False)
        result = json.loads(done.stdout)
        if done.returncode:
            raise RuntimeError(result.get("error", done.stderr))
        return result


def vec(flat, i):
    """Complex 3-vector i from a flattened [xr, xi, yr, yi, zr, zi, ...] array."""
    b = 6 * i
    return [complex(flat[b], flat[b + 1]), complex(flat[b + 2], flat[b + 3]), complex(flat[b + 4], flat[b + 5])]


def points_def(points):
    return "points %d\n" % len(points) + "".join("%.12g %.12g %.12g\n" % p for p in points)


def loop(radius, z=0.0, wire_radius=0.0005, extra=""):
    """A full circle in the xy plane at height z, using the elliptic-arc directive."""
    count = 1 + len(extra.strip().splitlines())
    return f"{wire_radius} {count}\ne 0 0 {z} 0 0 1 1 0 0 {radius} {radius} 0 360\n" + extra


def test_examples():
    for wire_def in sorted((ROOT / "examples" / "directives").glob("*/wire_def")):
        expected_file = wire_def.parent / "result.lmat"
        if not expected_file.exists():
            continue
        rows = [[float(v) for v in line.split()] for line in expected_file.read_text().splitlines() if line.strip()]
        res = run_json(wire_def.read_text())
        n = len(res["inductance"])
        worst = 0.0
        for i in range(n):
            for j in range(n):
                worst = max(worst, abs(res["inductance"][i][j] - rows[i][j]) / abs(rows[i][j]))
            worst = max(worst, abs(res["wires"][i]["length"] - rows[n][i]) / rows[n][i])
        check(f"example {wire_def.parent.name} matches result.lmat", worst < 1e-5, f"max rel diff {worst:.1e}")


def test_loop_axis_field():
    R, I = 0.1, 2.0
    zs = [0.0, 0.02, 0.05, 0.1, 0.2]
    res = run_json("1\n" + loop(R, extra=f"i {I} 0\n"), "-p", points_def([(0, 0, z) for z in zs]), ["--freq", "1000"])
    worst = 0.0
    for k, z in enumerate(zs):
        bz = vec(res["field"]["B"], k)[2].real
        analytic = MU0 * I * R * R / (2 * (R * R + z * z) ** 1.5)
        worst = max(worst, abs(bz - analytic) / analytic)
    check("loop on-axis B matches analytic formula", worst < 1e-3, f"max rel err {worst:.1e}")


def test_faraday():
    """Line integral of E around a coaxial circle equals -j w times the flux of B."""
    R, I, f, r, z = 0.1, 1.0, 50e3, 0.06, 0.03
    nr, nphi = 200, 64
    ring = [(r * math.cos(2 * math.pi * k / nphi), r * math.sin(2 * math.pi * k / nphi), z) for k in range(nphi)]
    radial = [((k + 0.5) * r / nr, 0.0, z) for k in range(nr)]
    res = run_json("1\n" + loop(R, extra=f"i {I} 0\n"), "-p", points_def(ring + radial), ["--freq", str(f)])
    E = res["field"]["E"]
    B = res["field"]["B"]
    circulation = 0j
    for k in range(nphi):
        phi = 2 * math.pi * k / nphi
        e = vec(E, k)
        circulation += (-math.sin(phi) * e[0] + math.cos(phi) * e[1]) * (2 * math.pi * r / nphi)
    flux = sum(vec(B, nphi + k)[2] * 2 * math.pi * radial[k][0] * (r / nr) for k in range(nr))
    expected = -1j * 2 * math.pi * f * flux
    err = abs(circulation - expected) / abs(expected)
    check("Faraday's law: closed integral of E = -jw flux", err < 2e-3, f"rel err {err:.1e}")


def box_faces(xmin, xmax, ymin, ymax, zmin, zmax, n):
    """Midpoint quadrature points (point, outward normal, area) on the faces of a box."""
    out = []

    def face(axis, value, sign, a_range, b_range):
        (a0, a1), (b0, b1) = a_range, b_range
        da, db = (a1 - a0) / n, (b1 - b0) / n
        for i in range(n):
            for j in range(n):
                a, b = a0 + (i + 0.5) * da, b0 + (j + 0.5) * db
                p = [0.0, 0.0, 0.0]
                others = [k for k in range(3) if k != axis]
                p[axis], p[others[0]], p[others[1]] = value, a, b
                normal = [0.0, 0.0, 0.0]
                normal[axis] = sign
                out.append((tuple(p), normal, da * db))

    ranges = [(xmin, xmax), (ymin, ymax), (zmin, zmax)]
    for axis in range(3):
        others = [ranges[k] for k in range(3) if k != axis]
        face(axis, ranges[axis][0], -1.0, *others)
        face(axis, ranges[axis][1], 1.0, *others)
    return out


def test_poynting_theorem():
    """Net Re{S} flowing into a box around a resistive receiver equals its I^2 R / 2 loss."""
    R, gap, f, R2 = 0.1, 0.1, 100e3, 0.5
    wire_def = "2\n" + loop(R, 0.0, extra="i 1 0\n") + loop(R, gap, extra=f"R {R2}\nT 0\n")
    quad = box_faces(-0.16, 0.16, -0.16, 0.16, gap - 0.04, gap + 0.04, 40)
    res = run_json(wire_def, "-p", points_def([q[0] for q in quad]), ["--freq", str(f)])
    S = res["field"]["S"]
    outflow = 0.0
    for k, (_, normal, area) in enumerate(quad):
        s = vec(S, k)
        outflow += sum(s[c].real * normal[c] for c in range(3)) * area
    i2 = complex(*res["circuit"]["current"][1])
    loss = 0.5 * abs(i2) ** 2 * R2
    reported = res["circuit"]["resistive_loss"][1]
    err = abs(-outflow - loss) / loss
    check("Poynting theorem around the receiver", err < 1e-2, f"inflow {-outflow:.4e} W vs loss {loss:.4e} W, rel err {err:.1e}")
    check("circuit resistive loss reported consistently", abs(reported - loss) / loss < 1e-9)

    # Resonant receiver: I2 = -j w M I1 / R2 exactly.
    M = res["inductance"][0][1]
    expected = -1j * 2 * math.pi * f * M / R2
    err = abs(i2 - expected) / abs(expected)
    check("tuned receiver current equals -jwM I1 / R2", err < 1e-9, f"rel err {err:.1e}")

    # Energy conservation in the circuit: source power = total resistive loss.
    p_in = sum(p[0] for p in res["circuit"]["source_power"])
    p_r = sum(res["circuit"]["resistive_loss"])
    check("circuit power balance", abs(p_in - p_r) <= 1e-9 * max(p_r, 1e-30), f"in {p_in:.4e} W, loss {p_r:.4e} W")


def test_threads_deterministic():
    wire_def = "2\n" + loop(0.1, 0.0, extra="i 1 0\n") + loop(0.08, 0.05, extra="i 0 1\n")
    obs = "grid -0.1 0.1 5 -0.1 0.1 5 -0.05 0.1 4\n"
    a = run_json(wire_def, "-p", obs, ["--freq", "1e5", "--circuit", "-t", "1"])
    b = run_json(wire_def, "-p", obs, ["--freq", "1e5", "--circuit", "-t", "7"])
    same = a["inductance"] == b["inductance"] and a["field"]["S"] == b["field"]["S"]
    check("results independent of thread count", same)


def test_kernels_agree_far_field():
    wire_def = "1\n" + loop(0.1, extra="i 1 0\n")
    obs = points_def([(0.3, 0.1, 0.2), (0.0, 0.25, -0.1)])
    ex = run_json(wire_def, "-p", obs, ["--freq", "1e4", "--kernel", "exact"])
    mp = run_json(wire_def, "-p", obs, ["--freq", "1e4", "--kernel", "midpoint"])
    worst = 0.0
    for q in ("E", "B"):
        for k in range(2):
            a, b = vec(ex["field"][q], k), vec(mp["field"][q], k)
            na = math.sqrt(sum(abs(c) ** 2 for c in a))
            worst = max(worst, math.sqrt(sum(abs(x - y) ** 2 for x, y in zip(a, b))) / na)
    check("exact and midpoint kernels agree away from the wire", worst < 1e-5, f"max rel diff {worst:.1e}")


def test_near_wire_finite():
    wire_def = "1\n" + loop(0.1, extra="i 1 0\n")
    res = run_json(wire_def, "-p", points_def([(0.1, 0.0, 0.0), (0.1002, 0.0, 0.0)]), ["--freq", "1e4"])
    values = res["field"]["B"] + res["field"]["E"]
    finite = all(v is not None and math.isfinite(v) for v in values)
    check("fields stay finite on and inside the conductor", finite)


def test_errors():
    for text, needle in (("1\n0.001 1\nq 1\n", "unknown directive"), ("1\n0.001 2\nl 2\n0 0 0\n", "unexpected end")):
        try:
            run_json(text)
            check(f"error reported for '{needle}'", False)
        except RuntimeError as error:
            check(f"error reported for '{needle}'", needle in str(error), str(error))


def test_legacy_cli():
    """The original command forms and output files keep working unchanged."""
    example = ROOT / "examples" / "python" / "b-field_plot"
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        wire_def = str(example / "wire_def")

        # ivsolver wire_def output            (inductance mode, 2 arguments)
        done = subprocess.run([str(BIN), wire_def, str(tmp / "a")], capture_output=True, text=True)
        files = sorted(p.name for p in tmp.iterdir())
        check("legacy 'ivsolver in out' writes .lmat and .trj", done.returncode == 0 and files == ["a.lmat", "a_0.trj"], str(files))
        lmat = (tmp / "a.lmat").read_bytes()
        check("legacy .lmat layout (n rows + lengths row, CRLF, tabs)", lmat.count(b"\r\n") == 1 and b"\t" in lmat and len(lmat.split(b"\r\n")[1].split()) == 1)

        # ivsolver wire_def output -l dummy    (the old 4-argument form)
        done = subprocess.run([str(BIN), wire_def, str(tmp / "b"), "-l", "unused"], capture_output=True, text=True)
        check("legacy 'ivsolver in out -l x' accepted", done.returncode == 0 and (tmp / "b.lmat").read_bytes() == lmat)

        # ivsolver wire_def output -b obs_def  (B-field mode)
        done = subprocess.run([str(BIN), wire_def, str(tmp / "c"), "-b", str(example / "obs_def")], capture_output=True, text=True)
        rows = [[float(v) for v in line.split()] for line in (tmp / "c.obs").read_text().splitlines()]
        ref = [[float(v) for v in line.split()] for line in (example / "wires.obs").read_text().splitlines()]
        same_shape = len(rows) == len(ref) and all(len(r) == 9 for r in rows)
        worst = max(abs(a - b) / max(abs(b), 1e-12) for r1, r2 in zip(rows, ref) for a, b in zip(r1, r2)) if same_shape else 1
        check("legacy -b output matches committed example wires.obs", done.returncode == 0 and same_shape and worst < 1e-5, f"{len(rows)} rows, max rel diff {worst:.1e}")
        check("legacy -b writes trajectories", (tmp / "c_0.trj").exists())


def test_cores():
    """Linear magnetic cores (volume integral method)."""
    R, wr, a = 1.0, 0.002, 0.02
    base = f"1\n{wr} 2\ne 0 0 0 0 0 1 1 0 0 {R} {R} 0 360\ni 1 0\n"
    L0 = run_json(base)["inductance"][0][0]
    B0 = MU0 / (2 * R)

    # mu_r = 1 is air: nothing changes.
    same = run_json(base + "cores 1\nbox 1 0.01 -0.01 -0.01 -0.01 0.01 0.01 0.01\n")["inductance"][0][0]
    check("core with mu_r = 1 leaves L unchanged", abs(same - L0) < 1e-18)

    # One cubic cell is the uniformly magnetised cube, whose demagnetising factor is exactly 1/3.
    for mur in (2, 100, 1e4):
        chi = mur - 1
        j = run_json(base + f"cores 1\nbox {mur} {a} {-a/2} {-a/2} {-a/2} {a/2} {a/2} {a/2}\n")
        expected = a ** 3 * chi / (1 + chi / 3) * B0 ** 2 / MU0
        got = j["inductance"][0][0] - L0
        check(f"single-cell cube dL, mu_r={mur:g}", abs(got / expected - 1) < 1e-3, f"ratio {got / expected:.5f}")

    # Small sphere at the centre of a big loop (nearly uniform field): dL = V 3(mu_r-1)/(mu_r+2) B0^2/mu0.
    r, h, mur = 0.01, 0.002, 2
    j = run_json(base + f"cores 1\nsphere {mur} {h} 0 0 0 {r}\n")
    V = j["core"]["cells"] * h ** 3
    expected = V * 3 * (mur - 1) / (mur + 2) * B0 ** 2 / MU0
    got = j["inductance"][0][0] - L0
    check("sphere in a loop, mu_r=2 (dL)", abs(got / expected - 1) < 0.02, f"ratio {got / expected:.4f}")
    check("core solve converged", j["core"]["residual"] < 1e-6, f"residual {j['core']['residual']:.1e}")

    # Inside the sphere B = 3 mu_r/(mu_r+2) B0; far away the core hardly matters.
    pts = [(0.0011, 0.0013, 0.0017), (0, 0, 0.5)]
    air = run_json(base, "-p", points_def(pts), ["--freq", "1000"])
    core = run_json(base + f"cores 1\nsphere {mur} {h} 0 0 0 {r}\n", "-p", points_def(pts), ["--freq", "1000"])
    ratio_in = core["field"]["B"][4] / air["field"]["B"][4]
    check("B inside a mu_r=2 sphere", abs(ratio_in / 1.5 - 1) < 0.03, f"ratio {ratio_in:.4f}")
    ratio_far = core["field"]["B"][10] / air["field"]["B"][10]
    check("B far from the core is nearly unchanged", abs(ratio_far - 1) < 1e-3, f"ratio {ratio_far:.5f}")

    # Two coaxial loops around a rod: L stays symmetric, mutual inductance increases.
    two = "2\n" + loop(0.05, 0.0, extra="i 1 0\n") + loop(0.05, 0.05, extra="i 0 0\n")
    bare = run_json(two)["inductance"]
    cored = run_json(two + "cores 1\ncylinder 200 0.01 0 0 0.025 0.01 0.2 z\n")["inductance"]
    check("core L symmetric", abs(cored[0][1] - cored[1][0]) < 1e-15 * max(1, abs(cored[0][1])) + 1e-18)
    check("core raises self and mutual inductance", cored[0][0] > bare[0][0] and cored[0][1] > bare[0][1])

    for bad, what in (("cores 1\nbox 0.5 0.01 0 0 0 1 1 1\n", "mu_r below 1"),
                      ("cores 1\nbox 2 0.01 0 0 0 1\n", "truncated box"),
                      ("cores 1\ncone 2 0.01 0 0 0 1\n", "unknown shape"),
                      ("cores 1\nbox 2 0.0001 0 0 0 1 1 1\n", "too many cells")):
        try:
            run_json(base + bad)
            check(f"core input error: {what}", False)
        except RuntimeError:
            check(f"core input error: {what}", True)


if __name__ == "__main__":
    if not BIN.is_file():
        sys.exit(f"solver binary not found: {BIN} (run make first)")
    test_examples()
    test_loop_axis_field()
    test_faraday()
    test_poynting_theorem()
    test_threads_deterministic()
    test_kernels_agree_far_field()
    test_near_wire_finite()
    test_errors()
    test_legacy_cli()
    test_cores()
    print(f"\n{'All tests passed' if not failures else f'{len(failures)} test(s) failed'}")
    sys.exit(1 if failures else 0)
