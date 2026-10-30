# `b`: three-turn cylindrical pickup

## Input fields and line order

| Line(s) | Fields, in parser order | Meaning |
|---|---|---|
| 1 | `n` | One wire. |
| 2 | `cs_radius m` | Cross-section radius in metres and two commands: `0.001 2`. |
| 3 | `s k` | First form a spiral from `k = 4` key points. |
| 4–7 | `x y z` on each of `k` lines | Key-point coordinates in metres. Equal angular positions make three complete turns. |
| 8 | `b radius` | Bend **all points formed so far** around the x-axis with cylinder radius `0.10 m`; this command has no following point lines. |

All coordinates, lengths, and wire radii are in metres. One line feed separates each numbered input line. The `m` value on line 2 counts directive commands, while `k` counts the point lines belonging to one directive.

## Design intent

Create three outward-growing turns at 40, 50, 60, and 70 mm radial key points, then conform the entire spiral to a cylindrical housing. The resulting path is curved in 3D rather than lying in one plane.

For each existing point `(x, y, z)`, `b` keeps `x` and maps the other coordinates to `((R+z) sin(y/R), (R+z) cos(y/R) − R)` in the new y/z plane, with `R = 0.10 m`. This is why the top view still resembles a spiral while the 3D view shows the winding wrapped onto a cylinder.

## Complete wire definition

```text
1
0.001 2
s 4
0.04 0 0
0.05 0 0
0.06 0 0
0.07 0 0
b 0.10
```

The same input is saved as [wire_def](wire_def). After `make`, run:

```bash
build/ivsolver examples/directives/cylindrical_bend_b/wire_def examples/directives/cylindrical_bend_b/result
```

## Generated winding

![Rendered `b`: three-turn cylindrical pickup trajectory](winding.png)

## Calculated results

| Solver output | Value |
|---|---:|
| Self-inductance | 1.52907 µH |
| Wire centreline length | 1.03709 m |

`b` changes existing points; it does not append a new wire segment. It preserves x and maps y and z onto a cylinder around the x-axis. See the [solver matrix and length output](result.lmat) and [sampled trajectory](result_0.trj).
