# `e`: partial elliptical pickup

## Input fields and line order

| Line(s) | Fields, in parser order | Meaning |
|---|---|---|
| 1 | `n` | One wire. |
| 2 | `cs_radius m` | Cross-section radius in metres and one command: `0.001 1`. |
| 3 (one line) | `e cx cy cz nx ny nz xx xy xz major minor start end` | All fields on the same line, in this exact order; no following key-point lines. |
| 3, fields 2–4 | `cx cy cz` | Ellipse centre coordinates in metres. |
| 3, fields 5–7 | `nx ny nz` | Plane normal vector. |
| 3, fields 8–10 | `xx xy xz` | Positive local x-direction vector, projected into the plane; local y is normal × x-direction. |
| 3, fields 11–12 | `major minor` | Major and minor **semi-radii** in metres. |
| 3, fields 13–14 | `start end` | Start and end angles in degrees, sampled in that order. |

All coordinates, lengths, and wire radii are in metres. One line feed separates each numbered input line. The `m` value on line 2 counts directive commands, while `k` counts the point lines belonging to one directive.

## Design intent

Create a partial XY-plane elliptical pickup with 100 mm and 60 mm semi-radii. The curve starts at 10° and ends at 350°, leaving a 20° opening for terminals.

## Complete wire definition

```text
1
0.001 1
e 0 0 0 0 0 1 1 0 0 0.10 0.06 10 350
```

The same input is saved as [wire_def](wire_def). After `make`, run:

```bash
build/ivsolver examples/directives/elliptic_loop_e/wire_def examples/directives/elliptic_loop_e/result
```

## Generated winding

![Rendered `e`: partial elliptical pickup trajectory](winding.png)

## Calculated results

| Solver output | Value |
|---|---:|
| Self-inductance | 0.454968 µH |
| Wire centreline length | 0.489402 m |

The x-direction is projected into the plane; the positive local y-direction is normal × x-direction. The path is intentionally open. See the [solver matrix and length output](result.lmat) and [sampled trajectory](result_0.trj).
