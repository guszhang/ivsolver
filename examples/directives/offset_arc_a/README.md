# `a`: spiral around a displaced axis

## Input fields and line order

| Line(s) | Fields, in parser order | Meaning |
|---|---|---|
| 1 | `n` | One wire. |
| 2 | `cs_radius m` | Cross-section radius in metres and one command: `0.001 1`. |
| 3 | `a axis_x axis_y direction k` | Vertical axis position in metres; `direction = 1` means increasing angle; `k = 3` key-point lines. |
| 4–6 | `x y z` on each of `k` lines | Coordinates in metres. The arc between points is measured around the stated axis. |

All coordinates, lengths, and wire radii are in metres. One line feed separates each numbered input line. The `m` value on line 2 counts directive commands, while `k` counts the point lines belonging to one directive.

## Design intent

Make a two-turn spiral centred 40 mm to the right of the global origin. The key points at x = 100, 110, and 120 mm correspond to local radii of 60, 70, and 80 mm. This places a coil within a larger assembly without translating the whole wire after formation.

## Complete wire definition

```text
1
0.001 1
a 0.04 0 1 3
0.10 0 0
0.11 0 0
0.12 0 0
```

The same input is saved as [wire_def](wire_def). After `make`, run:

```bash
build/ivsolver examples/directives/offset_arc_a/wire_def examples/directives/offset_arc_a/result
```

## Generated winding

![Rendered `a`: spiral around a displaced axis trajectory](winding.png)

## Calculated results

| Solver output | Value |
|---|---:|
| Self-inductance | 1.15898 µH |
| Wire centreline length | 0.87984 m |

Its self-inductance matches the `s` example because translating a winding does not change its self-inductance. See the [solver matrix and length output](result.lmat) and [sampled trajectory](result_0.trj).
