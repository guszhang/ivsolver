# `p`: back-to-back parabolic loop

## Input fields and line order

| Line(s) | Fields, in parser order | Meaning |
|---|---|---|
| 1 | `n` | One wire. |
| 2 | `cs_radius m` | Cross-section radius in metres and two commands: `0.001 2`. |
| 3 | `p ox oy oz nx ny nz xx xy xz a b c x1 x2` | First parabolic edge; all fields are on this one line. |
| 4 | `p ox oy oz nx ny nz xx xy xz a b c x1 x2` | Second edge; the reversed x interval returns to the first endpoint. There are no following key-point lines. |
| 3 and 4, fields 2–4 | `ox oy oz` | Local plane origin coordinates in metres. |
| 3 and 4, fields 5–7 | `nx ny nz` | Plane normal vector. |
| 3 and 4, fields 8–10 | `xx xy xz` | Positive local x-direction vector; local y is normal × x-direction. |
| 3 and 4, fields 11–13 | `a b c` | Coefficients of local `y = ax² + bx + c`; units are 1/m, unitless, and m, respectively. |
| 3 and 4, fields 14–15 | `x1 x2` | Start and end local x coordinates in metres, traversed in that order. |

All coordinates, lengths, and wire radii are in metres. One line feed separates each numbered input line. The `m` value on line 2 counts directive commands, while `k` counts the point lines belonging to one directive.

## Design intent

Make a closed lens-shaped winding from two parabolic edges in the XY plane. The lower edge is `y = 3x²` from −100 to 100 mm. The upper edge is `y = −3x² + 0.06` traversed back from 100 to −100 mm. The edges meet at both ends.

## Complete wire definition

```text
1
0.001 2
p 0 0 0 0 0 1 1 0 0 3 0 0 -0.10 0.10
p 0 0 0 0 0 1 1 0 0 -3 0 0.06 0.10 -0.10
```

The same input is saved as [wire_def](wire_def). After `make`, run:

```bash
build/ivsolver examples/directives/parabolic_curve_p/wire_def examples/directives/parabolic_curve_p/result
```

## Generated winding

![Rendered `p`: back-to-back parabolic loop trajectory](winding.png)

## Calculated results

| Solver output | Value |
|---|---:|
| Self-inductance | 0.326424 µH |
| Wire centreline length | 0.422845 m |

Here `a` has units 1/m, `b` is dimensionless, and `c` is metres. The two directives are consecutive segments of one wire. See the [solver matrix and length output](result.lmat) and [sampled trajectory](result_0.trj).
