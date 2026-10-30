# `c`: three sagging sides of a triangular loop

## Input fields and line order

| Line(s) | Fields, in parser order | Meaning |
|---|---|---|
| 1 | `n` | One wire. |
| 2 | `cs_radius m` | Cross-section radius in metres and one command: `0.001 1`. |
| 3 | `c k` | Length-specified catenary command; `k = 4` point lines. |
| 4 | `x y z length` | First vertex in metres; its length value is unused because no segment enters this point. |
| 5–7 | `x y z length` on each following line | Next vertex and the cable length from the **previous** vertex to this one. The last vertex repeats the first to close the triangle. |

All coordinates, lengths, and wire radii are in metres. One line feed separates each numbered input line. The `m` value on line 2 counts directive commands, while `k` counts the point lines belonging to one directive.

## Design intent

Suspend a triangular winding from three XY-plane anchor points. Each horizontal side spans 2 m while its cable length is 2.2 m, so all three sides sag below the anchor plane. The repeated first point closes the loop.

## Complete wire definition

```text
1
0.001 1
c 4
0 1.154700538 0 0
-1 -0.577350269 0 2.2
1 -0.577350269 0 2.2
0 1.154700538 0 2.2
```

The same input is saved as [wire_def](wire_def). After `make`, run:

```bash
build/ivsolver examples/directives/length_sag_c/wire_def examples/directives/length_sag_c/result
```

## Generated winding

![Rendered `c`: three sagging sides of a triangular loop trajectory](winding.png)

## Calculated results

| Solver output | Value |
|---|---:|
| Self-inductance | 8.34463 µH |
| Wire centreline length | 6.6 m |

The three 2.2 m lengths sum to the calculated 6.6 m wire length. The first `0` length field is a placeholder. See the [solver matrix and length output](result.lmat) and [sampled trajectory](result_0.trj).
