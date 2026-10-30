# `l`: three-turn square spiral

## Input fields and line order

| Line(s) | Fields, in parser order | Meaning |
|---|---|---|
| 1 | `n` | Number of wires; `1` here. |
| 2 | `cs_radius m` | Wire cross-section radius in metres and number of commands; `0.001 1` here. |
| 3 | `l k` | Straight-line command and count of following key-point lines; `k = 13`. |
| 4–16 | `x y z` on each of `k` lines | Key points in metres, in travel order. Each consecutive pair is joined by a straight segment. |

All coordinates, lengths, and wire radii are in metres. One line feed separates each numbered input line. The `m` value on line 2 counts directive commands, while `k` counts the point lines belonging to one directive.

## Design intent

Make a three-turn polygonal pickup with expanding square sides. The inner, middle, and outer turns use corner extents of 40, 50, and 60 mm. Diagonal transitions connect one turn to the next. The path ends at the outer starting corner; it is an open spiral with terminals at different radii.

## Complete wire definition

```text
1
0.001 1
l 13
-0.04 -0.04 0
0.04 -0.04 0
0.04 0.04 0
-0.04 0.04 0
-0.05 -0.05 0
0.05 -0.05 0
0.05 0.05 0
-0.05 0.05 0
-0.06 -0.06 0
0.06 -0.06 0
0.06 0.06 0
-0.06 0.06 0
-0.06 -0.06 0
```

The same input is saved as [wire_def](wire_def). After `make`, run:

```bash
build/ivsolver examples/directives/linear_polygon_l/wire_def examples/directives/linear_polygon_l/result
```

## Generated winding

![Rendered `l`: three-turn square spiral trajectory](winding.png)

## Calculated results

| Solver output | Value |
|---|---:|
| Self-inductance | 1.74479 µH |
| Wire centreline length | 1.22101 m |

The result is a single-wire self-inductance. The solver connects all 13 key points in order. See the [solver matrix and length output](result.lmat) and [sampled trajectory](result_0.trj).
