# `s`: origin-centred spiral

## Input fields and line order

| Line(s) | Fields, in parser order | Meaning |
|---|---|---|
| 1 | `n` | One wire. |
| 2 | `cs_radius m` | Cross-section radius in metres, then one command: `0.001 1`. |
| 3 | `s k` | Spiral command with `k = 3` key-point lines. |
| 4–6 | `x y z` on each of `k` lines | Coordinates in metres. Interpolation rotates around the global z-axis; equal angles at successive points produce a full turn. |

All coordinates, lengths, and wire radii are in metres. One line feed separates each numbered input line. The `m` value on line 2 counts directive commands, while `k` counts the point lines belonging to one directive.

## Design intent

Make two outward-growing circular turns centred on the origin, with key-point radii of 60, 70, and 80 mm. Use this directive for a simple planar spiral without relocating its axis.

## Complete wire definition

```text
1
0.001 1
s 3
0.06 0 0
0.07 0 0
0.08 0 0
```

The same input is saved as [wire_def](wire_def). After `make`, run:

```bash
build/ivsolver examples/directives/origin_spiral_s/wire_def examples/directives/origin_spiral_s/result
```

## Generated winding

![Rendered `s`: origin-centred spiral trajectory](winding.png)

## Calculated results

| Solver output | Value |
|---|---:|
| Self-inductance | 1.15898 µH |
| Wire centreline length | 0.87984 m |

The two intervals between three same-angle points each generate one full turn. See the [solver matrix and length output](result.lmat) and [sampled trajectory](result_0.trj).
