# `g`: catenary with a chosen shape parameter

## Input fields and line order

| Line(s) | Fields, in parser order | Meaning |
|---|---|---|
| 1 | `n` | One wire. |
| 2 | `cs_radius m` | Cross-section radius in metres and one command: `0.001 1`. |
| 3 | `g k a` | `k = 2` anchor-point lines; catenary shape parameter `a = 0.08 m`. |
| 4–5 | `x y z` on each of `k` lines | Anchor coordinates in metres. The solver forms a sagging segment between them. |

All coordinates, lengths, and wire radii are in metres. One line feed separates each numbered input line. The `m` value on line 2 counts directive commands, while `k` counts the point lines belonging to one directive.

## Design intent

Model a flexible lead across a 160 mm horizontal span with both supports at z = 0. Choose the catenary parameter directly and let the solver form the hanging profile.

## Complete wire definition

```text
1
0.001 1
g 2 0.08
-0.08 0 0
0.08 0 0
```

The same input is saved as [wire_def](wire_def). After `make`, run:

```bash
build/ivsolver examples/directives/parameter_sag_g/wire_def examples/directives/parameter_sag_g/result
```

## Generated winding

![Rendered `g`: catenary with a chosen shape parameter trajectory](winding.png)

## Calculated results

| Solver output | Value |
|---|---:|
| Self-inductance | 0.185438 µH |
| Wire centreline length | 0.188029 m |

The calculated centreline length exceeds the 160 mm horizontal span because the lead sags. See the [solver matrix and length output](result.lmat) and [sampled trajectory](result_0.trj).
