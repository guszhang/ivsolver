# Wire-forming directive gallery

Each case is a standalone one-wire `ivsolver` input with a rendered top and 3D view, a trajectory, and a `.lmat` result. All coordinates and wire radii are in metres. The gallery covers every path-forming directive (`l`, `s`, `a`, `g`, `c`, `e`, `p`) and the surface-bending directive (`b`). Placement transforms (`x`, `y`, `z`, `o`, `n`) and current setting (`i`) are outside this winding-shape gallery.

| Directive | Use case | Calculated self-inductance |
|---|---|---:|
| [`l`](linear_polygon_l/README.md) | Three-turn square spiral | 1.7448 µH |
| [`s`](origin_spiral_s/README.md) | Planar circular spiral about the origin | 1.1590 µH |
| [`a`](offset_arc_a/README.md) | Spiral around a displaced winding centre | 1.1590 µH |
| [`g`](parameter_sag_g/README.md) | Suspended lead with prescribed catenary parameter | 0.1854 µH |
| [`c`](length_sag_c/README.md) | Three sagging sides of a triangular loop | 8.3446 µH |
| [`e`](elliptic_loop_e/README.md) | Partial elliptical pickup loop | 0.4550 µH |
| [`p`](parabolic_curve_p/README.md) | Back-to-back parabolic loop | 0.3264 µH |
| [`b`](cylindrical_bend_b/README.md) | Three-turn cylindrical pickup | 1.5291 µH |

Build the plain solver with `make`, then regenerate all inputs, solver results, and PNGs with `python3 examples/directives/generate.py` (Matplotlib is needed for rendering). Each example also gives a direct solver command. The figures show centreline geometry; they do not represent conductor thickness. The single-wire self-inductance values are useful geometry checks, not full wireless-power-system performance estimates.
