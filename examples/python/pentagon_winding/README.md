# Parameterised pentagon spiral winding

This example builds a planar spiral winding whose turns follow a rounded pentagon rather
than a circle, and calculates its self-inductance and total wire length with the
`ivsolver` binary.

The base shape is an inner pentagon whose five vertices act as the centre poles for the
cornering arcs. The winding is parameterised by

- `r_pentagon` -- circumradius of the inner pentagon,
- `r_cs` -- wire cross-section radius,
- `pitch` -- radial advance per turn,
- `n_turns` -- number of turns,
- `r_corner_0` -- corner arc radius of the innermost turn.

The wire walks outwards turn by turn, alternating the `a` and `l` directives: an arc
about a corner pole, then a straight run along the offset edge to the next pole. Each
corner advances the radius by `pitch / 5`, so one lap of five corners advances exactly one
`pitch`.

This is the use case the `a` directive exists for. The `s` directive can only spiral about
the origin, so a winding with five separate corner poles would otherwise need the whole
pattern translated between corners.

The notebook cross-checks the solver's reported length against a closed-form expression
for the same geometry; the two agree to within the discretisation error of approximating
the arcs by chords of one wire diameter (about 0.01% for the default parameters).

![Pentagon winding](winding.png)

## Usage

Build the solver first (`make` in the repository root), then run `pentagon.ipynb`.
