# Wireless power winding studies

These examples use three different winding geometries. The [flower pad](compact_pad/README.md) uses the current MCP design tool. The [hexagonal pair](hexagon_pair/README.md) and [bowl pair](bowl_misalignment/README.md) use the base `ivsolver` wire language directly. The MCP tool does not yet design arbitrary hexagons or 3D bowls; those two studies show how the solver can support future design modules and an LLM workflow that generates and evaluates wire definitions.

| Study | Geometry | Main question |
|---|---|---|
| [Flower pad](compact_pad/README.md) | Three planar spirals | How can three pad windings meet low mutual-coupling limits? |
| [Hexagonal pair](hexagon_pair/README.md) | Two facing polygon spirals | How do a finite air gap and lateral offset affect a hexagonal transmitter/receiver pair? |
| [Bowl and receiver](bowl_misalignment/README.md) | 3D bowl transmitter, flat receiver | Does the bowl retain coupling better as the receiver shifts sideways? |

Build the plain solver with `make`. To reproduce the last two studies, run `python3 examples/wireless_power/generate_pairs.py` from the repository root; it needs Matplotlib for the figures. The script writes each `wire_def`, runs `build/ivsolver`, and renders the calculated trajectories. Build and connect the MCP binary for the flower pad using the [main README](../../README.md#mcp-coil-design-service).

These are air-core inductance studies. A complete wireless power design needs a receiver placement range, operating frequency, compensation network, AC losses, ferrite/shield models where relevant, and electrical constraints.
