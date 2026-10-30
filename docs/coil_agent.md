# MCP Coil Design Service

The separate `ivsolver-coil-agent` service exposes an MCP tool that searches a constrained three-coil planar flower layout and delegates electromagnetic calculations to the existing `ivsolver` backend.

## Current supported topology

- Three identical coils lie in the XY plane, with their centres spaced 120 degrees apart.
- Each coil is an outward-growing spiral with a shared turn count, inner radius, and turn pitch.
- The optimizer searches inner radius, turn pitch, and lobe-centre radius within caller-supplied bounds.
- The solver computes the inductance matrix and wire lengths. The service derives DC resistance from wire resistivity and cross-section area, then ranks feasible designs by mean self-inductance/resistance.
- A design is feasible only when its self-inductance and all pairwise mutual inductances meet the supplied tolerances. Zero mutual coupling is treated as a tolerance, since numerical geometry and solver discretization do not make exact zero a useful input target.
- The result contains the candidate geometry, inductance matrix, lengths, resistance estimates, a fabrication description, the generated solver wire definition, and a PNG image content block for clients that render tool images.

This is a bounded grid search with local refinement over a symmetric starter topology, capped at 150 solver evaluations per request. It is not a general-purpose shape optimizer, CAD system, or manufacturing sign-off tool. The candidate filter rejects overlapping outer bounding circles, but does not yet perform detailed collision, thermal, or AC-loss analysis.

## Inputs an agent must obtain

The natural-language request can provide target inductance and turn count, but optimization also requires:

- wire cross-section radius in millimetres;
- conductor resistivity in ohm-metres;
- search bounds for inner radius, turn pitch, and lobe-centre radius;
- acceptable self-inductance error and mutual-inductance magnitude;
- optionally, maximum footprint radius.

The MCP server exposes `design_requirements` to provide the agent with this checklist. The agent should ask for missing physical and geometric constraints rather than silently inventing them. If a user authorizes defaults, the returned `constraints` object records the selected values. The sample request's `10 uH` target, three turns, and zero-coupling goal are represented directly; material, wire size, geometry bounds, and numerical tolerances remain required.

The design tool returns the preview as MCP PNG image content. MCP defines image content blocks, but inline display depends on the client. For a guaranteed custom interface in ChatGPT, the remote server would need to be packaged as an MCP App with a UI; the current build is a stdio tool server.

## Modules

- `coil_agent/models.py` validates dimensional design inputs.
- `coil_agent/geometry.py` compiles the supported flower topology to the existing wire-definition language.
- `coil_agent/solver.py` runs the independent `ivsolver` backend and parses its matrix and trajectories.
- `coil_agent/optimizer.py` searches the supplied geometry bounds and ranks feasible candidates by L/R.
- `coil_agent/render.py` creates the SVG preview.
- `coil_agent/description.py` writes a human-readable build description.
- `coil_agent/server.py` exposes the workflow through MCP over stdio.
- `Makefile.agent` builds a separate backend and packages the MCP service as `ivsolver-coil-agent`.

Install Python dependencies with `python3 -m pip install -r requirements-agent.txt`, then run `make -f Makefile.agent`. The original `Makefile` and `build/ivsolver` target remain independent.

The README contains local stdio, Secure MCP Tunnel, ChatGPT, and Responses API setup steps, plus a worked 10 µH design example with its rendered output and generated solver wire definition.

Set `IVSOLVER_BIN` to select a different backend when running the MCP server directly. The packaged service uses its bundled backend by default.
