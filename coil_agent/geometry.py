"""Compile the first supported high-level layout into IVSolver directives."""

import math

from .models import Candidate, DesignRequest


def flower_wire_definition(request: DesignRequest, candidate: Candidate) -> str:
    """Build three equal outward-growing planar spirals at 120 degree intervals.

    In this first topology, each spiral is one lobe of the three-lobe assembly.
    """
    lines = ["3"]
    for index in range(3):
        center_angle = math.pi / 2 + index * 2 * math.pi / 3
        cx = candidate.lobe_center_radius_m * math.cos(center_angle)
        cy = candidate.lobe_center_radius_m * math.sin(center_angle)
        lines.append(f"{request.wire_radius_m:.12g} 1")
        lines.append(f"a {cx:.12g} {cy:.12g} 1 {request.turns + 1}")
        for turn in range(request.turns + 1):
            radius = candidate.inner_radius_m + turn * candidate.pitch_m
            lines.append(f"{cx + radius:.12g} {cy:.12g} 0")
    return "\n".join(lines) + "\n"
