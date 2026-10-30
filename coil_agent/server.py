"""MCP server exposing the bounded three-lobe coil design optimizer."""

from mcp.server.fastmcp import FastMCP
from mcp.server.fastmcp.utilities.types import Image
from mcp.types import CallToolResult, TextContent

from .models import DesignRequest
from .optimizer import optimize
from .render import render_png


mcp = FastMCP(
    "ivsolver-coil-design",
    instructions=(
        "Design and evaluate a symmetric three-coil flower layout. Use design_requirements to identify "
        "missing wire, material, geometry-bound, and tolerance inputs. Ask for missing values unless "
        "the user authorizes reasonable defaults; disclose every selected value. Results are solver-evaluated, "
        "and an infeasible result must be reported as such."
    ),
)


@mcp.tool()
def design_requirements() -> dict:
    """List the physical and geometric inputs needed before designing a coil."""
    return {
        "known_from_example": {
            "coil_count": 3,
            "turns": 3,
            "target_self_inductance_uH": 10.0,
            "layout": "three identical outward-growing coils at 120 degrees in one plane",
        },
        "required_from_user_or_explicit_assumptions": [
            {"parameter": "wire_radius_mm", "meaning": "wire cross-section radius"},
            {"parameter": "resistivity_ohm_m", "meaning": "conductor DC resistivity"},
            {"parameter": "inner_radius_min_mm and inner_radius_max_mm", "meaning": "search bounds for the innermost turn"},
            {"parameter": "turn_pitch_min_mm and turn_pitch_max_mm", "meaning": "search bounds for radial spacing between turns"},
            {"parameter": "lobe_center_radius_min_mm and lobe_center_radius_max_mm", "meaning": "search bounds for lobe placement"},
            {"parameter": "inductance_tolerance_uH", "meaning": "allowed difference from the target self-inductance"},
            {"parameter": "mutual_tolerance_nH", "meaning": "maximum accepted magnitude of each mutual inductance; exact zero is represented by a tolerance"},
        ],
        "optional": ["max_footprint_radius_mm"],
        "policy": "Ask for missing values. If the user authorizes defaults, state every selected value in the returned assumptions.",
    }


@mcp.tool()
def design_coil(
    target_self_inductance_uH: float,
    wire_radius_mm: float,
    resistivity_ohm_m: float,
    inner_radius_min_mm: float,
    inner_radius_max_mm: float,
    turn_pitch_min_mm: float,
    turn_pitch_max_mm: float,
    lobe_center_radius_min_mm: float,
    lobe_center_radius_max_mm: float,
    inductance_tolerance_uH: float,
    mutual_tolerance_nH: float,
    turns: int = 3,
    samples_per_axis: int = 3,
    refine_rounds: int = 1,
    max_footprint_radius_mm: float | None = None,
) -> CallToolResult:
    """Optimize a symmetric three-lobe planar coil system and return metrics plus SVG.

    Each lobe is one outward-growing spiral. The three spirals have equal geometry
    and their centers lie 120 degrees apart. The bounded search maximizes L/R while
    trying to meet the specified self-inductance and mutual-inductance tolerances.
    Length inputs are millimetres; resistivity is ohm-metres.
    """
    mm = 1e-3
    request = DesignRequest(
        target_self_inductance_h=target_self_inductance_uH * 1e-6,
        wire_radius_m=wire_radius_mm * mm,
        resistivity_ohm_m=resistivity_ohm_m,
        turns=turns,
        inner_radius_bounds_m=(inner_radius_min_mm * mm, inner_radius_max_mm * mm),
        pitch_bounds_m=(turn_pitch_min_mm * mm, turn_pitch_max_mm * mm),
        lobe_center_radius_bounds_m=(lobe_center_radius_min_mm * mm, lobe_center_radius_max_mm * mm),
        inductance_tolerance_h=inductance_tolerance_uH * 1e-6,
        mutual_tolerance_h=mutual_tolerance_nH * 1e-9,
        samples_per_axis=samples_per_axis,
        refine_rounds=refine_rounds,
        max_footprint_radius_m=None if max_footprint_radius_mm is None else max_footprint_radius_mm * mm,
    )
    result = optimize(request)
    png = result.pop("rendering_png")
    result.pop("rendering_svg", None)
    text = (
        f"Feasible: {result['feasible']}\n"
        f"Self inductances (uH): {result['self_inductances_uH']}\n"
        f"Mutual inductances (nH): {result['mutual_inductances_nH']}\n"
        f"L/R (H/ohm): {result['inductance_over_resistance_h_per_ohm']}\n"
        f"Coil description: {result['coil_description']}\n"
        f"Warnings: {result['warnings']}"
    )
    return CallToolResult(
        content=[TextContent(type="text", text=text), Image(data=png, format="png").to_image_content()],
        structuredContent=result,
    )


def main() -> None:
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
