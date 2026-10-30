"""Produce a human-readable fabrication summary from an evaluated candidate."""

from .models import Candidate, DesignRequest, Evaluation


def fabrication_description(request: DesignRequest, candidate: Candidate, evaluation: Evaluation) -> str:
    radii = ", ".join(f"{value:.6g} m" for value in evaluation.wire_lengths_m)
    return (
        f"Three planar coils form a three-lobe flower, with their centres spaced at 120-degree intervals "
        f"on a circle of radius {candidate.lobe_center_radius_m:.6g} m. Each coil has {request.turns} outward-growing turns, "
        f"an inner winding radius of {candidate.inner_radius_m:.6g} m, and a radial pitch of {candidate.pitch_m:.6g} m. "
        f"The wire cross-section radius is {request.wire_radius_m:.6g} m. The calculated wire lengths are {radii}."
    )
