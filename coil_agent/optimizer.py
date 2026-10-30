"""Bounded deterministic search for the symmetric three-lobe starter topology."""

import itertools
import math

from .description import fabrication_description
from .geometry import flower_wire_definition
from .models import Candidate, DesignRequest
from .render import render_png, render_svg
from .solver import evaluate

MAX_BACKEND_EVALUATIONS = 150


def _grid(low: float, high: float, count: int) -> list[float]:
    if count == 1 or low == high:
        return [low]
    return [low + (high - low) * i / (count - 1) for i in range(count)]


def optimize(request: DesignRequest) -> dict:
    request.validate()
    ranges = [
        list(request.inner_radius_bounds_m),
        list(request.pitch_bounds_m),
        list(request.lobe_center_radius_bounds_m),
    ]
    best = None
    evaluated = 0
    attempted = 0
    errors = []

    for round_index in range(request.refine_rounds + 1):
        axes = [_grid(low, high, request.samples_per_axis) for low, high in ranges]
        for inner_radius, pitch, center_radius in itertools.product(*axes):
            if attempted >= MAX_BACKEND_EVALUATIONS:
                break
            candidate = Candidate(inner_radius, pitch, center_radius)
            outer_radius = inner_radius + request.turns * pitch
            if inner_radius < request.wire_radius_m or pitch < 2 * request.wire_radius_m:
                continue
            if request.max_footprint_radius_m is not None and center_radius + outer_radius > request.max_footprint_radius_m:
                continue
            center_spacing = math.sqrt(3.0) * center_radius
            if center_spacing < 2 * (outer_radius + request.wire_radius_m):
                continue
            try:
                attempted += 1
                result = evaluate(request, candidate)
                if not all(math.isfinite(value) for row in result.inductance_matrix_h for value in row):
                    errors.append("solver returned non-finite inductance for a candidate")
                    continue
                evaluated += 1
                diagonal = [result.inductance_matrix_h[i][i] for i in range(3)]
                mutual = [result.inductance_matrix_h[i][j] for i in range(3) for j in range(i + 1, 3)]
                resistances = [
                    request.resistivity_ohm_m * length / (math.pi * request.wire_radius_m**2)
                    for length in result.wire_lengths_m
                ]
                ratios = [inductance / resistance for inductance, resistance in zip(diagonal, resistances)]
                inductance_error = max(abs(value - request.target_self_inductance_h) for value in diagonal)
                mutual_error = max(abs(value) for value in mutual)
                feasible = inductance_error <= request.inductance_tolerance_h and mutual_error <= request.mutual_tolerance_h
                violation = (inductance_error / request.inductance_tolerance_h) ** 2 + (mutual_error / request.mutual_tolerance_h) ** 2
                rank = (0, -sum(ratios) / len(ratios)) if feasible else (1, violation, -sum(ratios) / len(ratios))
                if best is None or rank < best["rank"]:
                    best = {
                        "rank": rank,
                        "candidate": candidate,
                        "evaluation": result,
                        "diagonal": diagonal,
                        "mutual": mutual,
                        "resistances": resistances,
                        "ratios": ratios,
                        "inductance_error": inductance_error,
                        "mutual_error": mutual_error,
                        "feasible": feasible,
                    }
            except (OSError, RuntimeError, ValueError, TimeoutError) as exc:
                errors.append(str(exc))

        if attempted >= MAX_BACKEND_EVALUATIONS:
            break

        if best is None:
            break
        if round_index < request.refine_rounds:
            point = (best["candidate"].inner_radius_m, best["candidate"].pitch_m, best["candidate"].lobe_center_radius_m)
            refined = []
            for (low, high), value in zip(ranges, point):
                half_width = (high - low) / max(2, request.samples_per_axis - 1)
                refined.append([max(low, value - half_width), min(high, value + half_width)])
            ranges = refined

    if best is None:
        raise RuntimeError("No candidate could be evaluated. Check the solver binary and geometry bounds. " + "; ".join(errors[-3:]))

    candidate = best["candidate"]
    result = best["evaluation"]
    return {
        "feasible": best["feasible"],
        "topology": "three identical outward-growing spiral coils arranged at 120-degree intervals in the XY plane",
        "candidate": {
            "inner_radius_m": candidate.inner_radius_m,
            "turn_pitch_m": candidate.pitch_m,
            "lobe_center_radius_m": candidate.lobe_center_radius_m,
            "outer_radius_m": candidate.inner_radius_m + request.turns * candidate.pitch_m,
        },
        "inductance_matrix_h": result.inductance_matrix_h,
        "self_inductances_uH": [value * 1e6 for value in best["diagonal"]],
        "mutual_inductances_nH": [value * 1e9 for value in best["mutual"]],
        "wire_lengths_m": result.wire_lengths_m,
        "resistances_ohm": best["resistances"],
        "inductance_over_resistance_h_per_ohm": best["ratios"],
        "max_self_inductance_error_h": best["inductance_error"],
        "max_mutual_inductance_abs_h": best["mutual_error"],
        "searched_candidates": evaluated,
        "search_truncated": attempted >= MAX_BACKEND_EVALUATIONS,
        "constraints": {
            "target_self_inductance_h": request.target_self_inductance_h,
            "inductance_tolerance_h": request.inductance_tolerance_h,
            "mutual_tolerance_h": request.mutual_tolerance_h,
            "wire_radius_m": request.wire_radius_m,
            "resistivity_ohm_m": request.resistivity_ohm_m,
            "turns": request.turns,
            "inner_radius_bounds_m": list(request.inner_radius_bounds_m),
            "pitch_bounds_m": list(request.pitch_bounds_m),
            "lobe_center_radius_bounds_m": list(request.lobe_center_radius_bounds_m),
            "max_footprint_radius_m": request.max_footprint_radius_m,
            "samples_per_axis": request.samples_per_axis,
            "refine_rounds": request.refine_rounds,
        },
        "coil_description": fabrication_description(request, candidate, result),
        "wire_definition": flower_wire_definition(request, candidate),
        "rendering_svg": render_svg(result),
        "rendering_png": render_png(result),
        "warnings": [] if best["feasible"] else ["No candidate met all inductance and mutual-coupling tolerances in the searched bounds."],
    }
