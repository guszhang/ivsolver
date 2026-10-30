"""Validated input and result models for the coil design workflow."""

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class DesignRequest:
    target_self_inductance_h: float
    wire_radius_m: float
    resistivity_ohm_m: float
    turns: int
    inner_radius_bounds_m: tuple[float, float]
    pitch_bounds_m: tuple[float, float]
    lobe_center_radius_bounds_m: tuple[float, float]
    inductance_tolerance_h: float
    mutual_tolerance_h: float
    samples_per_axis: int = 3
    refine_rounds: int = 1
    max_footprint_radius_m: float | None = None

    def validate(self) -> None:
        positive = {
            "target_self_inductance_h": self.target_self_inductance_h,
            "wire_radius_m": self.wire_radius_m,
            "resistivity_ohm_m": self.resistivity_ohm_m,
            "inductance_tolerance_h": self.inductance_tolerance_h,
            "mutual_tolerance_h": self.mutual_tolerance_h,
        }
        for name, value in positive.items():
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be a finite positive number")
        if not 1 <= self.turns <= 20:
            raise ValueError("turns must be between 1 and 20")
        if not 3 <= self.samples_per_axis <= 5:
            raise ValueError("samples_per_axis must be between 3 and 5")
        if not 0 <= self.refine_rounds <= 2:
            raise ValueError("refine_rounds must be between 0 and 2")
        for name, bounds in (
            ("inner_radius_bounds_m", self.inner_radius_bounds_m),
            ("pitch_bounds_m", self.pitch_bounds_m),
            ("lobe_center_radius_bounds_m", self.lobe_center_radius_bounds_m),
        ):
            if len(bounds) != 2 or not all(math.isfinite(value) for value in bounds) or bounds[0] <= 0 or bounds[1] < bounds[0]:
                raise ValueError(f"{name} must be an increasing pair of positive metres")
        if self.max_footprint_radius_m is not None and (not math.isfinite(self.max_footprint_radius_m) or self.max_footprint_radius_m <= 0):
            raise ValueError("max_footprint_radius_m must be finite and positive when supplied")


@dataclass(frozen=True)
class Candidate:
    inner_radius_m: float
    pitch_m: float
    lobe_center_radius_m: float


@dataclass
class Evaluation:
    inductance_matrix_h: list[list[float]]
    wire_lengths_m: list[float]
    trajectories_m: list[list[tuple[float, float, float]]]
