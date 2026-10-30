"""Run the existing solver binary and parse its inductance/trajectory outputs."""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

from .geometry import flower_wire_definition
from .models import Candidate, DesignRequest, Evaluation


def default_solver_path() -> Path:
    configured = os.environ.get("IVSOLVER_BIN")
    if configured:
        return Path(configured).resolve()
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        bundled = Path(sys._MEIPASS) / "ivsolver-agent-backend"
        if bundled.exists():
            return bundled
    root = Path(__file__).resolve().parents[1]
    agent_backend = root / "build" / "ivsolver-agent-backend"
    return agent_backend if agent_backend.exists() else root / "build" / "ivsolver"


def evaluate(request: DesignRequest, candidate: Candidate, solver_path: Path | None = None) -> Evaluation:
    binary = (solver_path or default_solver_path()).resolve()
    if not binary.is_file():
        raise FileNotFoundError(f"IVSolver backend binary not found: {binary}")

    with tempfile.TemporaryDirectory(prefix="ivsolver-design-") as temp_dir:
        work = Path(temp_dir)
        definition = work / "wires.def"
        output_prefix = work / "result"
        definition.write_text(flower_wire_definition(request, candidate), encoding="utf-8")
        completed = subprocess.run(
            [str(binary), str(definition), str(output_prefix)],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        if completed.returncode:
            detail = completed.stderr.strip() or completed.stdout.strip() or "solver failed without a diagnostic"
            raise RuntimeError(f"IVSolver exited with {completed.returncode}: {detail}")

        rows = [[float(value) for value in line.split()] for line in Path(f"{output_prefix}.lmat").read_text().splitlines()]
        matrix = [row[:3] for row in rows[:3]]
        lengths = rows[3][:3]
        trajectories = []
        for index in range(3):
            points = [tuple(float(value) for value in line.split()) for line in Path(f"{output_prefix}_{index}.trj").read_text().splitlines()]
            trajectories.append(points)
        return Evaluation(matrix, lengths, trajectories)
