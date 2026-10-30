"""Render generated planar trajectories to a self-contained SVG."""

from io import BytesIO

from PIL import Image, ImageDraw

from .models import Evaluation


COLORS = ("#2563eb", "#dc2626", "#16a34a")


def render_svg(evaluation: Evaluation, width: int = 900, height: int = 700) -> str:
    all_points = [point for trajectory in evaluation.trajectories_m for point in trajectory]
    min_x = min(point[0] for point in all_points)
    max_x = max(point[0] for point in all_points)
    min_y = min(point[1] for point in all_points)
    max_y = max(point[1] for point in all_points)
    span = max(max_x - min_x, max_y - min_y, 1e-9)
    pad = span * 0.1
    scale = min((width - 2 * pad / span * width) / span, (height - 2 * pad / span * height) / span)
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2

    paths = []
    for index, trajectory in enumerate(evaluation.trajectories_m):
        coords = [
            (width / 2 + (x - center_x) * scale, height / 2 - (y - center_y) * scale)
            for x, y, _ in trajectory
        ]
        path = " ".join(("M" if i == 0 else "L") + f"{x:.2f},{y:.2f}" for i, (x, y) in enumerate(coords))
        paths.append(f'<path d="{path}" fill="none" stroke="{COLORS[index]}" stroke-width="1.5"/>')
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
        '<rect width="100%" height="100%" fill="white"/>'
        + "".join(paths)
        + "</svg>"
    )


def render_png(evaluation: Evaluation, width: int = 1200, height: int = 900) -> bytes:
    """Rasterize trajectories as a PNG for MCP clients that display image content."""
    all_points = [point for trajectory in evaluation.trajectories_m for point in trajectory]
    min_x = min(point[0] for point in all_points)
    max_x = max(point[0] for point in all_points)
    min_y = min(point[1] for point in all_points)
    max_y = max(point[1] for point in all_points)
    span = max(max_x - min_x, max_y - min_y, 1e-9)
    scale = min(width * 0.8 / span, height * 0.8 / span)
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    for index, trajectory in enumerate(evaluation.trajectories_m):
        coordinates = [
            (round(width / 2 + (x - center_x) * scale), round(height / 2 - (y - center_y) * scale))
            for x, y, _ in trajectory
        ]
        if len(coordinates) > 1:
            draw.line(coordinates, fill=COLORS[index], width=4, joint="curve")
    output = BytesIO()
    image.save(output, format="PNG", optimize=True)
    return output.getvalue()
