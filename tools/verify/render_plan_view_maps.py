from __future__ import annotations

import argparse
import json
import math
import re
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def normalize(value: float, vmin: float, vmax: float) -> float:
    if vmax <= vmin:
        return 0.5
    return max(0.0, min(1.0, (value - vmin) / (vmax - vmin)))


def colormap_rgb(frac: float) -> tuple[int, int, int]:
    frac = max(0.0, min(1.0, frac))
    if frac < 0.25:
        local = frac / 0.25
        return (0, int(255 * local), 255)
    if frac < 0.5:
        local = (frac - 0.25) / 0.25
        return (0, 255, int(255 * (1.0 - local)))
    if frac < 0.75:
        local = (frac - 0.5) / 0.25
        return (int(255 * local), 255, 0)
    local = (frac - 0.75) / 0.25
    return (255, int(255 * (1.0 - local)), 0)


def sanitize_stem(name: str) -> str:
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return stem.strip("._") or "field"


def flattened_values(rows: list[list[float]]) -> list[float]:
    return [value for row in rows for value in row]


def percentile(values: list[float], frac: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    frac = max(0.0, min(1.0, frac))
    pos = frac * (len(ordered) - 1)
    lower = int(math.floor(pos))
    upper = int(math.ceil(pos))
    if lower == upper:
        return ordered[lower]
    alpha = pos - lower
    return ordered[lower] * (1.0 - alpha) + ordered[upper] * alpha


def bilinear_upsample(rows: list[list[float]], scale: int) -> list[list[float]]:
    ny = len(rows)
    nx = len(rows[0]) if ny else 0
    if nx == 0 or ny == 0 or scale <= 1:
        return [list(row) for row in rows]

    out_nx = nx * scale
    out_ny = ny * scale
    output: list[list[float]] = []
    for out_j in range(out_ny):
        src_y = out_j / max(scale, 1)
        y0 = min(int(math.floor(src_y)), ny - 1)
        y1 = min(y0 + 1, ny - 1)
        ty = src_y - y0
        out_row: list[float] = []
        for out_i in range(out_nx):
            src_x = out_i / max(scale, 1)
            x0 = min(int(math.floor(src_x)), nx - 1)
            x1 = min(x0 + 1, nx - 1)
            tx = src_x - x0
            v00 = rows[y0][x0]
            v10 = rows[y0][x1]
            v01 = rows[y1][x0]
            v11 = rows[y1][x1]
            top = v00 * (1.0 - tx) + v10 * tx
            bottom = v01 * (1.0 - tx) + v11 * tx
            out_row.append(top * (1.0 - ty) + bottom * ty)
        output.append(out_row)
    return output


def choose_field_style(field: dict[str, Any], rows: list[list[float]]) -> dict[str, Any]:
    name = str(field.get("name", "")).lower()
    values = flattened_values(rows)
    vmin = min(values) if values else 0.0
    vmax = max(values) if values else 0.0
    style: dict[str, Any] = {
        "cmap": "turbo",
        "interpolation": "bilinear",
        "vmin": None,
        "vmax": None,
    }

    if any(key in name for key in ("u_velocity", "v_velocity", "w_velocity", "wind_component")):
        limit = max(percentile([abs(value) for value in values], 0.98), 1.0e-6)
        style.update({"cmap": "RdBu_r", "vmin": -limit, "vmax": limit})
        return style
    if name in {"terrain_height", "z_center", "geopotential_height"}:
        style.update({"cmap": "terrain"})
        return style
    if "reflectivity" in name:
        style.update({"cmap": "turbo", "vmin": -10.0, "vmax": 75.0})
        return style
    if "column_rain_water" in name:
        style.update({"cmap": "viridis", "vmin": 0.0, "vmax": percentile(values, 0.98)})
        return style
    if "relative_humidity" in name:
        style.update({"cmap": "Blues", "vmin": 0.0, "vmax": max(vmax, 100.0)})
        return style
    if "specific_humidity" in name:
        style.update({"cmap": "YlGnBu", "vmin": max(0.0, vmin)})
        return style
    if "cloud_water" in name or "rain_water" in name or "condensate" in name:
        style.update({"cmap": "PuBuGn", "vmin": 0.0, "vmax": percentile(values, 0.99)})
        return style
    if "dewpoint" in name or "temperature" in name or "theta" in name:
        style.update({"cmap": "coolwarm"})
        return style
    if "pressure" in name or "rho" in name:
        style.update({"cmap": "cividis"})
        return style
    if "wind_speed" in name:
        style.update({"cmap": "magma", "vmin": 0.0, "vmax": percentile(values, 0.98)})
        return style
    if values and all(math.isfinite(value) for value in values) and vmin >= 0.0:
        style["vmin"] = percentile(values, 0.02) if any(value < 0.0 for value in values) else max(0.0, percentile(values, 0.02))
        style["vmax"] = percentile(values, 0.98)
        if style["vmax"] is not None and style["vmin"] is not None and style["vmax"] <= style["vmin"]:
            style["vmin"] = 0.0
            style["vmax"] = None
    return style


def reshape_field(field: dict[str, Any]) -> list[list[float]]:
    nx = int(field["nx"])
    ny = int(field["ny"])
    values = field["values"]
    if len(values) != nx * ny:
        raise ValueError(f"Field {field['name']} has {len(values)} values; expected {nx * ny}")
    rows: list[list[float]] = []
    for j in range(ny):
        start = j * nx
        rows.append([float(value) for value in values[start : start + nx]])
    return rows


def render_with_matplotlib(
    rows: list[list[float]], output_path: Path, title: str, field: dict[str, Any]
) -> str:
    import matplotlib.pyplot as plt  # type: ignore

    style = choose_field_style(field, rows)
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    image = ax.imshow(
        rows,
        origin="lower",
        cmap=style["cmap"],
        interpolation=style["interpolation"],
        vmin=style["vmin"],
        vmax=style["vmax"],
    )
    ax.set_title(title)
    ax.set_xlabel("i")
    ax.set_ylabel("j")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.85)
    if field.get("units"):
        colorbar.set_label(str(field["units"]))
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return "png"


def render_ppm(rows: list[list[float]], output_path: Path, scale: int) -> str:
    upsampled = bilinear_upsample(rows, scale)
    ny = len(upsampled)
    nx = len(upsampled[0]) if ny else 0
    values = flattened_values(upsampled)
    vmin = min(values) if values else 0.0
    vmax = max(values) if values else 0.0

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("P3\n")
        handle.write(f"{nx} {ny}\n255\n")
        for src_j in reversed(range(ny)):
            row_pixels: list[str] = []
            for src_i in range(nx):
                frac = normalize(upsampled[src_j][src_i], vmin, vmax)
                r, g, b = colormap_rgb(frac)
                row_pixels.append(f"{r} {g} {b}")
            handle.write(" ".join(row_pixels))
            handle.write("\n")
    return "ppm"


def render_field(
    field: dict[str, Any], output_dir: Path, allow_matplotlib: bool, scale: int, title_prefix: str
) -> dict[str, Any]:
    rows = reshape_field(field)
    title = f"{title_prefix} {field['name']} [{field.get('units', '')}]".strip()
    stem = sanitize_stem(str(field["name"]))
    if allow_matplotlib:
        try:
            output_path = output_dir / f"{stem}.png"
            fmt = render_with_matplotlib(rows, output_path, title, field)
            return {
                "field": field.get("name", stem),
                "file_stem": stem,
                "units": field.get("units", ""),
                "path": str(output_path.resolve()),
                "format": fmt,
            }
        except ModuleNotFoundError:
            pass

    output_path = output_dir / f"{stem}.ppm"
    fmt = render_ppm(rows, output_path, scale)
    return {
        "field": field.get("name", stem),
        "file_stem": stem,
        "units": field.get("units", ""),
        "path": str(output_path.resolve()),
        "format": fmt,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render actual map images from a gwm-next dry plan-view JSON bundle."
    )
    parser.add_argument("--input", required=True, help="Path to gwm-next plan-view JSON")
    parser.add_argument("--output-dir", default=None, help="Directory for rendered image products")
    parser.add_argument("--fields", nargs="*", default=None, help="Optional subset of field names")
    parser.add_argument("--title-prefix", default="", help="Optional title prefix for rendered figures")
    parser.add_argument("--no-matplotlib", action="store_true", help="Force fallback PPM rendering")
    parser.add_argument("--scale", type=int, default=4, help="Pixel upscale for fallback PPM output")
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    payload = load_json(input_path)
    if payload.get("schema_version") != "gwm-next-plan-view/v1":
        raise SystemExit(f"Unsupported plan-view schema: {payload.get('schema_version', '')}")

    fields = payload.get("fields", [])
    if args.fields:
        wanted = set(args.fields)
        fields = [field for field in fields if field.get("name") in wanted]
    if not fields:
        raise SystemExit("No fields selected for rendering")

    default_output_dir = input_path.with_suffix("")
    output_dir = Path(args.output_dir).resolve() if args.output_dir else default_output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    rendered = [
        render_field(field, output_dir, not args.no_matplotlib, args.scale, args.title_prefix)
        for field in fields
    ]

    manifest = {
        "schema_version": "gwm-next-map-manifest/v1",
        "generated_at_utc": utc_now(),
        "input_plan_view_path": str(input_path),
        "case": payload.get("case", ""),
        "slice_k": payload.get("slice_k", 0),
        "grid": payload.get("grid", {}),
        "images": rendered,
    }

    manifest_path = output_dir / "map_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"Rendered {len(rendered)} maps from: {input_path}")
    print(f"Manifest: {manifest_path.resolve()}")
    for item in rendered:
        print(f"- {item['field']}: {item['path']}")


if __name__ == "__main__":
    main()
