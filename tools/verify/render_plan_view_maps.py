from __future__ import annotations

import argparse
import json
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
    rows: list[list[float]], output_path: Path, title: str
) -> str:
    import matplotlib.pyplot as plt  # type: ignore

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    image = ax.imshow(rows, origin="lower", cmap="turbo")
    ax.set_title(title)
    ax.set_xlabel("i")
    ax.set_ylabel("j")
    fig.colorbar(image, ax=ax, shrink=0.85)
    fig.savefig(output_path, dpi=140)
    plt.close(fig)
    return "png"


def render_ppm(rows: list[list[float]], output_path: Path, scale: int) -> str:
    ny = len(rows)
    nx = len(rows[0]) if ny else 0
    values = [value for row in rows for value in row]
    vmin = min(values) if values else 0.0
    vmax = max(values) if values else 0.0
    out_nx = max(1, nx * scale)
    out_ny = max(1, ny * scale)

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("P3\n")
        handle.write(f"{out_nx} {out_ny}\n255\n")
        for src_j in reversed(range(ny)):
            for _ in range(scale):
                row_pixels: list[str] = []
                for src_i in range(nx):
                    frac = normalize(rows[src_j][src_i], vmin, vmax)
                    r, g, b = colormap_rgb(frac)
                    pixel = f"{r} {g} {b}"
                    row_pixels.extend([pixel] * scale)
                handle.write(" ".join(row_pixels))
                handle.write("\n")
    return "ppm"


def render_field(
    field: dict[str, Any], output_dir: Path, allow_matplotlib: bool, scale: int, title_prefix: str
) -> dict[str, Any]:
    rows = reshape_field(field)
    title = f"{title_prefix} {field['name']} [{field.get('units', '')}]".strip()
    stem = str(field["name"])
    if allow_matplotlib:
        try:
            output_path = output_dir / f"{stem}.png"
            fmt = render_with_matplotlib(rows, output_path, title)
            return {
                "field": stem,
                "units": field.get("units", ""),
                "path": str(output_path.resolve()),
                "format": fmt,
            }
        except ModuleNotFoundError:
            pass

    output_path = output_dir / f"{stem}.ppm"
    fmt = render_ppm(rows, output_path, scale)
    return {
        "field": stem,
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
