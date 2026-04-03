from __future__ import annotations

import json
from pathlib import Path


KNOWN_TOOLS = {
    "cfrust": [Path(r"C:\Users\drew\cfrust")],
    "ecrust": [Path(r"C:\Users\drew\ecrust")],
    "wrf-rust": [Path(r"C:\Users\drew\wrf-rust")],
    "wrf-rust-plots": [Path(r"C:\Users\drew\wrf-rust-plots")],
    "metrust-or-rustmet": [
        Path(r"C:\Users\drew\metrust"),
        Path(r"C:\Users\drew\rustmet"),
        Path(r"C:\Users\drew\metrust-py"),
    ],
    "rusbie-or-rustbie": [
        Path(r"C:\Users\drew\rusbie"),
        Path(r"C:\Users\drew\rustbie"),
    ],
}


def main() -> None:
    status = {}
    for name, candidates in KNOWN_TOOLS.items():
        matches = [str(path) for path in candidates if path.exists()]
        status[name] = {
            "candidates": [str(path) for path in candidates],
            "matches": matches,
            "exists": bool(matches),
        }
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
