from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Verification skeleton for gpu-wm-next outputs."
    )
    parser.add_argument("--input", required=True, help="Model output path")
    parser.add_argument("--mode", default="screen", choices=["screen", "aloft", "products"])
    args = parser.parse_args()

    model_output = Path(args.input)
    print(f"Verification skeleton only. input={model_output.resolve()} mode={args.mode}")
    print("Expected downstream tools:")
    print("- wrf-rust for diagnostics / verification")
    print("- wrf-rust-plots for products")
    print("- Rust MetPy replacement for generic thermo / severe math")


if __name__ == "__main__":
    main()
