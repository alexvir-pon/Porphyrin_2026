#!/usr/bin/env python3
"""
Render fig5_chromatin_iron_heme_model.svg to print-resolution PNG and PDF.

Figure 5 for the "Chromatin regulator perturbation causes heme pathway
imbalance in yeast and human cells" manuscript.

Usage:
    python render_fig5.py                 # 300 and 600 dpi PNG + PDF
    python render_fig5.py --dpi 1200      # custom DPI for the PNG
    python render_fig5.py --svg other.svg # render a different SVG

Dependencies:
    pip install cairosvg

Font note (IMPORTANT):
    The SVG requests "Arial". cairosvg renders using the fonts installed on
    THIS machine. If Arial is not present (common on Linux), the renderer
    falls back to a default sans-serif, which changes text metrics.
      - On Windows/macOS Arial is installed by default -> nothing to do.
      - On Linux, either install Arial (e.g. `ttf-mscorefonts-installer`)
        or install Liberation Sans, which is metric-compatible with Arial
        (`sudo apt-get install fonts-liberation`). The script checks for a
        usable Arial/Liberation face and warns if neither is found.
    For a font-independent, fully faithful vector you can also just submit
    the SVG itself, or convert text-to-outlines in Inkscape.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

DEFAULT_SVG = "fig5_chromatin_iron_heme_model.svg"


def check_fonts():
    """Warn if neither Arial nor a metric-compatible substitute is available."""
    try:
        if shutil.which("fc-list"):
            out = subprocess.run(
                ["fc-list"], capture_output=True, text=True, check=False
            ).stdout.lower()
            if "arial" in out:
                print("[fonts] Arial found — text will render as Arial.")
                return
            if "liberation sans" in out:
                print(
                    "[fonts] Arial not found, but Liberation Sans is installed "
                    "(metric-compatible with Arial) — layout will match."
                )
                return
            print(
                "[fonts] WARNING: neither Arial nor Liberation Sans found. "
                "Text will fall back to a default sans-serif and metrics may "
                "shift. Install one of them for a faithful render:\n"
                "        Debian/Ubuntu: sudo apt-get install fonts-liberation\n"
                "        (or ttf-mscorefonts-installer for real Arial)"
            )
        else:
            print("[fonts] fontconfig (fc-list) not available; skipping font check.")
    except Exception as exc:  # noqa: BLE001
        print(f"[fonts] font check skipped ({exc}).")


def render(svg_path: Path, dpis):
    try:
        import cairosvg
    except ImportError:
        sys.exit(
            "cairosvg is not installed. Install it with:\n"
            "    pip install cairosvg"
        )

    svg_bytes = svg_path.read_bytes()
    stem = svg_path.with_suffix("")  # drop .svg

    # PNG at each requested DPI.
    # The SVG has explicit width/height in px, so we derive a scale factor
    # relative to 96 dpi (SVG screen resolution) to get true DPI output.
    for dpi in dpis:
        scale = dpi / 96.0
        png_out = f"{stem}_{dpi}dpi.png"
        cairosvg.svg2png(
            bytestring=svg_bytes,
            write_to=png_out,
            scale=scale,
            background_color="white",
        )
        print(f"[png ] wrote {png_out}  ({dpi} dpi, scale={scale:.2f}x)")

    # PDF (vector — resolution-independent, ideal for submission)
    pdf_out = f"{stem}.pdf"
    cairosvg.svg2pdf(bytestring=svg_bytes, write_to=pdf_out)
    print(f"[pdf ] wrote {pdf_out}  (vector)")


def main():
    ap = argparse.ArgumentParser(description="Render Fig. 5 SVG to PNG and PDF.")
    ap.add_argument("--svg", default=DEFAULT_SVG, help="input SVG path")
    ap.add_argument(
        "--dpi",
        type=int,
        nargs="*",
        default=[300, 600],
        help="one or more DPI values for PNG output (default: 300 600)",
    )
    args = ap.parse_args()

    svg_path = Path(args.svg)
    if not svg_path.exists():
        sys.exit(f"SVG not found: {svg_path}")

    check_fonts()
    render(svg_path, args.dpi)
    print(
        "\nDone. For journal submission, the vector PDF (or the SVG itself) is best;\n"
        "the 600 dpi PNG is a safe raster fallback."
    )


if __name__ == "__main__":
    main()
