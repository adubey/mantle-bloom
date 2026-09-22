#!/usr/bin/env python3
"""Render the issue #213 paired sweep with Pillow fonts and labeled axes."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

DEFAULT_SUMMARY = Path(__file__).resolve().parent / "results" / "issue213_hc_hm.summary.json"
FONT_PATH = "/System/Library/Fonts/Supplemental/Arial.ttf"


def _font(size: int):
    try:
        return ImageFont.truetype(FONT_PATH, size)
    except OSError:
        return ImageFont.load_default()


def _nice_step(span: float) -> float:
    raw = max(span / 4, 1e-9)
    magnitude = 10 ** math.floor(math.log10(raw))
    for scale in (1, 2, 2.5, 5, 10):
        if raw <= scale * magnitude:
            return scale * magnitude
    return 10 * magnitude


def render(summary: dict, output: Path) -> None:
    panels = (
        ("hc_mean_m", "Mean Hc", "km", 0.001),
        ("hc_p95_m", "95th percentile Hc", "km", 0.001),
        ("hc_at_max_fraction", "Hc at maximum", "% of nodes", 100),
        ("hm_mean_m", "Mean Hm", "km", 0.001),
        ("hm_p95_m", "95th percentile Hm", "km", 0.001),
        ("hm_at_max_fraction", "Hm at maximum", "% of nodes", 100),
    )
    points = [p for p in summary["checkpoints"] if p["paired_seeds"]]
    if not points:
        raise ValueError("summary has no paired checkpoints")
    ages = [p["checkpoint_years"] / 1e6 for p in points]
    max_age = max(ages)
    scale = 2
    width, height = 1500, 1180
    img = Image.new("RGB", (width * scale, height * scale), "#f6f8fc")
    draw = ImageDraw.Draw(img)
    title_font = _font(31 * scale)
    subtitle_font = _font(19 * scale)
    panel_font = _font(22 * scale)
    axis_font = _font(16 * scale)
    note_font = _font(16 * scale)
    old_color, new_color = "#2167b2", "#d66a21"

    def text(x, y, value, font, fill="#16263c", anchor=None):
        draw.text((round(x * scale), round(y * scale)), value, font=font, fill=fill, anchor=anchor)

    def line(coords, fill, w=1):
        draw.line([(round(x * scale), round(y * scale)) for x, y in coords], fill=fill, width=w * scale, joint="curve")

    text(52, 30, "Fault relief: old additive vs conserved", title_font)
    text(52, 73, f"30 paired seeds  •  10 My steps  •  0–{max_age:g} My", subtitle_font, "#53647b")
    line([(1050, 84), (1090, 84)], old_color, 4)
    text(1100, 73, "Old additive", subtitle_font)
    line([(1280, 84), (1320, 84)], new_color, 4)
    text(1330, 73, "Conserved", subtitle_font)

    for i, (metric, title, unit, factor) in enumerate(panels):
        col, row = i % 2, i // 2
        left = 52 + col * 735
        top = 138 + row * 338
        right = left + 678
        bottom = top + 276
        plot_left, plot_right = left + 82, right - 26
        plot_top, plot_bottom = top + 40, bottom - 58
        draw.rounded_rectangle((left * scale, top * scale, right * scale, bottom * scale),
                               radius=13 * scale, fill="white", outline="#dce3ed", width=2)
        text(left + 20, top + 12, f"{title} ({unit})", panel_font)

        old = [p["metrics"][metric]["old_mean"] * factor for p in points]
        new = [p["metrics"][metric]["new_mean"] * factor for p in points]
        lo, hi = min(old + new), max(old + new)
        pad = max((hi - lo) * 0.08, 0.03 if unit == "% of nodes" else 0.3)
        step = _nice_step((hi - lo) + 2 * pad)
        y_min = math.floor((lo - pad) / step) * step
        y_max = math.ceil((hi + pad) / step) * step
        if unit == "% of nodes":
            y_min = max(0.0, y_min)
        if y_max <= y_min:
            y_max = y_min + step

        def xy(age, value):
            return (plot_left + age / max_age * (plot_right - plot_left),
                    plot_bottom - (value - y_min) / (y_max - y_min) * (plot_bottom - plot_top))

        ticks = int(round((y_max - y_min) / step))
        for j in range(ticks + 1):
            value = y_min + j * step
            y = xy(0, value)[1]
            line([(plot_left, y), (plot_right, y)], "#e8edf4")
            label = f"{value:.2f}" if step < 0.1 else (f"{value:.1f}" if step < 1 else f"{value:g}")
            text(plot_left - 12, y, label, axis_font, "#53647b", anchor="rm")
        for age in (0, 60, 120, 180, 240, 320):
            if age > max_age:
                continue
            x = xy(age, y_min)[0]
            line([(x, plot_bottom), (x, plot_bottom + 5)], "#7e8da2")
            text(x, plot_bottom + 11, str(age), axis_font, "#53647b", anchor="mt")
        text((plot_left + plot_right) / 2, bottom - 12, "Age (My)", axis_font, "#53647b", anchor="mm")
        for values, color in ((old, old_color), (new, new_color)):
            line([xy(age, value) for age, value in zip(ages, values)], color, 4)
            for age, value in zip(ages, values):
                x, y = xy(age, value)
                r = 3.2 * scale
                draw.ellipse((x * scale - r, y * scale - r, x * scale + r, y * scale + r), fill=color)
        delta = points[-1]["metrics"][metric]["paired_delta_mean"] * factor
        sign = "+" if delta >= 0 else "−"
        text(right - 18, top + 20, f"Δ320: {sign}{abs(delta):.2f} {unit if unit != '% of nodes' else '%'}", note_font,
             "#53647b", anchor="ra")

    output.parent.mkdir(parents=True, exist_ok=True)
    img.resize((width, height), Image.Resampling.LANCZOS).save(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    output = args.out or args.summary.with_suffix(".png")
    render(json.loads(args.summary.read_text()), output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
