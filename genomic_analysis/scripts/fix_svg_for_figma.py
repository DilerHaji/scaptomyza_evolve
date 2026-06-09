#!/usr/bin/env python3

import re
import sys
import argparse


def fix_svg(svg_text):
    n_stroke_added = 0
    n_rect_fixed = 0
    n_circle_fixed = 0

    def fix_line_style(m):
        nonlocal n_stroke_added
        style = m.group(1)
        if re.search(r"(?<![\w-])stroke\s*:", style):
            return m.group(0)
        if "stroke-width" in style:
            style = style.rstrip().rstrip(";") + "; stroke: black;"
            n_stroke_added += 1
            return f"style='{style}'"
        return m.group(0)
    svg_text = re.sub(r"style='([^']*)'", fix_line_style, svg_text)

    def fix_bare_rect(m):
        nonlocal n_rect_fixed
        tag = m.group(0)
        if "style" in tag or 'fill=' in tag or 'stroke=' in tag:
            return tag
        n_rect_fixed += 1
        return tag.rstrip("/>").rstrip() + " style='fill: none; stroke: none;' />"
    svg_text = re.sub(r"<rect\s+[^>]*?/>", fix_bare_rect, svg_text)

    def fix_circle(m):
        nonlocal n_circle_fixed
        tag = m.group(0)
        style_m = re.search(r"style='([^']*)'", tag)
        if not style_m:
            return tag
        style = style_m.group(1)
        if re.search(r"(?<![\w-])fill\s*:", style):
            return tag
        if re.search(r"(?<![\w-])stroke\s*:", style):
            new_style = style.rstrip().rstrip(";") + "; fill: none;"
            n_circle_fixed += 1
            return tag.replace(f"style='{style}'", f"style='{new_style}'")
        return tag
    svg_text = re.sub(r"<circle\s+[^>]*?/>", fix_circle, svg_text)

    sys.stderr.write(f"  Added 'stroke: black' to {n_stroke_added} lines/polylines\n")
    sys.stderr.write(f"  Added 'fill: none; stroke: none' to {n_rect_fixed} bare rects\n")
    sys.stderr.write(f"  Added 'fill: none' to {n_circle_fixed} open circles\n")

    return svg_text


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("output", nargs="?")
    p.add_argument("--in-place", action="store_true")
    args = p.parse_args()

    with open(args.input) as f:
        svg = f.read()

    fixed = fix_svg(svg)

    if args.in_place:
        out = args.input
    elif args.output:
        out = args.output
    else:
        sys.stderr.write("Need output path or --in-place\n")
        sys.exit(1)

    with open(out, "w") as f:
        f.write(fixed)

    sys.stderr.write(f"  Bytes: {len(svg)} → {len(fixed)}\n")


if __name__ == "__main__":
    main()
