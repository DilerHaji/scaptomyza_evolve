#!/usr/bin/env python3
import re
import sys
import argparse


def strip_clippaths(svg_text):
    svg_text = re.sub(r"<clipPath\b[^>]*>.*?</clipPath>", "", svg_text, flags=re.DOTALL)
    svg_text = re.sub(r"\s+clip-path\s*=\s*['\"]url\(#[^)]+\)['\"]", "", svg_text)
    return svg_text


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("output", nargs="?", help="Output path; omit with --in-place")
    p.add_argument("--in-place", action="store_true",
                   help="Overwrite the input file")
    args = p.parse_args()

    with open(args.input) as f:
        svg = f.read()

    cleaned = strip_clippaths(svg)

    if args.in_place:
        out = args.input
    elif args.output:
        out = args.output
    else:
        sys.stderr.write("Need either output path or --in-place\n")
        sys.exit(1)

    with open(out, "w") as f:
        f.write(cleaned)

    n_cp = svg.count("<clipPath")
    n_attr = len(re.findall(r"clip-path\s*=", svg))
    sys.stderr.write(f"  Removed {n_cp} clipPath definitions and {n_attr} clip-path attributes\n")
    sys.stderr.write(f"  Input:  {len(svg)} bytes\n")
    sys.stderr.write(f"  Output: {len(cleaned)} bytes\n")


if __name__ == "__main__":
    main()
