#!/usr/bin/env python3
"""
Convert an IFU center CSV (e.g., temp/IFU_cen_unit2G.csv) into the
ASCII table format used by virus2_ifucen_<UNIT>.txt.

Usage:
  python scripts/antigen_convert_ifucen_csv.py -i temp/IFU_cen_unit2G.csv -o temp/virus2_ifucen_D2G.txt

CSV expected columns (header row):
  slitID, headID, type, px, subpx, ra, dec,

Output columns (fixed-width, matching Antigen config files):
  fiber_id head_id   ifu_x    ifu_y trace_row exclude_fiber

Mapping rules:
- fiber_id: sequential index starting at 1, preserving CSV order (after header)
- head_id: CSV headID as-is
- ifu_x, ifu_y:
  * type == 'obj' -> (ra, dec)
  * type == 'sky' -> (0.0, -600.0)
  * otherwise (e.g., blank/-/_S*_): (nan, nan)
- trace_row: prefer px if finite; else subpx if finite; else nan
- exclude_fiber: 1.0 for rows that are placeholders/blanks (headID '-' or starts/ends with '_' or type == 'blank' or if (ifu_x, ifu_y) are NaN); 0.0 otherwise

This script is intentionally strict about parsing but lenient about whitespace in the CSV.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Optional


def is_finite_number(value: str) -> bool:
    try:
        f = float(value)
        return math.isfinite(f)
    except Exception:
        return False


def parse_float(value: str) -> float:
    try:
        f = float(value)
        if math.isfinite(f):
            return f
    except Exception:
        pass
    return float('nan')


def should_exclude(head_id: str, typ: str, ifu_x: float, ifu_y: float) -> float:
    head = (head_id or '').strip()
    t = (typ or '').strip().lower()
    if head == '-' or head.startswith('_') and head.endswith('_'):
        return 1.0
    if t == 'blank':
        return 1.0
    if math.isnan(ifu_x) or math.isnan(ifu_y):
        return 1.0
    return 0.0


def format_value_or_nan(value: float, width: int, prec: int) -> str:
    if math.isnan(value):
        # right-align the literal 'nan' within the field width
        return f"{'nan':>{width}}"
    fmt = f"{{:{width}.{prec}f}}"
    return fmt.format(value)


def convert(input_csv: Path, output_txt: Path) -> None:
    rows = []
    with input_csv.open('r', newline='') as f:
        reader = csv.DictReader(f)
        # Normalize fieldnames (strip spaces and trailing commas)
        field_map = {name: name.strip().strip(',') for name in reader.fieldnames or []}
        reader.fieldnames = list(field_map.values())
        for raw in reader:
            row = {field_map.get(k, k): (v.strip() if isinstance(v, str) else v) for k, v in raw.items()}
            rows.append(row)

    # Prepare output directory
    output_txt.parent.mkdir(parents=True, exist_ok=True)

    # Write header similar to existing config files
    header1 = "fiber_id head_id   ifu_x    ifu_y trace_row exclude_fiber"
    header2 = "-------- ------- ------- -------- --------- -------------"

    with output_txt.open('w') as out:
        out.write(header1 + "\n")
        out.write(header2 + "\n")

        fiber_id = 1
        for r in rows:
            head_id = r.get('headID') or r.get('headId') or r.get('head_id') or ''
            typ = r.get('type', '')
            # Numbers
            px = r.get('px', '')
            subpx = r.get('subpx', '')
            ra = r.get('ra', '')
            dec = r.get('dec', '')

            # Determine values per rules
            if (typ or '').strip().lower() == 'obj':
                ifu_x = parse_float(ra)
                ifu_y = parse_float(dec)
            elif (typ or '').strip().lower() == 'sky':
                ifu_x = 0.0
                ifu_y = -600.0
            else:
                ifu_x = float('nan')
                ifu_y = float('nan')

            # trace row preference: px, else subpx, else nan
            tr = parse_float(px)
            if math.isnan(tr):
                tr = parse_float(subpx)

            exclude = should_exclude(head_id, typ, ifu_x, ifu_y)

            # Build formatted line (match widths used in existing files)
            # widths inferred from antigen/config_files/virus2/D2G/virus2_ifucen_D2G.txt
            fiber_s = f"{fiber_id:7d}"
            head_s = f"{str(head_id).strip():>7}"
            ifu_x_s = format_value_or_nan(ifu_x, 7, 4)
            ifu_y_s = format_value_or_nan(ifu_y, 8, 4)
            tr_s = format_value_or_nan(tr, 9, 1)
            excl_s = format_value_or_nan(exclude, 13, 1)

            out.write(f"{fiber_s} {head_s} {ifu_x_s} {ifu_y_s} {tr_s} {excl_s}\n")
            fiber_id += 1

    print(f"Wrote {fiber_id-1} records to {output_txt}")


def main():
    parser = argparse.ArgumentParser(description="Convert IFU center CSV to Antigen ASCII format.")
    parser.add_argument('-i', '--input', required=True, help='Path to input CSV (e.g., temp/IFU_cen_unit2G.csv)')
    parser.add_argument('-o', '--output', required=True, help='Path to output TXT (e.g., temp/virus2_ifucen_D2G.txt)')
    args = parser.parse_args()

    input_csv = Path(args.input).expanduser().resolve()
    output_txt = Path(args.output).expanduser().resolve()

    if not input_csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    convert(input_csv, output_txt)


if __name__ == '__main__':
    main()
