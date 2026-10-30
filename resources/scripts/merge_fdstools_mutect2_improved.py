"""The report: FDSTOOLS' calls (the table written by process_fdstools_output_improved_better.py)
with one row per variant and one LOW row per amplicon below --depth, written as a coloured Excel file."""
import argparse
import pandas as pd
import re
import sys
import traceback
from openpyxl import load_workbook
from openpyxl.styles import PatternFill, Border, Side

import iupac


white_border = Border(
    left=Side(border_style="thin", color="FFFFFF"),
    right=Side(border_style="thin", color="FFFFFF"),
    top=Side(border_style="thin", color="FFFFFF"),
    bottom=Side(border_style="thin", color="FFFFFF")
)


_SUBSTITUTION = re.compile(r'^([ACGT])(\d+)([ACGTRYMKSWVHDBN])$')
_DELETION = re.compile(r'^([ACGT])(\d+)(-|[acgt])$')


def substitution_parts(variant_str):
    """(ref, position, code) of a point substitution label: the base when it is the
    only one present (T16189C), else the IUPAC code of all bases present (T16189Y,
    T16189H). None for anything else."""
    if not isinstance(variant_str, str):
        return None
    match = _SUBSTITUTION.match(variant_str.strip())
    if not match or match.group(3) == match.group(1):
        return None
    return match.groups()


def order_key(variant_str):
    """Key that orders the rows before they are sorted by position: a deletion at a
    position (A523- or A523a, keyed "A523-") comes before the substitution there
    (T16189C, T16189Y, ... keyed "T16189~", as "~" sorts after "-"); other labels
    are keyed upper-cased."""
    parts = substitution_parts(variant_str)
    if parts:
        ref, pos, _code = parts
        return f"{ref}{pos}~"
    deletion = _DELETION.match(str(variant_str).strip())
    if deletion and deletion.group(3) in ("-", deletion.group(1).lower()):
        return f"{deletion.group(1)}{deletion.group(2)}-"
    return str(variant_str).upper()


def as_number(value):
    """A value from the table as a number where it is one; values per base
    ("C 40, A 10") stay text."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def one_row_per_position(df, min_vf, deleted):
    """FDSTOOLS' substitution rows at one position as one row: T16189Y (C 40%) and
    T16189W (A 10%) become T16189H with vf "C 40, A 10" and reads "C 200, A 50".
    rCRS counts as present when what is left of the reads (100% minus the other
    bases and the molecules with a deletion there, `deleted`) reaches min_vf, so
    C 70% and A 30% become T16189M. The FDSTOOLS step already writes one row per
    position, so this only joins rows it left apart."""
    parts = df["FDSTOOLS"].apply(substitution_parts)
    groups = {}
    for i, p in parts.items():
        if p:
            groups.setdefault(p[:2], []).append(i)
    rows, combined = [], []
    for (ref, pos), idx in groups.items():
        if len(idx) < 2:
            continue
        vfs, reads = {}, {}
        for i in idx:
            bases = iupac.alts(parts[i][2], ref)
            values = iupac.parse_values(df.at[i, "vf_FDS"], bases)
            n = as_number(df.at[i, "rd_FDS"])
            if values is None or not isinstance(n, float):  # already given per base: left as it is
                break
            vfs.update(values)
            reads[next(iter(bases))] = n
        else:
            ref_share = 100 - sum(vfs.values()) - deleted.get((ref, pos), 0)
            present = set(vfs) | ({ref} if ref_share >= min_vf else set())
            row = {c: " | ".join(dict.fromkeys(str(v) for v in df.loc[idx, c].dropna())) or None
                   for c in df.columns}
            order = sorted(vfs, key=lambda b: (-vfs[b], b))
            row["FDSTOOLS"] = f"{ref}{pos}{iupac.code(present)}"
            row["vf_FDS"] = iupac.format_values(vfs, 2)
            row["rd_FDS"] = ", ".join(f"{b} {reads[b]:g}" for b in order)
            rows.append(row)
            combined += idx
    if not rows:
        return df
    return pd.concat([df.drop(index=combined), pd.DataFrame(rows)], ignore_index=True)


def deletion_shares(labels, vfs):
    """{(ref, position): percent} of FDSTOOLS' deletion rows (A523-, A523a)."""
    out = {}
    for label, vf in zip(labels, vfs):
        m = _DELETION.match(str(label))
        if m and m.group(3) in ("-", m.group(1).lower()):
            value = as_number(vf)
            if isinstance(value, float):
                out[(m.group(1), m.group(2))] = value
    return out


def load_amplicon_ranges(library_path):
    """{amplicon: (start, end)} from the FDSTOOLS library's [genome_position] block."""
    ranges, in_block = {}, False
    for line in open(library_path):
        line = line.strip()
        if line.startswith("[genome_position]"):
            in_block = True
            continue
        if in_block and line.startswith("["):
            break
        if in_block and "=" in line:
            name, values = line.split("=", 1)
            parts = [v.strip() for v in values.split(",")]
            if len(parts) >= 3:
                ranges[name.strip()] = (int(parts[1]), int(parts[2]))
    return ranges


def build_low_rows(fds_low, ranges):
    """One LOW row per amplicon below the depth threshold, in the library's order,
    with FDSTOOLS' reads of the amplicon."""
    depths = dict(zip(fds_low["marker"], fds_low["rd_FDS"]))
    return pd.DataFrame([{"FMP": "LOW", "rd_FDS": depths[amplicon], "marker": amplicon,
                          "marker_range": f"{start}-{end}"}
                         for amplicon, (start, end) in ranges.items() if amplicon in depths])


def fdstools_report(file_fdstools: str, min_vf: float = 5.0, marker_map: str = None) -> pd.DataFrame:
    df = pd.read_csv(file_fdstools, sep="\t")
    if df.empty:
        print("The FDSTOOLS table is empty (negative control / H2O). Writing empty output.", file=sys.stderr)
        return pd.DataFrame()

    for c in ("vf_FDS", "rd_FDS"):
        if c in df.columns and df[c].dtype == object:
            df[c] = df[c].map(as_number)

    # LOW rows are rebuilt per amplicon at the end
    fds_low = df[df["FDSTOOLS"] == "LOW"]
    df = df[df["FDSTOOLS"] != "LOW"].copy()
    if marker_map:
        ranges = load_amplicon_ranges(marker_map)
    else:
        ranges = {m: tuple(map(int, str(r).split("-"))) for m, r in zip(fds_low["marker"], fds_low["marker_range"])}

    df = one_row_per_position(df, min_vf, deletion_shares(df["FDSTOOLS"], df["vf_FDS"]))
    df = df.rename(columns={"FDSTOOLS": "FMP"})

    # by position; at one position a deletion before a substitution
    df["_order_key"] = df["FMP"].apply(order_key)
    df = df.sort_values(by="_order_key", kind="mergesort").drop(columns=["_order_key"]).reset_index(drop=True)
    df["variant_position"] = df["FMP"].apply(
        lambda seq: float(m.group(1)) if (m := re.search(r"(\d+\.?\d*)", str(seq))) else float("inf"))
    df = df.sort_values(by="variant_position").drop(columns=["variant_position"])

    front = [c for c in ("FMP", "vf_FDS", "rd_FDS") if c in df.columns]
    df = df[front + [c for c in df.columns if c not in front]]

    low_rows = build_low_rows(fds_low, ranges)
    if not low_rows.empty:
        df = pd.concat([df, low_rows], ignore_index=True)
    return df


def apply_excel_styles(excel_path: str):
    try:
        wb = load_workbook(excel_path)
        ws = wb.active

        fill_green = PatternFill(start_color="C6EFCE", end_color="C6EFCE", fill_type="solid")  # FDSTOOLS' columns

        # The first column (the call) by kind
        fill_low = PatternFill(start_color="FEFE01", end_color="FEFE01", fill_type="solid")
        fill_iupac = PatternFill(start_color="FAC000", end_color="FAC000", fill_type="solid")
        fill_lowercase = PatternFill(start_color="3CB0F1", end_color="3CB0F1", fill_type="solid")
        fill_dash = PatternFill(start_color="92D14F", end_color="92D14F", fill_type="solid")
        fill_default = PatternFill(start_color="36B150", end_color="36B150", fill_type="solid")

        header = [cell.value for cell in ws[1]]
        for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
            for idx, cell in enumerate(row):
                col_name = header[idx]
                cell.border = white_border

                if idx == 0:
                    val = str(cell.value)
                    if val == "LOW":
                        cell.fill = fill_low
                    elif re.search(r"[a-z]", val):
                        cell.fill = fill_lowercase
                    elif "-" in val:
                        cell.fill = fill_dash
                    elif re.match(r"^[ACGT]\d+[RYMKSWVHDBN]$", val):
                        cell.fill = fill_iupac
                    else:
                        cell.fill = fill_default
                elif col_name.endswith("_FDS") or col_name in ["marker", "marker_range", "num_markers", "variant_note", "dominant_molecule"]:
                    cell.fill = fill_green

        for col in ws.columns:
            max_length = 0
            column = col[0].column_letter  # Get column name (e.g., 'A', 'B', etc.)
            for cell in col:
                if cell.value:
                    max_length = max(max_length, len(str(cell.value)))
            ws.column_dimensions[column].width = max_length + 1

        wb.save(excel_path)
    except Exception as e:
        print(f"Error styling Excel file: {e}", file=sys.stderr)
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write the report from FDSTOOLS' calls.")
    parser.add_argument("fdstools_table", help="FDSTOOLS' calls (TSV, with an 'FDSTOOLS' column), from process_fdstools_output_improved_better.py")
    parser.add_argument("output_file", help="Path to save the report (XLSX format).")
    parser.add_argument("--min_vf", type=float, default=5.0, help="Minor allele frequency threshold: where FDSTOOLS writes a position's bases on separate rows, rCRS counts as present when its share reaches min_vf")
    parser.add_argument("--marker_map", help="FDSTOOLS library file, for amplicon names and ranges")
    parser.add_argument("--variant_note", choices=["true", "false"], default="true", help="false: the report has no variant_note column")

    args = parser.parse_args()

    try:
        report = fdstools_report(args.fdstools_table, args.min_vf, args.marker_map)
        if report.empty:
            pd.DataFrame([["No variants detected (negative control / H2O)"]]).to_excel(args.output_file, index=False, header=False, engine="openpyxl")
            print(f"Empty output written to: {args.output_file}")
            sys.exit(0)
        report.drop(columns=["ID", "is_noise_or_low_frq"], errors="ignore", inplace=True)
        if args.variant_note == "false":
            report = report.drop(columns=["variant_note"], errors="ignore")
        report.to_excel(args.output_file, index=False, engine="openpyxl")
        print(f"Report saved to: {args.output_file}")
    except Exception:
        traceback.print_exc()
        print("Writing the report failed. Please check the input file and its format.", file=sys.stderr)
        sys.exit(1)

    try:
        apply_excel_styles(args.output_file)
        print("Excel styling completed successfully.")
    except Exception:
        traceback.print_exc()
        print("Styling failed. The Excel file was created, but formatting could not be applied.", file=sys.stderr)
        sys.exit(1)
