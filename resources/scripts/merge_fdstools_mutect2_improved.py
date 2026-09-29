import argparse
import pandas as pd
import re
import sys
import traceback
from openpyxl import load_workbook
from openpyxl.styles import PatternFill, Border, Side, Alignment


white_border = Border(
    left=Side(border_style="thin", color="FFFFFF"),
    right=Side(border_style="thin", color="FFFFFF"),
    top=Side(border_style="thin", color="FFFFFF"),
    bottom=Side(border_style="thin", color="FFFFFF")
)


def lh_bounds_pct(threshold_pct):
    """Symmetric floor/ceiling (percentage scale, 0-100) around a single
    length-heteroplasmy threshold. Accepts either side (e.g. 10 or 90) and
    always returns (floor, ceiling) with floor <= ceiling."""
    return min(threshold_pct, 100 - threshold_pct), max(threshold_pct, 100 - threshold_pct)


def is_major_format(variant_str):
    """Whether an already-formatted LHP/DEL/INS label is in major
    (uppercase- or dash-terminated) or minor (lowercase-terminated) form.
    Returns None if it can't be determined (missing/empty)."""
    if not isinstance(variant_str, str) or not variant_str.strip():
        return None
    last_char = variant_str.strip().split()[-1][-1]
    if last_char == "-":
        return True
    if last_char.isalpha():
        return last_char.isupper()
    return None


IUPAC_PAIRS = {
    "R": {"A", "G"}, "Y": {"C", "T"}, "M": {"A", "C"},
    "K": {"G", "T"}, "S": {"G", "C"}, "W": {"A", "T"},
}
_SUBSTITUTION = re.compile(r'^([ACGT])(\d+)([ACGTRYMKSW])$')
_DELETION = re.compile(r'^([ACGT])(\d+)(-|[acgt])$')


def substitution_parts(variant_str):
    """(ref, position, alt, is_major) for a point substitution label.

    A heteroplasmic call carries an IUPAC code for {ref, alt} (T16189Y),
    a homoplasmic one carries the alt base itself (T16189C). Returns None
    for anything that isn't a plain point substitution.
    """
    if not isinstance(variant_str, str):
        return None
    match = _SUBSTITUTION.match(variant_str.strip())
    if not match:
        return None
    ref, pos, code = match.groups()
    if code in IUPAC_PAIRS:
        others = IUPAC_PAIRS[code] - {ref}
        if len(others) != 1:
            return None
        return ref, pos, others.pop(), False
    return ref, pos, code, True


def merge_key(variant_str):
    """Key that puts the same locus on one row across both callers.

    Length variants only differ by letter case between a major and minor
    call ("-309.1C" vs "-309.1c"), so upper-casing is enough. Point
    substitutions don't: the heteroplasmic form is spelled with an IUPAC
    code ("T16189Y") and the homoplasmic one with the alt base
    ("T16189C"), which never match. Both are collapsed to the alt-base
    form here so the two callers meet on one row and get reconciled,
    instead of the same variant being reported twice. Deletions likewise:
    major "A523-" and minor "A523a" both key as "A523-".
    """
    parts = substitution_parts(variant_str)
    if parts:
        ref, pos, alt, _is_major = parts
        return f"{ref}{pos}{alt}"
    deletion = _DELETION.match(str(variant_str).strip())
    if deletion and deletion.group(3) in ("-", deletion.group(1).lower()):
        return f"{deletion.group(1)}{deletion.group(2)}-"
    return str(variant_str).upper()


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


def mutect2_amplicon_depths(depth_file, ranges):
    """Mutect2-side depth per amplicon: reads at the amplicon's middle position
    in the BAM Mutect2 runs on (p09's samtools depth output)."""
    depths = {}
    for line in open(depth_file):
        parts = line.split()
        if len(parts) < 3:
            continue
        pos, depth = int(parts[1]), int(float(parts[2]))
        for amplicon, (start, end) in ranges.items():
            if start <= pos <= end:
                depths[amplicon] = depth
    return depths


def build_low_rows(fds_low, mt2_depths, ranges, depth_threshold):
    """One LOW row per amplicon where either caller's depth is below
    depth_threshold. FDSTOOLS/MUTECT2 name the amplicon for the caller(s) that
    are low; called_by_* read "low" or "ok" per caller."""
    fds_depths = dict(zip(fds_low["marker"], fds_low["rd_FDS"]))
    rows = []
    for amplicon, (start, end) in ranges.items():
        fds_is_low = amplicon in fds_depths
        mt2_is_low = amplicon in mt2_depths and mt2_depths[amplicon] < depth_threshold
        if not (fds_is_low or mt2_is_low):
            continue
        rows.append({
            "FMP": "LOW",
            "FDSTOOLS": amplicon if fds_is_low else None,
            "rd_FDS": fds_depths.get(amplicon),
            "MUTECT2": amplicon if mt2_is_low else None,
            "rd_MT2": mt2_depths.get(amplicon),
            "called_by_FDSTOOLS": "low" if fds_is_low else "ok",
            "called_by_MUTECT2": "low" if mt2_is_low else "ok",
            "marker": amplicon,
            "marker_range": f"{start}-{end}",
        })
    return pd.DataFrame(rows)


def caller_average_pct(row, vf_fds, vf_mt2, weighted=False):
    """Average of the two callers' frequencies (percent). With weighted=True each
    caller is weighted by its depth for the call (FDSTOOLS: amplicon reads,
    Mutect2: ref+alt reads), falling back to the plain average if a depth is missing."""
    if not weighted:
        return (vf_fds + vf_mt2 * 100) / 2
    d_fds = pd.to_numeric(row.get("interpolated_total_coverage"), errors="coerce")
    try:
        d_mt2 = sum(float(x) for x in str(row.get("rd_MT2")).split(","))
    except ValueError:
        d_mt2 = float("nan")
    if pd.notna(d_fds) and pd.notna(d_mt2) and d_fds + d_mt2 > 0:
        return (vf_fds * d_fds + vf_mt2 * 100 * d_mt2) / (d_fds + d_mt2)
    return (vf_fds + vf_mt2 * 100) / 2


def merge_variant_callers(file_fdstools: str, file_mutect2: str, lh_thresh: float = 10.0,
                          min_vf: float = 5.0, mutect2_depth_file: str = None,
                          marker_map: str = None, depth_threshold: int = 10,
                          weighted_average: bool = False) -> pd.DataFrame:
    try:
        df1 = pd.read_csv(file_fdstools, sep="\t")
        df2 = pd.read_csv(file_mutect2, sep="\t")
    except Exception as e:
        print(f"Error reading input files: {e}", file=sys.stderr)
        raise

    if df1.empty and df2.empty:
        print("Both inputs are empty (negative control / H2O). Writing empty output.", file=sys.stderr)
        return pd.DataFrame()

    try:
        # LOW rows are rebuilt per amplicon for both callers after the merge
        fds_low = df1[df1["FDSTOOLS"] == "LOW"]
        df1 = df1[df1["FDSTOOLS"] != "LOW"].copy()
        if marker_map:
            ranges = load_amplicon_ranges(marker_map)
        else:
            ranges = {m: tuple(map(int, str(r).split("-"))) for m, r in zip(fds_low["marker"], fds_low["marker_range"])}
        mt2_depths = mutect2_amplicon_depths(mutect2_depth_file, ranges) if mutect2_depth_file else {}

        # Rename original variant columns before merge to avoid conflict
        df1["FMP"] = df1["FDSTOOLS"]
        df2["FMP"] = df2["MUTECT2"]

        # Join on a case-normalized key so the same locus is recognized as
        # the same event even when the two callers' own individual variant
        # frequencies land on opposite sides of the major/minor case
        # boundary for a length-heteroplasmy call (e.g. "-309.1C" vs
        # "-309.1c"). Each caller's own original-case label is preserved
        # in FMP_FDSTOOLS/FMP_MUTECT2 below and reconciled afterward.
        df1["_merge_key"] = df1["FMP"].apply(merge_key)
        df2["_merge_key"] = df2["FMP"].apply(merge_key)

        merged = pd.merge(df1, df2, on="_merge_key", how="outer", suffixes=("_FDSTOOLS", "_MUTECT2"))
        merged.drop(columns=["_merge_key"], inplace=True)

        # Where both callers agree on a length-heteroplasmy locus but
        # disagree on major-vs-minor classification (case mismatch above),
        # decide from the *averaged* variant frequency instead of trusting
        # either caller's individual estimate, and use whichever caller's
        # already-shifted/formatted label matches that averaged decision.
        # fds_confirms/mt2_confirms record whether each caller's *own* call
        # actually matches the reported classification (None when no
        # reconciliation was needed/possible, i.e. not a both-called length
        # variant) - used below to correct called_by_FDSTOOLS/MUTECT2 so
        # they mean "this caller supports what's reported", not just
        # "this caller called something at this locus".
        lh_floor, lh_ceiling = lh_bounds_pct(lh_thresh)

        def reconcile_row(row):
            fmp_fds = row.get("FMP_FDSTOOLS")
            fmp_mt2 = row.get("FMP_MUTECT2")
            current_type = row.get("Type")
            both_called = pd.notna(fmp_fds) and pd.notna(fmp_mt2)
            is_length_type = current_type in ("DEL", "INS", "LHP")

            # Point substitutions get the same treatment as length variants:
            # when both callers hit the locus but land on opposite sides of
            # the homoplasmy boundary - one spelling it T16189C, the other
            # T16189Y - decide from the averaged frequency rather than
            # trusting either estimate, and let called_by_* flag whichever
            # caller disagrees with what gets reported. Without this the two
            # spellings never matched, so the same variant appeared twice in
            # the report with no indication the callers disagreed.
            if both_called and not is_length_type:
                fds_parts = substitution_parts(fmp_fds)
                mt2_parts = substitution_parts(fmp_mt2)
                vf_fds, vf_mt2 = row.get("vf_FDS"), row.get("vf_MT2")
                if (fds_parts and mt2_parts
                        and fds_parts[3] != mt2_parts[3]
                        and pd.notna(vf_fds) and pd.notna(vf_mt2)):
                    ref, pos, alt, _ = fds_parts
                    avg_pct = caller_average_pct(row, vf_fds, vf_mt2, weighted_average)
                    desired_major = avg_pct >= (100 - min_vf)
                    if desired_major:
                        fmp = f"{ref}{pos}{alt}"
                    else:
                        code = next((c for c, bases in IUPAC_PAIRS.items()
                                     if bases == {ref, alt}), alt)
                        fmp = f"{ref}{pos}{code}"
                    return pd.Series({
                        "FMP": fmp,
                        "Type": "SNP" if desired_major else "PHP",
                        "fds_confirms": fds_parts[3] == desired_major,
                        "mt2_confirms": mt2_parts[3] == desired_major,
                    })

            if not (both_called and is_length_type):
                fmp = fmp_fds if pd.notna(fmp_fds) else fmp_mt2
                return pd.Series({"FMP": fmp, "Type": current_type,
                                   "fds_confirms": None, "mt2_confirms": None})

            vf_fds, vf_mt2 = row.get("vf_FDS"), row.get("vf_MT2")
            if pd.isna(vf_fds) or pd.isna(vf_mt2):
                fmp = fmp_mt2 if pd.notna(fmp_mt2) else fmp_fds
                return pd.Series({"FMP": fmp, "Type": current_type,
                                   "fds_confirms": None, "mt2_confirms": None})

            avg_pct = caller_average_pct(row, vf_fds, vf_mt2, weighted_average)
            desired_major = avg_pct >= lh_ceiling

            fds_confirms = is_major_format(fmp_fds) == desired_major
            mt2_confirms = is_major_format(fmp_mt2) == desired_major

            if fds_confirms:
                fmp = fmp_fds
            elif mt2_confirms:
                fmp = fmp_mt2
            else:
                fmp = fmp_mt2  # fallback; shouldn't normally happen

            if desired_major == (current_type in ("DEL", "INS")):
                var_type = current_type
            elif desired_major:
                ref, var = str(row.get("Ref", "")), str(row.get("Variant", ""))
                var_type = "DEL" if len(ref) > len(var) else "INS"
            else:
                var_type = "LHP"

            return pd.Series({"FMP": fmp, "Type": var_type,
                               "fds_confirms": fds_confirms, "mt2_confirms": mt2_confirms})

        reconciled = merged.apply(reconcile_row, axis=1)
        merged["FMP"] = reconciled["FMP"]
        if "Type" in merged.columns:
            merged["Type"] = reconciled["Type"]
        merged.drop(columns=["FMP_FDSTOOLS", "FMP_MUTECT2"], inplace=True, errors="ignore")

        merged["called_by_FDSTOOLS"] = (~merged["FDSTOOLS"].isna()).astype(bool)
        merged["called_by_MUTECT2"] = (~merged["MUTECT2"].isna()).astype(bool)

        # fds_confirms/mt2_confirms are only non-null when BOTH callers made
        # a call at this locus (both_called above) - so a False here never
        # means "this caller didn't call it" (that's the plain False set
        # just above, untouched by this override). It means the caller DID
        # call it, but on the other side of the major/minor (or SNP/PHP)
        # threshold from the averaged-frequency decision that FMP reports -
        # a real disagreement between the two callers, not a non-call.
        # Flagging both cases identically as False reads as "MUTECT2/
        # FDSTOOLS missed this", which is misleading for the disagreement
        # case, so it gets its own label instead of False.
        merged["called_by_FDSTOOLS"] = merged["called_by_FDSTOOLS"].astype(object)
        merged["called_by_MUTECT2"] = merged["called_by_MUTECT2"].astype(object)

        fds_override = reconciled["fds_confirms"].notna()
        merged.loc[fds_override.values, "called_by_FDSTOOLS"] = reconciled.loc[fds_override, "fds_confirms"].apply(
            lambda confirms: True if confirms else "DISAGREEMENT").values
        mt2_override = reconciled["mt2_confirms"].notna()
        merged.loc[mt2_override.values, "called_by_MUTECT2"] = reconciled.loc[mt2_override, "mt2_confirms"].apply(
            lambda confirms: True if confirms else "DISAGREEMENT").values

        def extract_position(seq):
            match = re.search(r"(\d+\.?\d*)", str(seq))
            return float(match.group(1)) if match else float('inf')

        merged["variant_position"] = merged["FMP"].apply(extract_position)
        merged = merged.sort_values(by="variant_position").drop(columns=["variant_position"])
        # merged = merged.sort_values(by="marker")

        front = ["FMP"]
        # other = [col for col in merged.columns if col not in front]
        # return merged[front + other]

        # Reorder known columns if they exist
        priority_order = [
            "FDSTOOLS",
            "vf_FDS",
            "rd_FDS",
            "MUTECT2",
            "vf_MT2",
            "rd_MT2",
            "MBQ",
            "called_by_FDSTOOLS",
            "called_by_MUTECT2"
        ]
        existing_priority = [col for col in priority_order if col in merged.columns]
        remaining_columns = [col for col in merged.columns if col not in front and col not in existing_priority]
        merged = merged[front + existing_priority + remaining_columns]

        low_rows = build_low_rows(fds_low, mt2_depths, ranges, depth_threshold)
        if not low_rows.empty:
            merged = pd.concat([merged, low_rows], ignore_index=True)

        return merged


    except Exception as e:
        print(f"Error processing and merging data: {e}", file=sys.stderr)
        raise

def apply_excel_styles(excel_path: str):
    try:
        wb = load_workbook(excel_path)
        ws = wb.active

        # Define fills
        fill_green = PatternFill(start_color="C6EFCE", end_color="C6EFCE", fill_type="solid")  # D to M
        fill_blue = PatternFill(start_color="D9E1F2", end_color="D9E1F2", fill_type="solid")   # N to X
        fill_red = PatternFill(start_color="FFC7CE", end_color="FFC7CE", fill_type="solid")    # False flags (caller didn't call it at all)
        fill_disagree = PatternFill(start_color="FFD966", end_color="FFD966", fill_type="solid")  # DISAGREEMENT flags (caller called it, but on the other side of the major/minor threshold)

        # Define additional fill colors for column A
        fill_low = PatternFill(start_color="FEFE01", end_color="FEFE01", fill_type="solid")
        fill_iupac = PatternFill(start_color="FAC000", end_color="FAC000", fill_type="solid")
        fill_lowercase = PatternFill(start_color="3CB0F1", end_color="3CB0F1", fill_type="solid")
        fill_false_flag = PatternFill(start_color="F50003", end_color="F50003", fill_type="solid")
        fill_dash = PatternFill(start_color="92D14F", end_color="92D14F", fill_type="solid")
        fill_default = PatternFill(start_color="36B150", end_color="36B150", fill_type="solid")


        header = [cell.value for cell in ws[1]]
        for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
            for idx, cell in enumerate(row):
                col_name = header[idx]
                cell.border = white_border

                # Special coloring for first column "variant"
                if idx == 0:
                    val = str(cell.value)
                    if val == "LOW":
                        cell.fill = fill_low
                    elif row[8].value is False or row[9].value is False or row[8].value == "DISAGREEMENT" or row[9].value == "DISAGREEMENT":
                        cell.fill = fill_false_flag
                    elif re.search(r"[a-z]", val):
                        cell.fill = fill_lowercase
                    elif "-" in val:
                        cell.fill = fill_dash
                    elif re.search(r"[MRYWSK]", val):
                        cell.fill = fill_iupac
                    else:
                        cell.fill = fill_default


                # Red fill for False flags
                if col_name == "called_by_FDSTOOLS" and cell.value is False:
                    cell.fill = fill_green
                elif col_name == "called_by_MUTECT2" and cell.value is False:
                    cell.fill = fill_red

                if col_name.endswith("_MUTECT2") or col_name.endswith("_MT2") or col_name in ["MUTECT2", "Filter", "Pos", "Ref", "Variant", "GT", "Type", "MBQ", "Filter" ]:
                    cell.fill = fill_blue
                elif col_name.endswith("_FDSTOOLS") or col_name.endswith("_FDS") or col_name in ["FDSTOOLS", "interp_total", "marker", "marker_range", "num_markers", "variant_note"]:
                    cell.fill = fill_green


                if col_name in ("called_by_FDSTOOLS", "called_by_MUTECT2"):
                    # In a LOW row these say "low" or "ok" for each caller's depth
                    if row[0].value == "LOW":
                        if cell.value == "low":
                            cell.fill = fill_low
                    elif cell.value is False:
                        cell.fill = fill_red
                    elif cell.value == "DISAGREEMENT":
                        cell.fill = fill_disagree
                    # True/False are booleans and Excel centers those by
                    # default, but "DISAGREEMENT" is a plain string, which
                    # Excel left-aligns by default - so without an explicit
                    # alignment it visually stood out from True/False in
                    # the same column. Center both axes explicitly so all
                    # three values in this column line up the same way.
                    cell.alignment = Alignment(horizontal="center", vertical="center")
        for col in ws.columns:
            max_length = 0
            column = col[0].column_letter  # Get column name (e.g., 'A', 'B', etc.)
            for cell in col:
                try:
                    if cell.value:
                        max_length = max(max_length, len(str(cell.value)))
                except:
                    pass
            adjusted_width = (max_length + 1)
            ws.column_dimensions[column].width = adjusted_width

        wb.save(excel_path)
    except Exception as e:
        print(f"Error styling Excel file: {e}", file=sys.stderr)
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge mitochondrial variant calls from FDSTOOLS and MUTECT2.")
    parser.add_argument("caller1", help="Path to the FDSTOOLS file (TSV format, with 'FDSTOOLS' column).")
    parser.add_argument("caller2", help="Path to the MUTECT2 file (TSV format, with 'MUTECT2' column).")
    parser.add_argument("output_file", help="Path to save the merged output (XLSX format).")
    parser.add_argument("--lh_thresh", type=float, default=10.0, help="Symmetric length-heteroplasmy threshold, used to reconcile major/minor case disagreements between callers via their averaged variant frequency")
    parser.add_argument("--min_vf", type=float, default=5.0, help="Minor allele frequency threshold; its complement (100-min_vf) is the homoplasmy boundary used to reconcile SNP-vs-IUPAC disagreements between callers via their averaged variant frequency")
    parser.add_argument("--mutect2_depth", help="samtools depth at each amplicon's middle position in the BAM Mutect2 runs on (p09); amplicons below --depth are reported as LOW for Mutect2")
    parser.add_argument("--marker_map", help="FDSTOOLS library file, for amplicon names and ranges")
    parser.add_argument("--depth", type=int, default=10, help="Read depth below which an amplicon is reported as LOW")
    parser.add_argument("--disagreement_average", choices=["plain", "depth_weighted"], default="plain", help="How the two callers' frequencies are averaged when they disagree on major vs minor: plain average, or weighted by each caller's read depth")

    args = parser.parse_args()

    try:
        df_merged = merge_variant_callers(args.caller1, args.caller2, args.lh_thresh, args.min_vf,
                                          args.mutect2_depth, args.marker_map, args.depth,
                                          args.disagreement_average == "depth_weighted")
        if df_merged.empty:
            pd.DataFrame([["No variants detected (negative control / H2O)"]]).to_excel(args.output_file, index=False, header=False, engine="openpyxl")
            print(f"Empty output written to: {args.output_file}")
            sys.exit(0)
        df_merged.drop(columns=["ID","is_noise_or_low_frq"], errors="ignore", inplace=True)
        df_merged.rename(columns={
            "interpolated_total_coverage": "interp_total",
            "total_wo_noise_or_low_frq": "total_clean",
            "variant_frequency_wo_noise_or_low_frq": "variant_frequency_clean"
        }, inplace=True)
        if "vf_MT2" in df_merged.columns:
            df_merged["vf_MT2"] = (df_merged["vf_MT2"] * 100).round(2)
        df_merged.to_excel(args.output_file, index=False, engine="openpyxl")
        print(f"Merged Excel file saved to: {args.output_file}")
    except Exception:
        traceback.print_exc()
        print("Merging failed. Please check your input files and formats.", file=sys.stderr)
        sys.exit(1)

    try:
        apply_excel_styles(args.output_file)
        print("Excel styling completed successfully.")
    except Exception:
        traceback.print_exc()
        print("Styling failed. The Excel file was created, but formatting could not be applied.", file=sys.stderr)
        sys.exit(1)
