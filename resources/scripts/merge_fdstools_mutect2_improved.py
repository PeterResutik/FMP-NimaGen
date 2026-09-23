import argparse
import os
import pandas as pd
import re
import sys
import traceback
from openpyxl import load_workbook
from openpyxl.styles import PatternFill, Border, Side, Alignment

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from call_repeat_regions import MUTECT2_TARGET_REGIONS
from heteroplasmy_thresholds import lh_bounds_pct


white_border = Border(
    left=Side(border_style="thin", color="FFFFFF"),
    right=Side(border_style="thin", color="FFFFFF"),
    top=Side(border_style="thin", color="FFFFFF"),
    bottom=Side(border_style="thin", color="FFFFFF")
)


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
    instead of the same variant being reported twice.
    """
    parts = substitution_parts(variant_str)
    if parts:
        ref, pos, alt, _is_major = parts
        return f"{ref}{pos}{alt}"
    return str(variant_str).upper()


def mark_disabled_mutect2_calls(merged, reference_fasta):
    """When homopolymer_mutect2_reporting is 'disabled' (process_mutect2_
    output_improved.py has already dropped Mutect2's own indel/length-axis
    rows in MUTECT2_TARGET_REGIONS, e.g. chrM:16180-16193, with no
    replacement - see disable_homopolymer_length_calls there), mark
    called_by_MUTECT2 as "DISABLED" for the corresponding merged rows
    instead of leaving it as the plain False a genuine absence would
    produce, so a reader can tell "we didn't ask" apart from "Mutect2
    asked and found nothing" or a real DISAGREEMENT.

    Scoped to MUTECT2_TARGET_REGIONS - chrM:16180-16193 and chrM:303-315,
    matching disable_homopolymer_length_calls and process_fdstools_
    output_improved_better.py's own target_regions - both this and that
    used to be genome-wide across all 9 reference-derived boundary-run
    regions, but the FDSTOOLS side was narrowed this same session,
    leaving Mutect2's own genome-wide disabling (and this marking, if it
    stayed genome-wide too) touching regions FDSTOOLS no longer does
    anything special for at all (2026-09-22, user, chrM:6419: "I would
    restrict mutect2 as well same ways we do fdstools"; 2026-09-23,
    extended to 303-315 the same way). Within a TRUE boundary-run region
    (one with a "leading" entry, e.g. 16180-16193): any row positioned in
    the leading run, plus any length-axis row positioned in the extension
    run. A plain region with no leading run (e.g. 303-315) only has the
    length-axis-in-extension check - there's no separate leading run to
    mark. Ordinary substitution rows in the extension run (e.g. T16189C)
    are left alone either way, since Mutect2 was never disabled for those.

    Deliberately NOT keyed off a "Type" column: this runs on the already-
    merged table, where Type only ever came from Mutect2's own output (
    FDSTOOLS' has no such column) - so for exactly the rows this function
    needs to find, ones FDSTOOLS reported that Mutect2 no longer has after
    being disabled, Type is NaN, not "LHP"/"DEL"/"INS". Length-axis rows
    are identified from the FMP label's own format instead: a plain point
    substitution is always REF+POS+ALT with no decimal and no trailing
    -/lowercase suffix (e.g. "T16189C", "T16189Y"); everything else in the
    extension run - report_separate_frame's own decimal-anchored cumulative
    calls (e.g. "-16193.1c"), and shared_frame's plain-integer shifted
    deletion calls (e.g. "C16193c") alike - is a length claim.
    """
    regions = MUTECT2_TARGET_REGIONS
    if not regions:
        return merged

    def extract_position(seq):
        match = re.search(r"(\d+\.?\d*)", str(seq))
        return float(match.group(1)) if match else None

    plain_substitution_re = re.compile(r'^[ACGT]\d+[ACGTRYMKSW]$')

    positions = merged["FMP"].apply(extract_position)
    is_length_type = ~merged["FMP"].astype(str).str.match(plain_substitution_re)

    disabled_mask = pd.Series(False, index=merged.index)
    for region in regions:
        leading = region["leading"]
        extension_start, extension_end = region["extension"]
        if leading:
            leading_start, leading_end = leading
            in_leading = positions.between(leading_start, leading_end)
        else:
            in_leading = pd.Series(False, index=merged.index)
        # Not a plain between(extension_start, extension_end): a decimal
        # insertion label anchored at the extension run's own far end
        # (e.g. "-16193.1c") extracts to 16193.1, just past extension_end -
        # between() would miss it entirely.
        in_extension = (positions >= extension_start) & (positions < extension_end + 1)
        disabled_mask |= in_leading | (is_length_type & in_extension)

    merged.loc[disabled_mask, "called_by_MUTECT2"] = "DISABLED"
    return merged


def merge_variant_callers(file_fdstools: str, file_mutect2: str, lh_thresh: float = 10.0,
                          min_vf: float = 5.0, homopolymer_mutect2_reporting: str = "true_mutect2",
                          reference_fasta: str = None) -> pd.DataFrame:
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
                    avg_pct = (vf_fds + vf_mt2 * 100) / 2
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

            avg_pct = (vf_fds + vf_mt2 * 100) / 2
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

        if homopolymer_mutect2_reporting == "disabled" and reference_fasta:
            merged = mark_disabled_mutect2_calls(merged, reference_fasta)

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
        fill_disabled = PatternFill(start_color="D9D9D9", end_color="D9D9D9", fill_type="solid")  # DISABLED (caller deliberately not consulted here - not a problem, so kept distinct from red/yellow)

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

                # variant_note is now Mutect2-only (FDSTOOLS' own note
                # column was renamed to dominant_molecule_note - see
                # process_fdstools_output_improved_better.py), so the two
                # no longer collide and pd.merge no longer suffixes
                # either one with _FDSTOOLS/_MUTECT2 (2026-09-23).
                if col_name.endswith("_MUTECT2") or col_name.endswith("_MT2") or col_name in ["MUTECT2", "Filter", "Pos", "Ref", "Variant", "GT", "Type", "MBQ", "Filter", "variant_note"]:
                    cell.fill = fill_blue
                elif col_name.endswith("_FDSTOOLS") or col_name.endswith("_FDS") or col_name in ["FDSTOOLS", "interp_total", "marker", "marker_range", "num_markers", "dominant_molecule_note"]:
                    cell.fill = fill_green


                if col_name in ("called_by_FDSTOOLS", "called_by_MUTECT2"):
                    if cell.value is False:
                        cell.fill = fill_red
                    elif cell.value == "DISAGREEMENT":
                        cell.fill = fill_disagree
                    elif cell.value == "DISABLED":
                        cell.fill = fill_disabled
                    # True/False are booleans and Excel centers those by
                    # default, but "DISAGREEMENT"/"DISABLED" are plain
                    # strings, which Excel left-aligns by default - so
                    # without an explicit alignment they visually stood
                    # out from True/False in the same column. Center both
                    # axes explicitly so every value in this column lines
                    # up the same way.
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
    parser.add_argument("--homopolymer_mutect2_reporting", choices=["disabled", "true_mutect2"], default="true_mutect2",
                         help="Must match the same flag passed to process_mutect2_output_improved.py. When 'disabled', called_by_MUTECT2 is marked DISABLED (rather than the plain False an absence would produce) for rows in boundary-run regions where Mutect2's own calling was turned off - requires --reference.")
    parser.add_argument("--reference", help="Reference FASTA - required when --homopolymer_mutect2_reporting disabled")

    args = parser.parse_args()

    try:
        df_merged = merge_variant_callers(args.caller1, args.caller2, args.lh_thresh, args.min_vf,
                                           args.homopolymer_mutect2_reporting, args.reference)
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
