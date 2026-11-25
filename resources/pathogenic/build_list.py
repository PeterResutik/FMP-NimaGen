#!/usr/bin/env python3
"""The list of pathogenic variants (--pathogenic_list) from MITOMAP's database dump: every disease variant with the
status Cfrm (confirmed), from the tables of coding/control-region variants (mmutation) and of rRNA/tRNA variants
(rtmutation). Nothing else is read from the dump.

usage: build_list.py mitomap.dump.sql.gz DATE > mitomap_confirmed_DATE.tsv

The dump: https://mitomap.org/downloads/mitomap.dump.sql.gz (DATE: its date, YYYY-MM-DD). It also holds tables that
are not needed here; delete it once the list is built.
"""
import gzip
import sys

TABLES = ("mmutation", "rtmutation")


def rows(lines):
    """(table, fields) of the rows of TABLES in a PostgreSQL dump ('COPY ... FROM stdin;' blocks)"""
    table = None
    for line in lines:
        if table is None:
            if line.startswith("COPY mitomap."):
                name = line[len("COPY mitomap."):].split(" ", 1)[0]
                table = name if name in TABLES else None
                if table:
                    columns = line[line.index("(") + 1:line.index(")")].replace('"', "").split(", ")
            continue
        if line.startswith("\\."):
            table = None
            continue
        yield table, dict(zip(columns, line.rstrip("\n").split("\t")))


def confirmed(dump_lines):
    """[(position, ref, alt, allele, locus)] of the Cfrm variants, by position; a deletion has alt '-'"""
    out = {}
    for table, r in rows(dump_lines):
        if not r["status"].startswith("Cfrm"):
            continue
        ref, alt = r["refna"], r["regna"]
        alt = "-" if alt in (":", "") else alt
        out.setdefault((int(r["position"]), ref, alt), (r["allele"], r["locus"]))
    return [(p, ref, alt, *rest) for (p, ref, alt), rest in sorted(out.items())]


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    dump, date = sys.argv[1:]
    with gzip.open(dump, "rt", encoding="utf-8", errors="replace") as h:
        variants = confirmed(h)
    print(f"# MITOMAP disease variants with the status Cfrm (confirmed), {len(variants)} variants")
    print(f"# from MITOMAP's database dump of {date} (https://mitomap.org), built with build_list.py")
    print("position\tref\talt\tallele\tlocus")
    for v in variants:
        print("\t".join(map(str, v)))


if __name__ == "__main__":
    main()
