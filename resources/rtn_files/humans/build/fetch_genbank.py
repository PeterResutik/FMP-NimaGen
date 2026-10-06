#!/usr/bin/env python3
"""Download GenBank sequences by accession (NCBI E-utilities, FASTA), 200 per request,
at most two requests a second. Finished batches are kept, so the script can be
started again after an interruption.

usage: fetch_genbank.py accessions.txt outdir
"""
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ACC = [l.strip() for l in open(sys.argv[1]) if l.strip()]
OUT = Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
B = 200
failed = 0
for i in range(0, len(ACC), B):
    f = OUT / f"batch_{i // B:04d}.fa"
    if f.exists() and f.stat().st_size > 0:
        continue
    ids = ACC[i:i + B]
    data = urllib.parse.urlencode({"db": "nuccore", "id": ",".join(ids), "rettype": "fasta",
                                   "retmode": "text", "tool": "fmp_nimagen_humans_build"}).encode()
    for attempt in range(6):
        try:
            with urllib.request.urlopen(urllib.request.Request(URL, data=data), timeout=120) as r:
                txt = r.read().decode()
            n = txt.count(">")
            if n >= 0.9 * len(ids):
                f.write_text(txt)
                break
            raise ValueError(f"only {n} of {len(ids)} records")
        except Exception as e:
            print(f"batch {i // B}: {e}; retry {attempt + 1}", flush=True)
            time.sleep(5 * (attempt + 1))
    else:
        print(f"batch {i // B}: FAILED", flush=True)
        failed += 1
    time.sleep(0.5)
    if (i // B) % 25 == 0:
        print(f"batch {i // B} of {len(ACC) // B}", flush=True)
sys.exit(f"{failed} batches failed; start the script again" if failed else 0)
