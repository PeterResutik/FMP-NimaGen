# Pathogenic variants

`mitomap_confirmed_2026-11-18.tsv` lists the 130 disease variants that MITOMAP
gives the status Cfrm (confirmed), from its tables of coding/control-region
variants and of rRNA/tRNA variants, as of MITOMAP's database dump of
18 November 2026. The report marks a row `pathogenic = yes` when it carries one
of them (`--pathogenic_list`), so that such variants can be left out of a
haplotype, e.g. before it is shared or added to a database.

Each line gives the rCRS position, the rCRS bases there, the bases of the
variant (`-` for a deletion) and MITOMAP's name of it. The report spells each
variant by the same rules as its own rows, so a deletion in a run of equal bases
matches wherever MITOMAP placed it.

## Building the list again

```bash
curl -L -o mitomap.dump.sql.gz https://mitomap.org/downloads/mitomap.dump.sql.gz
python build_list.py mitomap.dump.sql.gz YYYY-MM-DD > mitomap_confirmed_YYYY-MM-DD.tsv
rm mitomap.dump.sql.gz
```

`YYYY-MM-DD` is the date of the dump. The script reads only the two variant
tables; the dump holds others too, so delete it afterwards. Give the new file to
the pipeline with `--pathogenic_list`.

## Source

MITOMAP: A Human Mitochondrial Genome Database. https://www.mitomap.org.
Lott MT, Leipzig JN, Derbeneva O, Xie HM, Chalkia D, Sarmady M, Procaccio V,
Wallace DC (2013). mtDNA variation and analysis using MITOMAP and MITOMASTER.
Current Protocols in Bioinformatics 44:1.23.1–26.

MITOMAP's content is shared under the Creative Commons Attribution 3.0 License
(https://creativecommons.org/licenses/by/3.0/).
