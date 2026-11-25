# NUMT sequences

`Calabrese_Dayama_Smart_Numts.fa` holds the NUMT sequences rtn compares reads
with (`--numts_index_dir`, `--numts_index_base`): the NUMT alleles of Dayama et
al., Calabrese et al. and Smart et al. (see RtN!'s README). It and its bwa and
samtools indexes (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`, `.fai`) are the files of
[RtN!](https://github.com/Ahhgust/RtN) as RtN! ships them (unchanged since 2019;
RtN! commit 4883c67). rtn uses them only to tag each read with its distance to
the nearest NUMT (`ZN`); which reads it removes depends on the human
mitogenomes (`../humans/`). The humans file's build uses them to leave out
genomes that look like a NUMT (`../humans/build/`).

RtN! is under the GNU General Public License, version 3; `LICENSE.txt` in this
folder is its licence text, as in RtN!'s repository. The rest of FMP-NimaGen is
under the MIT License.

Woerner AE, Cihlar JC, Smart U, Budowle B (2020). Numt identification and
removal with RtN!. Bioinformatics 36(20):5115–5116.
https://doi.org/10.1093/bioinformatics/btaa642
