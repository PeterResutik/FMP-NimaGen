"""mix_reads.py and evaluate_mixtures.py: mixtures of real samples."""
import gzip
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "evaluation"))
import evaluate_mixtures as em  # noqa: E402
import mix_reads as mr  # noqa: E402


def write_sample(path, reads):
    """reads: [(R1 start, number of pairs)]"""
    with gzip.open(f"{path}_R1_001.fastq.gz", "wt") as f1, gzip.open(f"{path}_R2_001.fastq.gz", "wt") as f2:
        n = 0
        for start, count in reads:
            for _ in range(count):
                n += 1
                seq = (start + "ACGT" * 40)[:100]
                f1.write(f"@r{n} 1:N:0:1\n{seq}\n+\n{'I' * len(seq)}\n")
                f2.write(f"@r{n} 2:N:0:1\n{seq}\n+\n{'I' * len(seq)}\n")


def test_every_amplicon_gets_the_exact_share(tmp_path, monkeypatch):
    index = mr.primer_index(mr.REPO / "resources/primers/left_primers.fasta")
    x = next(k for k, v in index.items() if v == "010")
    y = next(k for k, v in index.items() if v == "020")
    write_sample(tmp_path / "A", [(x, 50), (y, 30), ("N" * 35, 10)])
    write_sample(tmp_path / "B", [(x, 40), (y, 5), ("N" * 35, 10)])
    monkeypatch.setattr(sys, "argv", ["mix_reads.py", str(tmp_path / "A_R1_001.fastq.gz"),
                                      str(tmp_path / "B_R1_001.fastq.gz"), "20", str(tmp_path / "mix"), "1"])
    mr.main()
    taken = {l.split("\t")[0]: l.split("\t")[1:] for l in open(tmp_path / "mix_mix.tsv").read().splitlines()[1:]}
    # 010: A limits (50 / 0.8), 020: B limits (5 / 0.2), no primer: A limits (10 / 0.8)
    assert taken["010"] == ["50", "50", "12", "40"]
    assert taken["020"] == ["20", "30", "5", "5"]
    assert taken["no_primer"] == ["10", "10", "2", "10"]
    with gzip.open(tmp_path / "mix_R1_001.fastq.gz", "rt") as f:
        tags = Counter(l[1:6] for i, l in enumerate(f) if i % 4 == 0)
    assert tags == {"mixA_": 80, "mixB_": 19}


def test_a_mixture_name_gives_both_samples_and_the_share():
    m = em.MIXTURE.fullmatch("mix-HG03295-HG00117-m5_L001")
    assert (m[1], m[2], m[3]) == ("HG03295", "HG00117", "5")
    assert em.MIXTURE.fullmatch("mix-HG00376-HG01799-B02_L001")[3] == "02"
