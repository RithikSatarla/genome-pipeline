"""Generate a small synthetic bacterial-like genome for demo purposes.

Produces a ~30 kb FASTA with realistic GC content and planted ORFs of
varying lengths on both strands. This lets the pipeline run out of the
box without needing internet access to NCBI.

Run from the project root:
    python scripts/make_demo_genome.py
"""

import random
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


# Reproducible
random.seed(42)

GENOME_LEN = 30_000
TARGET_GC = 0.50
N_GENES = 35

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "demo_genome.fasta"


def random_dna(n: int, gc: float) -> str:
    bases = []
    for _ in range(n):
        if random.random() < gc:
            bases.append(random.choice("GC"))
        else:
            bases.append(random.choice("AT"))
    return "".join(bases)


def random_orf(length_aa: int, gc: float) -> str:
    """Generate an ATG-start, stop-ending ORF of a given AA length."""
    body_len = length_aa - 1  # exclude start codon
    codons = []
    stops = {"TAA", "TAG", "TGA"}
    for _ in range(body_len):
        while True:
            c = random_dna(3, gc)
            if c not in stops:
                codons.append(c)
                break
    return "ATG" + "".join(codons) + random.choice(list(stops))


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)

    # Start with random "intergenic" DNA.
    genome = list(random_dna(GENOME_LEN, TARGET_GC))

    # Pick non-overlapping windows for genes.
    gene_specs = []
    cursor = 200
    for _ in range(N_GENES):
        length_aa = random.choice([120, 150, 180, 220, 260, 300, 350, 400])
        length_nt = length_aa * 3 + 3
        gap = random.randint(80, 400)
        start = cursor + gap
        end = start + length_nt
        if end > GENOME_LEN - 200:
            break
        strand = random.choice(["+", "-"])
        gene_specs.append((start, end, strand, length_aa))
        cursor = end

    # Plant each ORF.
    for start, end, strand, length_aa in gene_specs:
        orf = random_orf(length_aa, TARGET_GC)
        if strand == "-":
            orf = str(Seq(orf).reverse_complement())
        for i, base in enumerate(orf):
            genome[start + i] = base

    seq = "".join(genome)
    record = SeqRecord(
        Seq(seq),
        id="DEMO_SYN_001",
        description=(
            f"Synthetic demo genome | length={len(seq)} bp | "
            f"target GC={TARGET_GC:.2f} | planted_genes={len(gene_specs)} | "
            "FOR PIPELINE TESTING ONLY"
        ),
    )
    SeqIO.write([record], OUT, "fasta")
    print(f"Wrote {OUT}")
    print(f"  Length: {len(seq):,} bp")
    print(f"  Planted genes: {len(gene_specs)}")


if __name__ == "__main__":
    main()
