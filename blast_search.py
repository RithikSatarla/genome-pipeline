"""Compute basic statistics on a genome sequence."""

from collections import Counter
from dataclasses import dataclass


@dataclass
class GenomeStats:
    length: int
    gc_content: float
    base_counts: dict
    n_count: int

    def summary(self) -> str:
        lines = [
            f"  Length:       {self.length:,} bp",
            f"  GC content:   {self.gc_content:.2%}",
            f"  A: {self.base_counts.get('A', 0):,}  "
            f"T: {self.base_counts.get('T', 0):,}  "
            f"G: {self.base_counts.get('G', 0):,}  "
            f"C: {self.base_counts.get('C', 0):,}",
            f"  Ambiguous (N): {self.n_count:,}",
        ]
        return "\n".join(lines)


def genome_stats(seq: str) -> GenomeStats:
    """Compute length, GC content, and base counts."""
    seq = str(seq).upper()
    counts = Counter(seq)
    g, c = counts.get("G", 0), counts.get("C", 0)
    a, t = counts.get("A", 0), counts.get("T", 0)
    gc = (g + c) / max(g + c + a + t, 1)
    return GenomeStats(
        length=len(seq),
        gc_content=gc,
        base_counts=dict(counts),
        n_count=counts.get("N", 0),
    )


def gc_content_window(seq: str, window: int = 1000, step: int = 200) -> tuple[list, list]:
    """Sliding-window GC content along the genome.

    Returns (positions, gc_values) for plotting.
    """
    seq = str(seq).upper()
    positions, gcs = [], []
    for i in range(0, len(seq) - window + 1, step):
        sub = seq[i : i + window]
        g = sub.count("G")
        c = sub.count("C")
        a = sub.count("A")
        t = sub.count("T")
        denom = max(g + c + a + t, 1)
        positions.append(i + window // 2)
        gcs.append((g + c) / denom)
    return positions, gcs


def gc_skew(seq: str, window: int = 1000, step: int = 200) -> tuple[list, list]:
    """Compute (G-C)/(G+C) skew in sliding windows.

    GC skew flips sign at the origin and terminus of replication
    in many bacterial genomes — a classic bioinformatics signature.
    """
    seq = str(seq).upper()
    positions, skews = [], []
    for i in range(0, len(seq) - window + 1, step):
        sub = seq[i : i + window]
        g = sub.count("G")
        c = sub.count("C")
        if g + c == 0:
            continue
        positions.append(i + window // 2)
        skews.append((g - c) / (g + c))
    return positions, skews
