"""Generate plots for the genome analysis report."""

from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")  # non-interactive backend


def plot_gc_content(positions, gc_values, mean_gc, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 3.5))
    ax.plot(positions, gc_values, color="#2563eb", lw=1.2)
    ax.axhline(mean_gc, color="#ef4444", ls="--", lw=1, label=f"Mean GC = {mean_gc:.2%}")
    ax.set_xlabel("Genome position (bp)")
    ax.set_ylabel("GC content")
    ax.set_title("GC content across the genome (1 kb sliding window)")
    ax.legend(loc="best", frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_gc_skew(positions, skew_values, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 3.5))
    ax.plot(positions, skew_values, color="#16a34a", lw=1.2)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xlabel("Genome position (bp)")
    ax.set_ylabel("(G−C) / (G+C)")
    ax.set_title("GC skew along the genome")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_orf_length_distribution(orfs, out_path: Path) -> None:
    lengths = [o.length_aa for o in orfs]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(lengths, bins=30, color="#7c3aed", edgecolor="white")
    ax.set_xlabel("ORF length (amino acids)")
    ax.set_ylabel("Count")
    ax.set_title(f"Predicted ORF length distribution (n = {len(orfs)})")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_orf_map(orfs, genome_length: int, out_path: Path) -> None:
    """Linear map of ORF positions, colored by strand."""
    fig, ax = plt.subplots(figsize=(11, 2.5))
    for orf in orfs:
        y = 0.6 if orf.strand == "+" else 0.2
        color = "#2563eb" if orf.strand == "+" else "#ef4444"
        ax.barh(y, orf.end - orf.start, left=orf.start, height=0.25,
                color=color, edgecolor="white", linewidth=0.4)
    ax.set_xlim(0, genome_length)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.2, 0.6])
    ax.set_yticklabels(["− strand", "+ strand"])
    ax.set_xlabel("Genome position (bp)")
    ax.set_title(f"ORF map ({len(orfs)} ORFs)")
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_aa_composition(comp: dict, out_path: Path) -> None:
    aas = list(comp.keys())
    freqs = [comp[a] for a in aas]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(aas, freqs, color="#0ea5e9", edgecolor="white")
    ax.set_xlabel("Amino acid")
    ax.set_ylabel("Frequency")
    ax.set_title("Amino acid composition of predicted proteome")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
