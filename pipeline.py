"""Find open reading frames (ORFs) in a genome.

An ORF is a stretch of DNA that starts with ATG and ends with a stop codon
(TAA/TAG/TGA), with no stop codons in between. ORFs above a length threshold
are good candidate genes.
"""

from dataclasses import dataclass
from Bio.Seq import Seq

START = "ATG"
STOPS = {"TAA", "TAG", "TGA"}


@dataclass
class ORF:
    start: int          # 0-based start on plus strand
    end: int            # exclusive end on plus strand
    strand: str         # '+' or '-'
    frame: int          # 0, 1, or 2
    length_nt: int
    length_aa: int
    protein: str        # translated amino acid sequence

    def to_dict(self) -> dict:
        return {
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "frame": self.frame,
            "length_nt": self.length_nt,
            "length_aa": self.length_aa,
        }


def _orfs_one_strand(seq_str: str, strand: str, min_aa: int, total_len: int) -> list[ORF]:
    """Find ORFs on a given strand by scanning all 3 reading frames."""
    orfs = []
    for frame in range(3):
        i = frame
        while i <= len(seq_str) - 3:
            codon = seq_str[i : i + 3]
            if codon == START:
                # Walk forward until a stop codon or end of sequence.
                j = i
                while j <= len(seq_str) - 3:
                    c = seq_str[j : j + 3]
                    if c in STOPS:
                        break
                    j += 3
                length_nt = j - i
                if length_nt // 3 >= min_aa:
                    protein_seq = str(Seq(seq_str[i:j]).translate())
                    if strand == "+":
                        start, end = i, j + 3
                    else:
                        # Map back to plus-strand coordinates.
                        start = total_len - (j + 3)
                        end = total_len - i
                    orfs.append(
                        ORF(
                            start=start,
                            end=end,
                            strand=strand,
                            frame=frame,
                            length_nt=length_nt,
                            length_aa=length_nt // 3,
                            protein=protein_seq,
                        )
                    )
                i = j + 3
            else:
                i += 3
    return orfs


def find_orfs(seq, min_aa: int = 100) -> list[ORF]:
    """Find all ORFs >= min_aa amino acids in 6 reading frames.

    Args:
        seq: A DNA string or Biopython Seq object.
        min_aa: Minimum protein length to keep (default 100 aa).
    """
    plus = str(seq).upper()
    minus = str(Seq(plus).reverse_complement())
    total = len(plus)
    orfs = _orfs_one_strand(plus, "+", min_aa, total)
    orfs += _orfs_one_strand(minus, "-", min_aa, total)
    orfs.sort(key=lambda o: o.start)
    return orfs
