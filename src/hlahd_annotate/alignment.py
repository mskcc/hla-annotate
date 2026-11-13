"""Exon 2+3 sequences from IPD-IMGT/HLA nucleotide alignments (alignments/<GENE>_nuc.txt).

Class I typing by HLA-HD scores exons 2 and 3 only, so the phase test compares
alleles over exactly those exons. The alignment files give every allele on the
same columns, so per-position comparison needs no aligner:

    ' A*01:01:01:01     ATG GCC ... G.|GC ...'   reference row (first allele)
    ' A*02:01:01:01     --- --- ... -.|-- ...'   '-' same as reference
                                                  '*' not sequenced, '.' gap
                                                  '|' exon boundary
Columns are returned as single characters: a base, '.' for a gap, '?' unknown.
"""
import re
from pathlib import Path

_NAME_RE = re.compile(r"^[A-Z0-9]+\*\d{2,3}(:\d{2,4}[A-Z]?){1,3}[NLSQCAG]?$")
TYPED_EXONS = (2, 3)  # 1-based exon numbers HLA-HD scores for class I


def _read_rows(path: Path) -> dict[str, str]:
    """Allele name -> concatenated aligned row (spaces removed, '|' kept)."""
    rows: dict[str, list[str]] = {}
    with open(path) as f:
        for line in f:
            if not line.startswith(' ') or not line.strip():
                continue
            parts = line.split()
            if not _NAME_RE.match(parts[0]):
                continue
            rows.setdefault(parts[0], []).append(''.join(parts[1:]))
    if not rows:
        raise ValueError(f"no allele rows found in IMGT alignment {path}")
    return {name: ''.join(chunks) for name, chunks in rows.items()}


def _resolve(row: str, ref: str) -> str:
    """Expand '-' to the reference base and '*' to '?'; keep '.' and '|'.

    A few expression variants (e.g. A*03:437Q, C*04:09L in IMGT 3.65.0) read
    past the reference stop codon, so their rows run longer than the reference
    row: the extra 3' columns have no reference base and are kept as written.
    A row shorter than the reference ends in unknown ('?') columns.
    """
    ref = ref.ljust(len(row), '.')
    row = row.ljust(len(ref), '*')
    out = []
    for c, r in zip(row, ref):
        if r == '|':  # exon boundaries always come from the reference row
            out.append('|')
            continue
        if c == '-' and r == '.':
            raise ValueError("alignment row has '-' where the reference has no base")
        out.append(r if c == '-' else '?' if c == '*' else c.upper())
    return ''.join(out)


def _typed_columns(resolved: str) -> str:
    """Columns of exons 2 and 3 (exon 1 is the segment before the first '|')."""
    segments = resolved.split('|')
    if len(segments) <= max(TYPED_EXONS) - 1:
        raise ValueError("alignment row has fewer exon segments than expected")
    return ''.join(segments[e - 1] for e in TYPED_EXONS)


def _match_name(allele: str, names: list[str]) -> str | None:
    """Exact IMGT name, else the first name extending it at a field boundary."""
    if allele in names:
        return allele
    for name in names:
        if name.startswith(allele + ':'):
            return name
    return None


def load_typed_exons(alignment_dir: Path, locus: str, alleles: set[str]) -> dict[str, str]:
    """Return {allele: exon 2+3 columns} for the requested alleles of one locus.

    Alleles absent from the alignment are left out; the caller counts them as
    not evaluable rather than treating them as identical.
    """
    path = Path(alignment_dir) / f"{locus}_nuc.txt"
    rows = _read_rows(path)
    names = list(rows)
    ref = rows[names[0]]
    out = {}
    for allele in alleles:
        name = _match_name(allele.removeprefix('HLA-'), names)
        if name is not None:
            out[allele] = _typed_columns(_resolve(rows[name], ref))
    return out


def phase_relation(x: tuple[str, str], y: tuple[str, str]) -> tuple[str, int]:
    """Compare two genotypes' exon 2+3 columns.

    Returns (relation, unknown_columns) where relation is:
      'identical'  each allele of x equals an allele of y over every known column
      'phase'      same unordered base pair at every known column, but the
                   chromosome assignment differs (a phase conjugate)
      'different'  some column carries different bases
    Columns where any of the four alleles is unknown ('?') are skipped and counted.
    """
    (x1, x2), (y1, y2) = x, y
    unknown = 0
    same_straight = same_cross = True
    for a, b, c, d in zip(x1, x2, y1, y2):
        if '?' in (a, b, c, d):
            unknown += 1
            continue
        if sorted((a, b)) != sorted((c, d)):
            return 'different', unknown
        same_straight &= (a == c and b == d)
        same_cross &= (a == d and b == c)
    if same_straight or same_cross:
        return 'identical', unknown
    return 'phase', unknown
