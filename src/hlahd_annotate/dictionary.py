"""Exon presence per allele from HLA-HD's dictionary (<dict>/<GENE>/<GENE>.list).

Each .list row is one allele with '+'/'-' for the genomic sequence, the CDS and
each exon, as HLA-HD built it from IPD-IMGT/HLA:

    #HLA name  Allele name        gene  nuc  exon1 ... exon8
    HLA00005   HLA-A*02:01:01:01  +     +    +     ... +

A '-' exon is absent from the official sequence (many alleles are known only
for exons 2-3). An allele whose every dictionary entry lacks an exon is an
Incomplete Allele.
"""
from pathlib import Path


def load_exon_presence(dict_dir: Path, locus: str) -> dict[str, list[int]]:
    """{allele name without 'HLA-': [missing exon numbers]} for one locus."""
    path = Path(dict_dir) / locus / f"{locus}.list"
    rows = [line.rstrip('\n').split('\t') for line in open(path)]
    header = next((r for r in rows if r and r[0].startswith('#HLA name')), None)
    if header is None:
        raise ValueError(f"no '#HLA name' header in HLA-HD dictionary list {path}")
    exons = [(i, int(h.removeprefix('exon'))) for i, h in enumerate(header) if h.startswith('exon')]
    if not exons:
        raise ValueError(f"no exon columns in HLA-HD dictionary list {path}")
    table = {}
    for r in rows:
        if not r or r[0].startswith('#') or len(r) < len(header):
            continue
        table[r[1].removeprefix('HLA-')] = [n for i, n in exons if r[i] != '+']
    if not table:
        raise ValueError(f"no allele rows in HLA-HD dictionary list {path}")
    return table


def missing_exons(allele: str, table: dict[str, list[int]]) -> list[int] | None:
    """Exons missing from a reported allele; [] if complete, None if not in the dictionary.

    A reported name may cover several dictionary entries (A*02:01:01 covers
    A*02:01:01:01, ...); it counts as incomplete only when every one of them
    lacks an exon, and the most complete entry's missing exons are returned.
    """
    a = allele.removeprefix('HLA-')
    hits = [v for k, v in table.items() if k == a or k.startswith(a + ':')]
    if not hits:
        return None
    return min(hits, key=len)


def fmt_exons(exons: list[int]) -> str:
    """[1, 5, 6, 7] -> '1, 5-7'."""
    out, start = [], None
    for i, n in enumerate(exons):
        if start is None:
            start = n
        if i + 1 == len(exons) or exons[i + 1] != n + 1:
            out.append(str(start) if start == n else f"{start}-{n}")
            start = None
    return ', '.join(out)
