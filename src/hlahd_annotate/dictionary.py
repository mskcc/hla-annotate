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


class ExonTable(dict):
    """{allele name without 'HLA-': [missing exon numbers]}; .exons holds the locus's exon numbers."""
    exons: tuple[int, ...] = ()


def load_exon_presence(dict_dir: Path, locus: str) -> ExonTable:
    """{allele name without 'HLA-': [missing exon numbers]} for one locus (an ExonTable)."""
    path = Path(dict_dir) / locus / f"{locus}.list"
    rows = [line.rstrip('\n').split('\t') for line in open(path)]
    header = next((r for r in rows if r and r[0].startswith('#HLA name')), None)
    if header is None:
        raise ValueError(f"no '#HLA name' header in HLA-HD dictionary list {path}")
    exons = [(i, int(h.removeprefix('exon'))) for i, h in enumerate(header) if h.startswith('exon')]
    if not exons:
        raise ValueError(f"no exon columns in HLA-HD dictionary list {path}")
    table = ExonTable()
    table.exons = tuple(n for _, n in exons)
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


def known_beyond_exons_2_3(allele: str, table: ExonTable) -> bool | None:
    """True if the dictionary has sequence for this allele outside exons 2-3, False if
    exons 2-3 are all it has, None if the allele is not in the dictionary.

    HLA-HD scores an allele on every exon it has (Tree.h: a branch extends to the exons
    of its alleles), so an '#Other' alternative with such an allele was tested there and lost.
    """
    missing = missing_exons(allele, table)
    if missing is None:
        return None
    return bool(set(table.exons) - {2, 3} - set(missing))


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
