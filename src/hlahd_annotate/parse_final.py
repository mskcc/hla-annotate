from pathlib import Path
import re
import pandas as pd

CLASS_I_LOCI = {'A', 'B', 'C'}
_NULL_VALUES = {'Not typed', '-', ''}


def pair_position(index: int) -> str:
    """allele_position for the index-th allele of a row: allele1, allele2,
    allele1_pair2, allele2_pair2, allele1_pair3, ..."""
    pair, slot = divmod(index, 2)
    base = f"allele{slot + 1}"
    return base if pair == 0 else f"{base}_pair{pair + 1}"


def position_pair_index(position: str) -> tuple[int, int]:
    """Inverse of pair_position: (pair index from 0, slot 0/1)."""
    base, _, pair = position.partition('_pair')
    return (int(pair) - 1 if pair else 0), int(base.removeprefix('allele')) - 1
_ALLELE_RE = re.compile(
    r"(?:HLA-)?(?P<locus>[A-Za-z0-9]+)\*[0-9]{2,3}"
    r"(?::[0-9A-Za-z]+){1,3}"
)


def is_missing_allele(value) -> bool:
    """Return whether *value* is a supported missing allele scalar.

    This helper is deliberately scalar-only.  In particular, it does not call
    ``bool(pd.NA)`` and it does not turn an arbitrary non-string into text.
    """
    if value is None or value is pd.NA:
        return True
    if isinstance(value, str):
        return value.strip() in _NULL_VALUES
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        return False
    # numpy.bool_ is scalar too, but avoid importing numpy just for this check.
    return type(missing).__name__ in {'bool', 'bool_'} and bool(missing)


def normalize_allele(value, locus: str | None = None) -> str | None:
    """Normalize one final-result allele while preserving its spelling.

    Missing sentinels become ``None``.  Non-null values must look like an HLA
    allele and, when *locus* is supplied, belong to that locus.
    """
    if is_missing_allele(value):
        return None
    if not isinstance(value, str):
        raise ValueError(f"allele must be a string or missing scalar, got {type(value).__name__}")
    allele = value.strip()
    match = _ALLELE_RE.fullmatch(allele)
    if not match or (locus is not None and match.group('locus').upper() != locus.upper()):
        expected = f" for locus {locus}" if locus else ''
        raise ValueError(f"malformed HLA allele{expected}: {value!r}")
    return allele


def normalize_allele_pair(first, second, locus: str) -> tuple[str | None, str | None]:
    """Return a canonical allele pair for one final-result pair.

    HLA-HD uses a dash in the second slot to denote a homozygous call.  The
    dash is duplicated only when the first slot is a valid allele; a dash in
    the first slot remains a no-call.
    """
    first = normalize_allele(first, locus)
    second_is_dash = isinstance(second, str) and second.strip() == '-'
    second = normalize_allele(second, locus)
    if first is not None and second is None and second_is_dash:
        second = first
    return first, second


def parse_final_result(result_file: Path) -> pd.DataFrame:
    """
    Parse *_final.result.txt. Returns only class I loci (A, B, C).

    Returns DataFrame with columns:
        locus, allele, allele_position, resolution, multiple_best_pairs
    where:
        - allele: str with HLA- prefix, or None/NaN for a true no-call
        - allele_position: 'allele1', 'allele2', 'allele1_pair2', 'allele2_pair2'
        - resolution: number of colon-separated fields (0 if null)
        - multiple_best_pairs: True if row has >2 allele columns
    """
    rows = []
    with open(result_file) as f:
        for line in f:
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            parts = line.split('\t')
            locus = parts[0]
            if locus not in CLASS_I_LOCI:
                continue

            alleles = parts[1:]
            multiple_best_pairs = len(alleles) > 2

            # Every reported pair is kept (HLA-HD prints all frequency-tied
            # survivors). Keep an explicitly present trailing empty field, and
            # supply a missing second slot when a row ends mid-pair.
            normalized_pairs = [
                normalize_allele_pair(
                    alleles[i], alleles[i + 1] if i + 1 < len(alleles) else '', locus)
                for i in range(0, len(alleles), 2)
            ]

            normalized_alleles = [allele for pair in normalized_pairs for allele in pair]
            positions = [pair_position(i) for i in range(len(normalized_alleles))]
            for pos, allele_val in zip(positions, normalized_alleles):
                if allele_val is not None:
                    fields = allele_val.removeprefix('HLA-').split(':')
                    resolution = len(fields)
                else:
                    resolution = 0

                rows.append({
                    'locus': locus,
                    'allele': allele_val,
                    'allele_position': pos,
                    'resolution': resolution,
                    'multiple_best_pairs': multiple_best_pairs,
                })

    return pd.DataFrame(rows, columns=['locus', 'allele', 'allele_position',
                                       'resolution', 'multiple_best_pairs'])
