from pathlib import Path

from .parse_final import is_missing_allele, normalize_allele

_NULL_ALLELE = {'-', ''}

_EXONS = ('exon2', 'exon3')


def _parse_coverage(field: str) -> dict:
    """
    Parse a coverage field like 'exon2:262.426:comp.0,exon3:303.069:incomp.7'
    into {exonN_depth, exonN_incomp, exonN_incomplete} for exon2 and exon3.

    exonN_incomp is the count of incomplete bases (the N in comp.N/incomp.N),
    or None if absent/unparseable. exonN_incomplete is True for anything other
    than a clean 'comp'/'comp.0' status: incomp in any form (incomp.0, no
    count, bad count), comp with a nonzero count, or an unknown status.
    Missing/null field -> None values, not incomplete.
    """
    cov = {}
    for e in _EXONS:
        cov[f'{e}_depth'] = None
        cov[f'{e}_incomp'] = None
        cov[f'{e}_incomplete'] = False

    if is_missing_allele(field) or (isinstance(field, str) and field.strip() in _NULL_ALLELE):
        return cov
    if not isinstance(field, str):
        raise ValueError(f"coverage field must be a string or missing scalar, got {type(field).__name__}")

    for part in field.split(','):
        bits = part.split(':')
        if len(bits) != 3:
            continue
        exon, depth, status = bits[0].strip(), bits[1].strip(), bits[2].strip()
        if exon not in _EXONS:
            continue
        try:
            cov[f'{exon}_depth'] = float(depth)
        except ValueError:
            cov[f'{exon}_depth'] = None
        name, _, count = status.partition('.')
        try:
            n = int(count) if count else None
        except ValueError:
            n = None
        if name == 'comp':
            cov[f'{exon}_incomp'] = 0 if n is None else n
            cov[f'{exon}_incomplete'] = bool(n)
        else:
            cov[f'{exon}_incomp'] = n
            cov[f'{exon}_incomplete'] = True

    return cov


def _est_allele(token: str) -> str | None:
    """Normalize one est.txt allele token, or None for non-allele text.

    est.txt is advisory (flags and cross-check only), so HLA-HD placeholders
    such as "No candidate." mean "no allele here" rather than failing the
    whole sample.
    """
    try:
        return normalize_allele(token)
    except ValueError:
        return None


def _first_allele(field: str) -> str | None:
    """First allele in a comma-separated est column, HLA- stripped; None if null."""
    if is_missing_allele(field):
        return None
    if not isinstance(field, str):
        raise ValueError(f"allele field must be a string or missing scalar, got {type(field).__name__}")
    field = field.strip()
    if field in _NULL_ALLELE:
        return None
    first = field.split(',')[0].strip().removeprefix('HLA-')
    return _est_allele(first)


def parse_est(est_file: Path) -> dict:
    """
    Parse *_{A,B,C}.est.txt. Used for flag extraction and cross-validation only.
    Never used as the authoritative source for allele calls.

    Returns dict with:
        multiple_best_pairs (bool): True if #Best allele pair count > 1
        has_ambiguous_pair (bool): True if #Other ambiguous pair section present
        ambiguous_alleles (str | None): 'allele1 / allele2' from ambiguous section
        best_pair_alleles (list[str]): all alleles from best-pair data rows, HLA- stripped
        best_pairs (list[dict]): per best pair, with representative names
            allele1/allele2 (str | None) and allele1_coverage/allele2_coverage dicts
        best_pair_slots (list[tuple[list, list]]): per best pair, every allele in
            each slot's comma list (HLA- stripped); the pre-tiebreak candidates
        other_pairs (list[tuple[str, str]]): every '#Other ambiguous pair' row
            (informational: HLA-HD ranked these below the best pair)
        other_rows_unparsed (int): '#Other' rows that yielded no allele
        est_found (bool): False when the file does not exist
    """
    result = {
        'multiple_best_pairs': False,
        'has_ambiguous_pair': False,
        'ambiguous_alleles': None,
        'best_pair_alleles': [],
        'best_pairs': [],
        'best_pair_slots': [],
        'other_pairs': [],
        'other_rows_unparsed': 0,
        'est_found': est_file.exists(),
    }

    if not result['est_found']:
        return result

    with open(est_file) as f:
        lines = [line.rstrip('\n') for line in f]

    in_ambiguous_section = False
    for line in lines:
        if not line:
            continue
        if line.startswith('#Best allele pair'):
            # Format: "#Best allele pair\t<count>"
            parts = line.split('\t')
            count = int(parts[1]) if len(parts) > 1 else 1
            result['multiple_best_pairs'] = count > 1
            in_ambiguous_section = False
        elif line.startswith('#Other ambiguous pair'):
            result['has_ambiguous_pair'] = True
            in_ambiguous_section = True
        elif line.startswith('#'):
            in_ambiguous_section = False
        elif in_ambiguous_section:
            _add_other_pair(line.split('\t'), result)
        else:
            _add_best_pair(line.split('\t'), result)

    if result['other_pairs']:
        result['ambiguous_alleles'] = '; '.join(f"{a} / {b}" for a, b in result['other_pairs'])
    return result


def _slot_alleles(col: str) -> list[str]:
    """Every allele in one comma-listed est.txt slot, HLA- stripped; [] if null."""
    col = col.strip()
    if is_missing_allele(col):
        return []
    alleles = []
    for a in col.split(','):
        a = a.strip().removeprefix('HLA-')
        allele = _est_allele(a) if a else None
        if allele:
            alleles.append(allele)
    return alleles


def _add_other_pair(cols: list[str], result: dict) -> None:
    """'#Other ambiguous pair' data line (no coverage columns: HLA-HD ranked it
    below the best pair). A slot may be a comma list or '-'; the first allele
    represents it. Rows that yield no allele are counted, never silently lost."""
    if len(cols) < 2:
        result['other_rows_unparsed'] += 1
        return
    a1 = _first_allele(cols[0])
    a2 = _first_allele(cols[1]) or a1  # '-' second slot = homozygous
    if a1:
        result['other_pairs'].append((a1, a2))
    else:
        result['other_rows_unparsed'] += 1


def _add_best_pair(cols: list[str], result: dict) -> None:
    """Best pair data line: col0/col1 comma lists, col2/col3 coverage.

    A homozygous pair has '-' for slot 2 and for its coverage (HLA-HD
    Select.h writes one coverage column), so allele 2 reuses allele 1's.
    """
    slot1 = _slot_alleles(cols[0]) if len(cols) > 0 else []
    slot2 = _slot_alleles(cols[1]) if len(cols) > 1 else []
    cov1 = _parse_coverage(cols[2]) if len(cols) > 2 else _parse_coverage('')
    cov2 = _parse_coverage(cols[3]) if len(cols) > 3 else _parse_coverage('')
    if slot1 and not slot2 and cov2['exon2_depth'] is None and cov2['exon3_depth'] is None:
        cov2 = dict(cov1)
    result['best_pair_alleles'].extend(slot1 + slot2)
    result['best_pair_slots'].append((slot1, slot2))
    result['best_pairs'].append({
        'allele1': _first_allele(cols[0]) if len(cols) > 0 else None,
        'allele2': _first_allele(cols[1]) if len(cols) > 1 else None,
        'allele1_coverage': cov1,
        'allele2_coverage': cov2,
    })


def cross_check_allele(
    allele: str | None,
    best_pair_alleles: list[str],
) -> tuple[bool, str | None]:
    """
    Check that an allele from final.result.txt appears in the est.txt best pair allele lists.
    Matching is performed at the field depth of the reported allele (2 or 3 fields).

    Returns (match_found, detail_string).
    If match_found is True, detail_string is None.
    Hard warning to stderr is the caller's responsibility.
    """
    allele = normalize_allele(allele)
    if is_missing_allele(allele) or not best_pair_alleles:
        return True, None

    a = allele.removeprefix('HLA-')
    fields = a.split(':')
    n_fields = len(fields)

    truncated = set()
    for bp in best_pair_alleles:
        bp = normalize_allele(bp)
        if is_missing_allele(bp):
            continue
        bp_fields = bp.split(':')
        truncated.add(':'.join(bp_fields[:n_fields]))

    if a in truncated:
        return True, None

    shown = sorted(truncated)[:10]
    detail = f"final={allele}; est_best_alleles={','.join(shown)}"
    return False, detail
