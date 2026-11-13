import sys
from pathlib import Path
import pandas as pd

from .alignment import load_typed_exons
from .dictionary import fmt_exons, load_exon_presence, missing_exons
from .pgroup import load_pgroup_table, lookup_pgroup
from .parse_final import parse_final_result, CLASS_I_LOCI, is_missing_allele, position_pair_index
from .parse_est import parse_est, cross_check_allele
from .reporting import classify, two_field

DEFAULT_MIN_EXON_DEPTH = 75.0
ALL_HOMOZYGOUS_FLAG = 'All loci homozygous'
INCOMPLETE_ALLELE_FLAG = 'Incomplete Allele'

_EMPTY_COVERAGE = {'exon2_depth': None, 'exon2_incomp': None, 'exon2_incomplete': False,
                   'exon3_depth': None, 'exon3_incomp': None, 'exon3_incomplete': False}

ANNOTATED_COLUMNS = [
    'sample', 'locus', 'allele', 'allele_position', 'resolution',
    'p_group', 'p_group_found', 'multiple_best_pairs', 'has_ambiguous_pair',
    'ambiguous_alleles', 'est_mismatch', 'est_mismatch_detail',
    'exon2_depth', 'exon2_incomp', 'exon2_incomplete',
    'exon3_depth', 'exon3_incomp', 'exon3_incomplete',
    'incomplete_bp', 'incomplete_coverage',
    # Reporting model (locus-level values, repeated on each allele row).
    'two_field', 'two_field_genotype', 'call_status', 'flags',
    'frequency_resolution_details', 'provenance', 'phase_check', 'qc_evaluable',
    'incomplete_allele_detail',
]
# multiple_best_pairs / has_ambiguous_pair are raw est.txt indicators ('#Best allele
# pair' > 1, '#Other ambiguous pair' present). Report call_status and flags instead.


def _incomplete_bp(coverage: dict) -> int | None:
    """Total incomplete bases across exons 2 and 3.

    None when there is no coverage, or when an exon is flagged incomplete
    but its count is unknown (the total would understate it).
    """
    if coverage['exon2_depth'] is None and coverage['exon3_depth'] is None:
        return None
    total = 0
    for exon in ('exon2', 'exon3'):
        n = coverage[f'{exon}_incomp']
        if n is None:
            if coverage[f'{exon}_incomplete']:
                return None
            continue
        total += n
    return total


def _in_slot(allele: str, slot: list[str]) -> bool:
    """Field-boundary match of a reported allele against an est.txt slot list."""
    a = allele.removeprefix('HLA-')
    return any(x == a or x.startswith(a + ':') for x in slot)


def _coverage_for(allele, mate, est: dict, slot: int = 0) -> dict:
    """Coverage of `allele` from the est.txt best-pair line holding it and its mate.

    pick_up_allele sorts candidates by frequency, so final pair k need not come
    from est line k; match by allele names instead of by position.
    """
    if is_missing_allele(allele):
        return dict(_EMPTY_COVERAGE)
    lines = list(zip(est.get('best_pair_slots', []), est.get('best_pairs', [])))
    for want_mate in (mate, None):  # the pair's own line first, then any line
        for (s1, s2), bp in lines:
            s2 = s2 or s1  # homozygous line
            if _in_slot(allele, s1) and (want_mate is None or _in_slot(want_mate, s2)):
                # homozygous: both rows name the same allele; keep each row's own column
                key = f'allele{slot + 1}_coverage' if _in_slot(allele, [mate or '']) else 'allele1_coverage'
                return dict(bp[key])
            if _in_slot(allele, s2) and (want_mate is None or _in_slot(want_mate, s1)):
                return dict(bp['allele2_coverage'])
    return dict(_EMPTY_COVERAGE)


def _allele_record(sample, locus, row, est, pgroup_lookup, coverage) -> dict:
    allele = row['allele']
    p_group, p_group_found = lookup_pgroup(allele, pgroup_lookup)
    est_mismatch, est_mismatch_detail = False, None
    # Cross-check primary pair alleles only.
    if not is_missing_allele(allele) and 'pair' not in row['allele_position']:
        match, detail = cross_check_allele(allele, est.get('best_pair_alleles', []))
        if not match:
            est_mismatch, est_mismatch_detail = True, detail
            print(f"WARNING: est_mismatch for sample={sample} locus={locus} "
                  f"allele_position={row['allele_position']}: {detail}", file=sys.stderr)
    return {
        'sample': sample, 'locus': locus, 'allele': allele,
        'allele_position': row['allele_position'], 'resolution': row['resolution'],
        'p_group': p_group, 'p_group_found': p_group_found,
        # prefer est.txt (more precise) over the final.result.txt column count
        'multiple_best_pairs': est.get('multiple_best_pairs', row['multiple_best_pairs']),
        'has_ambiguous_pair': est.get('has_ambiguous_pair', False),
        'ambiguous_alleles': est.get('ambiguous_alleles'),
        'est_mismatch': est_mismatch, 'est_mismatch_detail': est_mismatch_detail,
        **coverage,
        'incomplete_bp': _incomplete_bp(coverage),
        'incomplete_coverage': coverage['exon2_incomplete'] or coverage['exon3_incomplete'],
        'two_field': None if is_missing_allele(allele) else two_field(allele),
    }


def _locus_pairs(rows: list[dict]) -> list[tuple]:
    """Reported (a1, a2) pairs of one locus, HLA- stripped, None for no-calls."""
    pairs: dict[int, list] = {}
    for r in rows:
        idx, slot = position_pair_index(r['allele_position'])
        a = None if is_missing_allele(r['allele']) else r['allele'].removeprefix('HLA-')
        pairs.setdefault(idx, [None, None])[slot] = a
    return [tuple(pairs[i]) for i in sorted(pairs)]


def _typed_exons(alignment_dir, locus, est, pairs):
    if alignment_dir is None:
        return None
    alleles = {a for p in pairs for a in p if a}
    alleles |= set(est.get('best_pair_alleles', []))
    alleles |= {a for p in est.get('other_pairs', []) for a in p}  # phase test needs these
    try:
        return load_typed_exons(alignment_dir, locus, alleles)
    except FileNotFoundError:
        print(f"WARNING: no IMGT alignment {locus}_nuc.txt in {alignment_dir}; "
              f"phase not evaluated for locus {locus}", file=sys.stderr)
        return None


def _exon_table(dictionary_dir, locus):
    if dictionary_dir is None:
        return None
    try:
        return load_exon_presence(dictionary_dir, locus)
    except FileNotFoundError:
        print(f"WARNING: no HLA-HD dictionary list {locus}/{locus}.list in {dictionary_dir}; "
              f"Incomplete Allele not evaluated for locus {locus}", file=sys.stderr)
        return None


def _incomplete_allele(pairs, table) -> tuple[bool, str]:
    """(flag, detail) for the reported alleles: an allele lacking exons in the
    HLA-HD dictionary is flagged; one absent from it is named, never passed."""
    alleles = sorted({a for p in pairs for a in p if a})
    if not alleles:
        return False, ''
    if table is None:
        return False, 'not evaluated'
    notes, flag = [], False
    for a in alleles:
        missing = missing_exons(a, table)
        if missing is None:
            notes.append(f"{a} not in HLA-HD dictionary")
        elif missing:
            flag = True
            notes.append(f"{a} missing exon{'s' if len(missing) > 1 else ''} {fmt_exons(missing)}")
    return flag, '; '.join(notes)


_FAILED_ROWS = [{'allele': None, 'allele_position': p, 'resolution': 0,
                 'multiple_best_pairs': False} for p in ('allele1', 'allele2')]


def _annotate_locus(sample, locus, rows, result_dir, pgroup_lookup, alignment_dir, min_depth,
                    exon_table=None):
    """(allele records, reporting model dict) for one locus."""
    est = parse_est(result_dir / f"{sample}_{locus}.est.txt")
    pairs = _locus_pairs(rows)
    recs = []
    for r in rows:
        idx, slot = position_pair_index(r['allele_position'])
        mate = pairs[idx][1 - slot]
        cov = _coverage_for(r['allele'], mate, est, slot)
        recs.append(_allele_record(sample, locus, r, est, pgroup_lookup, cov))
    reported_cov = [r for r in recs if not is_missing_allele(r['allele'])]
    model = classify(pairs, est, pgroup_lookup, _typed_exons(alignment_dir, locus, est, pairs),
                     reported_cov, min_depth)
    flag, model['incomplete_allele_detail'] = _incomplete_allele(pairs, exon_table)
    if flag:
        model['flags'].append(INCOMPLETE_ALLELE_FLAG)
    return recs, model


def _flag_all_homozygous(models: dict) -> None:
    """All class I loci typed and homozygous at two fields -> QC review."""
    if set(models) != set(CLASS_I_LOCI):
        return
    for m in models.values():
        g = m['two_field_genotype']
        if 'Failed' in m['call_status'] or not g or ';' in g or len(set(g.split('/'))) != 1:
            return
    for m in models.values():
        m['flags'].append(ALL_HOMOZYGOUS_FLAG)
        if 'QC review' not in m['call_status']:
            m['call_status'].append('QC review')


def annotate_sample(
    result_dir: Path,
    sample: str,
    pgroup_file: Path,
    alignment_dir: Path | None = None,
    min_exon_depth: float = DEFAULT_MIN_EXON_DEPTH,
    dictionary_dir: Path | None = None,
) -> pd.DataFrame:
    """
    Combine parse_final, parse_est, P-group lookup and the reporting model into
    one row per reported allele (ANNOTATED_COLUMNS). alignment_dir holds IMGT
    alignments/<GENE>_nuc.txt for the phase test; without it phase_check is
    'not_evaluated'. dictionary_dir is HLA-HD's dictionary/ (holds <GENE>/<GENE>.list)
    for the Incomplete Allele flag; without it incomplete_allele_detail is
    'not evaluated'. A class I locus missing from final.result.txt is reported
    as Failed (with a warning), never dropped.
    """
    result_dir = Path(result_dir)
    df = parse_final_result(result_dir / f"{sample}_final.result.txt")
    pgroup_lookup = load_pgroup_table(pgroup_file)

    per_locus = {}
    for locus in sorted(CLASS_I_LOCI):
        rows = df[df['locus'] == locus].to_dict('records')
        if not rows:
            print(f"WARNING: locus {locus} missing from final.result.txt for sample={sample}; "
                  f"reported as Failed", file=sys.stderr)
            rows = [dict(r) for r in _FAILED_ROWS]
        per_locus[locus] = _annotate_locus(sample, locus, rows, result_dir, pgroup_lookup,
                                           alignment_dir, min_exon_depth,
                                           _exon_table(dictionary_dir, locus))

    _flag_all_homozygous({locus: model for locus, (_, model) in per_locus.items()})
    records = []
    for recs, model in per_locus.values():
        model = {**model, 'call_status': ';'.join(model['call_status']),
                 'flags': ';'.join(model['flags'])}
        records += [{**r, **model} for r in recs]

    df = pd.DataFrame(records, columns=ANNOTATED_COLUMNS)
    # nullable int so a missing est.txt doesn't turn bp counts into floats
    df['incomplete_bp'] = df['incomplete_bp'].astype('Int64')
    return df
