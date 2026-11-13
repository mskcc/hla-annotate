import sys
from pathlib import Path
import pandas as pd

from .alignment import load_typed_exons
from .pgroup import load_pgroup_table, lookup_pgroup
from .parse_final import parse_final_result, CLASS_I_LOCI, is_missing_allele, position_pair_index
from .parse_est import parse_est, cross_check_allele
from .reporting import classify, two_field

DEFAULT_MIN_EXON_DEPTH = 75.0

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


def _coverage(est: dict, position: str) -> dict:
    pair_idx, slot = position_pair_index(position)
    best_pairs = est.get('best_pairs', [])
    if pair_idx < len(best_pairs):
        return dict(best_pairs[pair_idx][f'allele{slot + 1}_coverage'])
    return dict(_EMPTY_COVERAGE)


def _allele_record(sample, locus, row, est, pgroup_lookup) -> dict:
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
    coverage = _coverage(est, row['allele_position'])
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


def _locus_pairs(records: list[dict]) -> list[tuple]:
    """Reported (a1, a2) pairs of one locus, HLA- stripped, None for no-calls."""
    pairs: dict[int, list] = {}
    for r in records:
        idx, slot = position_pair_index(r['allele_position'])
        a = None if is_missing_allele(r['allele']) else r['allele'].removeprefix('HLA-')
        pairs.setdefault(idx, [None, None])[slot] = a
    return [tuple(pairs[i]) for i in sorted(pairs)]


def _typed_exons(alignment_dir, locus, est, pairs):
    if alignment_dir is None:
        return None
    alleles = {a for p in pairs for a in p if a}
    alleles |= set(est.get('best_pair_alleles', []))
    alleles |= {a for p in est.get('other_pairs', []) for a in p}
    return load_typed_exons(alignment_dir, locus, alleles)


def annotate_sample(
    result_dir: Path,
    sample: str,
    pgroup_file: Path,
    alignment_dir: Path | None = None,
    min_exon_depth: float = DEFAULT_MIN_EXON_DEPTH,
) -> pd.DataFrame:
    """
    Combine parse_final, parse_est, P-group lookup and the reporting model into
    one row per reported allele (ANNOTATED_COLUMNS). alignment_dir holds IMGT
    alignments/<GENE>_nuc.txt for the phase test; without it phase_check is
    'not_evaluated'. Emits a WARNING to stderr for each est cross-check failure.
    """
    result_dir = Path(result_dir)
    df = parse_final_result(result_dir / f"{sample}_final.result.txt")
    pgroup_lookup = load_pgroup_table(pgroup_file)

    records = []
    for locus in sorted(CLASS_I_LOCI):
        locus_rows = df[df['locus'] == locus]
        if locus_rows.empty:
            continue
        est = parse_est(result_dir / f"{sample}_{locus}.est.txt")
        recs = [_allele_record(sample, locus, row, est, pgroup_lookup)
                for _, row in locus_rows.iterrows()]
        pairs = _locus_pairs(recs)
        primary_cov = [r for r in recs if r['allele_position'] in ('allele1', 'allele2')
                       and not is_missing_allele(r['allele'])]
        model = classify(pairs, est, pgroup_lookup, _typed_exons(alignment_dir, locus, est, pairs),
                         primary_cov, min_exon_depth)
        model['call_status'] = ';'.join(model['call_status'])
        model['flags'] = ';'.join(model['flags'])
        records += [{**r, **model} for r in recs]

    df = pd.DataFrame(records, columns=ANNOTATED_COLUMNS)
    # nullable int so a missing est.txt doesn't turn bp counts into floats
    df['incomplete_bp'] = df['incomplete_bp'].astype('Int64')
    return df
