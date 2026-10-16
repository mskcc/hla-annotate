import sys
from pathlib import Path
import pandas as pd

from .pgroup import load_pgroup_table, lookup_pgroup
from .parse_final import parse_final_result, CLASS_I_LOCI, is_missing_allele
from .parse_est import parse_est, cross_check_allele

# Map allele_position to (best_pairs index, coverage key within that pair)
_POSITION_TO_COVERAGE = {
    'allele1': (0, 'allele1_coverage'),
    'allele2': (0, 'allele2_coverage'),
    'allele1_pair2': (1, 'allele1_coverage'),
    'allele2_pair2': (1, 'allele2_coverage'),
}

ANNOTATED_COLUMNS = [
    'sample', 'locus', 'allele', 'allele_position', 'resolution',
    'p_group', 'p_group_found', 'multiple_best_pairs', 'has_ambiguous_pair',
    'ambiguous_alleles', 'est_mismatch', 'est_mismatch_detail',
    'exon2_depth', 'exon2_incomp', 'exon2_incomplete',
    'exon3_depth', 'exon3_incomp', 'exon3_incomplete',
    'incomplete_bp', 'incomplete_coverage',
]


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


def annotate_sample(
    result_dir: Path,
    sample: str,
    pgroup_file: Path,
) -> pd.DataFrame:
    """
    Combine parse_final, parse_est, and P-group lookup into the annotated DataFrame.
    Emits a WARNING to stderr for each allele that fails est cross-check.

    Returns DataFrame with columns:
        sample, locus, allele, allele_position, resolution,
        p_group, p_group_found, multiple_best_pairs, has_ambiguous_pair,
        ambiguous_alleles, est_mismatch, est_mismatch_detail,
        exon2_depth, exon2_incomp, exon2_incomplete,
        exon3_depth, exon3_incomp, exon3_incomplete,
        incomplete_bp, incomplete_coverage
    """
    result_dir = Path(result_dir)
    final_file = result_dir / f"{sample}_final.result.txt"

    df = parse_final_result(final_file)
    pgroup_lookup = load_pgroup_table(pgroup_file)

    # Load est data for each class I locus (returns empty-flag dict if file missing)
    est_data: dict[str, dict] = {}
    for locus in CLASS_I_LOCI:
        est_file = result_dir / f"{sample}_{locus}.est.txt"
        est_data[locus] = parse_est(est_file)

    records = []
    for _, row in df.iterrows():
        locus = row['locus']
        allele = row['allele']
        est = est_data.get(locus, {})
        is_pair2 = 'pair2' in row['allele_position']

        # P-group lookup
        p_group, p_group_found = lookup_pgroup(allele, pgroup_lookup)

        # Cross-check (primary pair alleles only; pair2 alleles not cross-checked)
        est_mismatch = False
        est_mismatch_detail = None
        if not is_missing_allele(allele) and not is_pair2:
            match, detail = cross_check_allele(allele, est.get('best_pair_alleles', []))
            if not match:
                est_mismatch = True
                est_mismatch_detail = detail
                print(
                    f"WARNING: est_mismatch for sample={sample} locus={locus} "
                    f"allele_position={row['allele_position']}: {detail}",
                    file=sys.stderr,
                )

        # multiple_best_pairs: prefer est.txt (more precise) over final.result.txt count
        multiple_best_pairs = est.get('multiple_best_pairs', row['multiple_best_pairs'])

        # Per-allele exon coverage from est best_pairs
        coverage = {'exon2_depth': None, 'exon2_incomp': None, 'exon2_incomplete': False,
                    'exon3_depth': None, 'exon3_incomp': None, 'exon3_incomplete': False}
        best_pairs = est.get('best_pairs', [])
        # Unrecognized allele_position falls through to all-None coverage (intended silent-null)
        mapping = _POSITION_TO_COVERAGE.get(row['allele_position'])
        if mapping:
            idx, key = mapping
            if idx < len(best_pairs):
                coverage = dict(best_pairs[idx][key])
        incomplete_coverage = coverage['exon2_incomplete'] or coverage['exon3_incomplete']

        records.append({
            'sample': sample,
            'locus': locus,
            'allele': allele,
            'allele_position': row['allele_position'],
            'resolution': row['resolution'],
            'p_group': p_group,
            'p_group_found': p_group_found,
            'multiple_best_pairs': multiple_best_pairs,
            'has_ambiguous_pair': est.get('has_ambiguous_pair', False),
            'ambiguous_alleles': est.get('ambiguous_alleles'),
            'est_mismatch': est_mismatch,
            'est_mismatch_detail': est_mismatch_detail,
            'exon2_depth': coverage['exon2_depth'],
            'exon2_incomp': coverage['exon2_incomp'],
            'exon2_incomplete': coverage['exon2_incomplete'],
            'exon3_depth': coverage['exon3_depth'],
            'exon3_incomp': coverage['exon3_incomp'],
            'exon3_incomplete': coverage['exon3_incomplete'],
            'incomplete_bp': _incomplete_bp(coverage),
            'incomplete_coverage': incomplete_coverage,
        })

    df = pd.DataFrame(records, columns=ANNOTATED_COLUMNS)
    # nullable int so a missing est.txt doesn't turn bp counts into floats
    df['incomplete_bp'] = df['incomplete_bp'].astype('Int64')
    return df
