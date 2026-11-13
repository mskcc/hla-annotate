import shutil
import pandas as pd
from hlahd_annotate.annotate import annotate_sample

EXPECTED_COLUMNS = [
    'sample', 'locus', 'allele', 'allele_position', 'resolution',
    'p_group', 'p_group_found', 'multiple_best_pairs', 'has_ambiguous_pair',
    'ambiguous_alleles', 'est_mismatch', 'est_mismatch_detail',
    'exon2_depth', 'exon2_incomp', 'exon2_incomplete',
    'exon3_depth', 'exon3_incomp', 'exon3_incomplete',
    'incomplete_bp', 'incomplete_coverage',
    'two_field', 'two_field_genotype', 'call_status', 'flags',
    'frequency_resolution_details', 'provenance', 'phase_check', 'qc_evaluable',
]


def test_annotate_returns_dataframe(fixtures_dir, pgroup_file):
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    assert isinstance(df, pd.DataFrame)


def test_annotate_has_expected_columns(fixtures_dir, pgroup_file):
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    for col in EXPECTED_COLUMNS:
        assert col in df.columns, f"Missing column: {col}"


def test_annotate_sample_column_populated(fixtures_dir, pgroup_file):
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    assert (df['sample'] == 'sample3').all()


def test_annotate_class_i_only(fixtures_dir, pgroup_file):
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    assert set(df['locus'].unique()).issubset({'A', 'B', 'C'})


def test_annotate_pgroup_found(fixtures_dir, pgroup_file):
    # sample3 A locus allele1 is HLA-A*02:01:01 -> A*02:01P in mini pgroup fixture
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    a1 = df[(df['locus'] == 'A') & (df['allele_position'] == 'allele1')].iloc[0]
    assert a1['p_group_found'] == True  # noqa: E712 — numpy bool comparison
    assert a1['p_group'] == 'A*02:01P'


def test_annotate_no_flags_for_normal_sample(fixtures_dir, pgroup_file):
    # sample3: B has multiple columns in final but est.txt is missing -> defaults False
    # A has est.txt with 1 best pair -> multiple_best_pairs=False, has_ambiguous_pair=False
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    assert not df['multiple_best_pairs'].any()
    assert not df['has_ambiguous_pair'].any()
    assert not df['est_mismatch'].any()


def test_annotate_est_mismatch_raises_warning(fixtures_dir, pgroup_file, tmp_path, capsys):
    # Create a final.result.txt with alleles NOT in the est.txt best pair alleles
    f = tmp_path / "bad_final.result.txt"
    f.write_text("A\tHLA-A*99:99:99\tHLA-A*99:99:88\n")
    # Copy sample1_A.est.txt as the est.txt for locus A
    shutil.copy(fixtures_dir / "sample1_A.est.txt", tmp_path / "bad_A.est.txt")
    df = annotate_sample(tmp_path, 'bad', pgroup_file)
    captured = capsys.readouterr()
    assert 'WARNING' in captured.err
    assert df['est_mismatch'].any()
    mismatch_rows = df[df['est_mismatch']]
    assert mismatch_rows['est_mismatch_detail'].notna().all()


def test_annotate_tsv_roundtrip(fixtures_dir, pgroup_file, tmp_path):
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    tsv_path = tmp_path / "sample3_annotated.tsv"
    df.to_csv(tsv_path, sep='\t', index=False)
    df2 = pd.read_csv(tsv_path, sep='\t')
    assert list(df2.columns) == list(df.columns)
    assert len(df2) == len(df)


COVERAGE_COLUMNS = ['exon2_depth', 'exon2_incomp', 'exon2_incomplete',
                    'exon3_depth', 'exon3_incomp', 'exon3_incomplete',
                    'incomplete_bp', 'incomplete_coverage']


def test_annotate_has_coverage_columns(fixtures_dir, pgroup_file):
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    for col in COVERAGE_COLUMNS:
        assert col in df.columns, f"Missing column: {col}"


def test_annotate_coverage_values_complete(fixtures_dir, pgroup_file):
    # sample3 A best pair has comp.0 on both exons -> not incomplete
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    a1 = df[(df['locus'] == 'A') & (df['allele_position'] == 'allele1')].iloc[0]
    assert a1['exon2_depth'] == 262.426
    assert a1['exon2_incomp'] == 0
    assert bool(a1['incomplete_coverage']) is False


def test_annotate_incomplete_coverage_flag(fixtures_dir, pgroup_file):
    # sample5 A: allele1 exon3 incomp.7, allele2 exon2 incomp.4
    df = annotate_sample(fixtures_dir, 'sample5', pgroup_file)
    a1 = df[(df['locus'] == 'A') & (df['allele_position'] == 'allele1')].iloc[0]
    a2 = df[(df['locus'] == 'A') & (df['allele_position'] == 'allele2')].iloc[0]
    assert a1['exon3_incomp'] == 7
    assert bool(a1['incomplete_coverage']) is True
    assert a2['exon2_incomp'] == 4
    assert bool(a2['incomplete_coverage']) is True


def test_annotate_coverage_none_when_est_missing(fixtures_dir, pgroup_file):
    # sample3 B locus has no est.txt -> coverage columns are None, flag False
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    b_rows = df[df['locus'] == 'B']
    assert b_rows['exon2_depth'].isna().all()
    assert not b_rows['incomplete_coverage'].any()


# --- homozygous + "No candidate." (ported from hlaHD_workflow 1.0.1 / 1.0.2) ---

def test_annotate_column_order_is_fixed(fixtures_dir, pgroup_file):
    from hlahd_annotate.annotate import ANNOTATED_COLUMNS
    assert ANNOTATED_COLUMNS == EXPECTED_COLUMNS
    df = annotate_sample(fixtures_dir, 'homozyg', pgroup_file)
    assert list(df.columns) == EXPECTED_COLUMNS


def test_annotate_missing_class_i_loci_are_failed_not_dropped(pgroup_file, tmp_path, capsys):
    # e.g. a bare "Couldn't read result file." line replaced the locus line
    (tmp_path / 'empty_final.result.txt').write_text("DRB1\tHLA-DRB1*15:01:01\t-\n")
    df = annotate_sample(tmp_path, 'empty', pgroup_file)
    assert list(df.columns) == EXPECTED_COLUMNS
    assert sorted(df['locus'].unique()) == ['A', 'B', 'C']
    assert set(df['call_status']) == {'Failed'} and df['allele'].isna().all()
    assert capsys.readouterr().err.count('reported as Failed') == 3


def test_annotate_homozygous_keeps_both_alleles(fixtures_dir, pgroup_file):
    df = annotate_sample(fixtures_dir, 'homozyg', pgroup_file)
    a = df[df['locus'] == 'A'].set_index('allele_position')
    for pos in ('allele1', 'allele2'):
        assert a.loc[pos, 'allele'] == 'HLA-A*02:01:01'
        assert a.loc[pos, 'p_group'] == 'A*02:01P'
        assert bool(a.loc[pos, 'p_group_found']) is True
        assert bool(a.loc[pos, 'est_mismatch']) is False
    assert a.loc['allele2', 'exon2_depth'] == 257.163


def test_annotate_no_candidate_est_does_not_fail_sample(fixtures_dir, pgroup_file):
    df = annotate_sample(fixtures_dir, 'homozyg', pgroup_file)
    assert len(df) == 6
    c1 = df[(df['locus'] == 'C') & (df['allele_position'] == 'allele1')].iloc[0]
    assert bool(c1['est_mismatch']) is False
    assert c1['exon2_depth'] == 262.426


# --- incomplete coverage: flag any incompleteness, report incomplete bp ---

def test_annotate_incomplete_bp_sums_exons(fixtures_dir, pgroup_file):
    # sample5 A: allele1 exon3 incomp.7 (exon2 comp.0), allele2 exon2 incomp.4
    df = annotate_sample(fixtures_dir, 'sample5', pgroup_file)
    a = df[df['locus'] == 'A'].set_index('allele_position')
    assert a.loc['allele1', 'incomplete_bp'] == 7
    assert bool(a.loc['allele1', 'exon3_incomplete']) is True
    assert bool(a.loc['allele1', 'exon2_incomplete']) is False
    assert a.loc['allele2', 'incomplete_bp'] == 4
    assert bool(a.loc['allele2', 'exon2_incomplete']) is True


def test_annotate_complete_has_zero_incomplete_bp(fixtures_dir, pgroup_file):
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    a1 = df[(df['locus'] == 'A') & (df['allele_position'] == 'allele1')].iloc[0]
    assert a1['incomplete_bp'] == 0
    assert bool(a1['incomplete_coverage']) is False


def test_annotate_incomplete_bp_missing_when_est_missing(fixtures_dir, pgroup_file):
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    b_rows = df[df['locus'] == 'B']
    assert b_rows['incomplete_bp'].isna().all()


def _write_sample(tmp_path, cov1, cov2):
    (tmp_path / 'inc_final.result.txt').write_text(
        "A\tHLA-A*02:01:01\tHLA-A*03:01:01\n")
    (tmp_path / 'inc_A.est.txt').write_text(
        "#Best allele pair\n"
        f"HLA-A*02:01:01\tHLA-A*03:01:01\t{cov1}\t{cov2}\n")


def test_annotate_flags_incomp_zero(pgroup_file, tmp_path):
    # "anything at all": an incomp status is flagged even when its count is 0
    _write_sample(tmp_path, 'exon2:200.0:incomp.0,exon3:210.0:comp.0',
                  'exon2:200.0:comp.0,exon3:210.0:comp.0')
    df = annotate_sample(tmp_path, 'inc', pgroup_file)
    a = df[df['locus'] == 'A'].set_index('allele_position')
    assert bool(a.loc['allele1', 'incomplete_coverage']) is True
    assert a.loc['allele1', 'incomplete_bp'] == 0
    assert bool(a.loc['allele2', 'incomplete_coverage']) is False


def test_annotate_flags_incomp_without_count(pgroup_file, tmp_path):
    # unparseable count: still flagged, bp unknown for that exon
    _write_sample(tmp_path, 'exon2:200.0:incomp,exon3:210.0:incomp.x',
                  'exon2:200.0:comp.0,exon3:210.0:incomp.3')
    df = annotate_sample(tmp_path, 'inc', pgroup_file)
    a = df[df['locus'] == 'A'].set_index('allele_position')
    assert bool(a.loc['allele1', 'exon2_incomplete']) is True
    assert bool(a.loc['allele1', 'exon3_incomplete']) is True
    assert bool(a.loc['allele1', 'incomplete_coverage']) is True
    assert pd.isna(a.loc['allele1', 'incomplete_bp'])
    assert a.loc['allele2', 'incomplete_bp'] == 3


# --- review findings 2026-11-13: real HLA-HD shapes ---

COV_GOOD = 'exon2:212.4:comp.0,exon3:198.7:comp.0'
COV_LOW = 'exon2:30.0:comp.0,exon3:31.0:comp.0'


def _write(tmp_path, sample, final, ests):
    (tmp_path / f'{sample}_final.result.txt').write_text(final)
    for locus, text in ests.items():
        (tmp_path / f'{sample}_{locus}.est.txt').write_text(text)


def _est(rows, other=()):
    out = [f"#Pair count\t{len(rows)}", f"#Best allele pair\t{len(rows)}", *rows]
    if other:
        out += [f"#Other ambiguous pair\t{len(other)}", *other]
    return '\n'.join(out) + '\n'


HET_B = _est([f"HLA-B*07:02:01:01\tHLA-B*44:02:01:01\t{COV_GOOD}\t{COV_GOOD}"])
HET_C = _est([f"HLA-C*07:02:01:01\tHLA-C*05:01:01:01\t{COV_GOOD}\t{COV_GOOD}"])
HET_FINAL_BC = "B\tHLA-B*07:02:01\tHLA-B*44:02:01\nC\tHLA-C*07:02:01\tHLA-C*05:01:01\n"


def test_real_homozygous_est_is_not_sent_to_review(pgroup_file, tmp_path):
    # HLA-HD writes '-' for allele 2 AND for its coverage column.
    _write(tmp_path, 'S', "A\tHLA-A*02:01:01\t-\n" + HET_FINAL_BC,
           {'A': _est([f"HLA-A*02:01:01:01\t-\t{COV_GOOD}\t-"]), 'B': HET_B, 'C': HET_C})
    a = annotate_sample(tmp_path, 'S', pgroup_file).set_index(['locus', 'allele_position'])
    assert a.loc[('A', 'allele2'), 'exon2_depth'] == 212.4
    assert a.loc[('A', 'allele1'), 'call_status'] == 'Unique'
    assert a.loc[('A', 'allele1'), 'flags'] == ''


def test_coverage_taken_from_the_est_line_of_the_reported_pair(pgroup_file, tmp_path):
    # Frequency picked line 2's pair; line 1 has low depth and must not leak in.
    rows = [f"HLA-A*68:01:01:01\tHLA-A*03:01:01:01\t{COV_LOW}\t{COV_LOW}",
            f"HLA-A*02:01:01:01\tHLA-A*03:01:01:01\t{COV_GOOD}\t{COV_GOOD}"]
    _write(tmp_path, 'S', "A\tHLA-A*02:01:01\tHLA-A*03:01:01\n" + HET_FINAL_BC,
           {'A': _est(rows), 'B': HET_B, 'C': HET_C})
    a = annotate_sample(tmp_path, 'S', pgroup_file)
    a = a[a['locus'] == 'A'].set_index('allele_position')
    assert a.loc['allele1', 'exon2_depth'] == 212.4
    assert a.loc['allele1', 'call_status'] == 'Frequency-resolved'
    assert 'Low exon depth' not in a.loc['allele1', 'flags']


def test_all_loci_homozygous_triggers_qc_review(pgroup_file, tmp_path):
    ests = {g: _est([f"HLA-{g}*{x}\t-\t{COV_GOOD}\t-"]) for g, x in
            (('A', '02:01:01:01'), ('B', '07:02:01:01'), ('C', '07:02:01:01'))}
    final = "A\tHLA-A*02:01:01\t-\nB\tHLA-B*07:02:01\t-\nC\tHLA-C*07:02:01\t-\n"
    _write(tmp_path, 'S', final, ests)
    df = annotate_sample(tmp_path, 'S', pgroup_file)
    assert set(df['call_status']) == {'Unique;QC review'}
    assert set(df['flags']) == {'All loci homozygous'}


def test_two_homozygous_loci_are_normal(pgroup_file, tmp_path):
    ests = {g: _est([f"HLA-{g}*{x}\t-\t{COV_GOOD}\t-"]) for g, x in
            (('A', '02:01:01:01'), ('B', '07:02:01:01'))}
    ests['C'] = HET_C
    final = "A\tHLA-A*02:01:01\t-\nB\tHLA-B*07:02:01\t-\nC\tHLA-C*07:02:01\tHLA-C*05:01:01\n"
    _write(tmp_path, 'S', final, ests)
    df = annotate_sample(tmp_path, 'S', pgroup_file)
    assert set(df['call_status']) == {'Unique'} and set(df['flags']) == {''}


def test_missing_alignment_file_degrades_locus_not_sample(pgroup_file, tmp_path, capsys):
    rows = [f"HLA-A*02:01:01:01\tHLA-A*03:01:01:01\t{COV_GOOD}\t{COV_GOOD}"]
    _write(tmp_path, 'S', "A\tHLA-A*02:01:01\tHLA-A*03:01:01\n" + HET_FINAL_BC,
           {'A': _est(rows, other=["HLA-A*02:02:01\tHLA-A*03:02:01"]), 'B': HET_B, 'C': HET_C})
    empty_align = tmp_path / 'align'
    empty_align.mkdir()
    df = annotate_sample(tmp_path, 'S', pgroup_file, alignment_dir=empty_align)
    assert df[df['locus'] == 'A']['phase_check'].iloc[0] == 'not_evaluated'
    assert 'no IMGT alignment A_nuc.txt' in capsys.readouterr().err
