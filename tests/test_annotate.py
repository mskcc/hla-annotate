import shutil
import pandas as pd
from hlahd_annotate.annotate import annotate_sample

EXPECTED_COLUMNS = [
    'sample', 'locus', 'allele', 'allele_position', 'resolution',
    'p_group', 'p_group_found', 'multiple_best_pairs', 'has_ambiguous_pair',
    'ambiguous_alleles', 'est_mismatch', 'est_mismatch_detail',
    'exon2_depth', 'exon2_incomp', 'exon3_depth', 'exon3_incomp',
    'incomplete_coverage',
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


COVERAGE_COLUMNS = ['exon2_depth', 'exon2_incomp', 'exon3_depth', 'exon3_incomp',
                    'incomplete_coverage']


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
