import pytest
import pandas as pd
from hlahd_annotate.parse_final import parse_final_result, CLASS_I_LOCI

# sample3_final.result.txt:
#   A  HLA-A*02:01:01  HLA-A*03:01:01          (3 cols, normal)
#   B  HLA-B*18:11  HLA-B*50:58  HLA-B*18:164  HLA-B*50:58  (5 cols, multiple pairs)
#   C  HLA-C*12:02:02  HLA-C*17:01:01          (3 cols, normal)
#   DRB1, DQA1, DQB1, DPA1, DPB1               (non-class-I, skipped)
#
# samplename2_final.result.txt:
#   A  HLA-A*02:01:01  HLA-A*03:01:01          (3 cols, normal)
#   B  HLA-B*52:01:01  HLA-B*41:01:01          (3 cols, normal)
#   C  HLA-C*12:02:02  HLA-C*17:01:01          (3 cols, normal)
#   DRB1..V                                     (many non-class-I, skipped)


def test_parse_allele_values_sample3(sample3_final):
    df = parse_final_result(sample3_final)
    a_alleles = df[df['locus'] == 'A']['allele'].tolist()
    assert 'HLA-A*02:01:01' in a_alleles
    assert 'HLA-A*03:01:01' in a_alleles


def test_parse_allele_positions_normal(samplename2_final):
    # samplename2 has all normal rows — only allele1 and allele2 positions
    df = parse_final_result(samplename2_final)
    positions = set(df['allele_position'].unique())
    assert positions == {'allele1', 'allele2'}


def test_parse_resolution_3field(sample3_final):
    # A*02:01:01 is 3-field
    df = parse_final_result(sample3_final)
    a1 = df[(df['locus'] == 'A') & (df['allele_position'] == 'allele1')].iloc[0]
    assert a1['resolution'] == 3


def test_parse_resolution_2field(sample3_final):
    # B*18:11 is 2-field (in the multiple-pair row)
    df = parse_final_result(sample3_final)
    b1 = df[(df['locus'] == 'B') & (df['allele_position'] == 'allele1')].iloc[0]
    assert b1['resolution'] == 2


def test_parse_not_typed_is_null(tmp_path):
    f = tmp_path / "s_final.result.txt"
    f.write_text("A\tNot typed\tNot typed\n")
    df = parse_final_result(f)
    assert df[(df['locus'] == 'A')]['allele'].isna().all()


def test_parse_dash_allele2_is_null(tmp_path):
    f = tmp_path / "s_final.result.txt"
    f.write_text("A\tHLA-A*02:01:01\t-\n")
    df = parse_final_result(f)
    a2 = df[(df['locus'] == 'A') & (df['allele_position'] == 'allele2')].iloc[0]
    assert a2['allele'] is None or pd.isna(a2['allele'])


def test_parse_multiple_best_pairs_sample3(sample3_final):
    # sample3_final B row has 5 columns (2 pairs)
    df = parse_final_result(sample3_final)
    b_rows = df[df['locus'] == 'B']
    assert b_rows['multiple_best_pairs'].all()
    b_positions = set(b_rows['allele_position'].tolist())
    assert 'allele1_pair2' in b_positions
    assert 'allele2_pair2' in b_positions


def test_parse_multiple_pairs_row_count(sample3_final):
    # B with 5 cols produces 4 position rows
    df = parse_final_result(sample3_final)
    b_rows = df[df['locus'] == 'B']
    assert len(b_rows) == 4


def test_parse_normal_has_no_multiple_pairs(samplename2_final):
    # samplename2 has all normal rows
    df = parse_final_result(samplename2_final)
    assert not df['multiple_best_pairs'].any()


def test_parse_skips_non_class_i_sample3(sample3_final):
    # sample3_final contains DRB1, DQA1, DQB1, DPA1, DPB1
    df = parse_final_result(sample3_final)
    assert set(df['locus'].unique()).issubset(CLASS_I_LOCI)


def test_parse_skips_non_class_i_samplename2(samplename2_final):
    # samplename2_final contains many non-class-I loci including E, F, G, H, J, K, L, V
    df = parse_final_result(samplename2_final)
    assert set(df['locus'].unique()).issubset(CLASS_I_LOCI)


def test_parse_samplename2_row_count(samplename2_final):
    # All 3 loci normal: 2 alleles x 3 loci = 6 rows
    df = parse_final_result(samplename2_final)
    assert len(df) == 6


def test_parse_sample3_total_rows(sample3_final):
    # A: 2, B: 4 (multiple pairs), C: 2 = 8 rows
    df = parse_final_result(sample3_final)
    assert len(df) == 8


def test_parse_required_columns(sample3_final):
    df = parse_final_result(sample3_final)
    assert set(df.columns) == {'locus', 'allele', 'allele_position', 'resolution', 'multiple_best_pairs'}


def test_parse_null_allele_has_zero_resolution(tmp_path):
    f = tmp_path / "s_final.result.txt"
    f.write_text("A\tHLA-A*02:01:01\tNot typed\n")
    df = parse_final_result(f)
    null_row = df[(df['locus'] == 'A') & (df['allele_position'] == 'allele2')].iloc[0]
    assert null_row['resolution'] == 0
