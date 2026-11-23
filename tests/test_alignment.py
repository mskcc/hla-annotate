from pathlib import Path

import pytest

from hlahd_annotate.alignment import load_typed_exons, phase_relation

ALIGN = Path(__file__).parent / "fixtures" / "align"
ALL = {'A*01:01:01:01', 'A*02:01:01:01', 'A*03:01:01:01', 'A*02:02:01',
       'A*03:02:01', 'A*02:99:01', 'A*02:50'}


@pytest.fixture
def ex():
    return load_typed_exons(ALIGN, 'A', ALL)


def test_exons_2_and_3_only_across_blocks(ex):
    # exon2 ACGTAC + exon3 GGCC; exon1 ATG and exon4 TT excluded
    assert ex['A*01:01:01:01'] == 'ACGTACGGCC'


def test_dash_is_reference_base(ex):
    assert ex['A*02:01:01:01'] == 'ATGTGCGGCC'


def test_star_is_unknown(ex):
    assert ex['A*02:50'] == '???TGCGGCC'


def test_difference_outside_typed_exons_is_invisible(ex):
    assert ex['A*02:99:01'] == ex['A*02:01:01:01']


def test_missing_allele_left_out_and_prefix_matches_field_boundary():
    out = load_typed_exons(ALIGN, 'A', {'A*02:02', 'A*02:0', 'A*99:01'})
    assert set(out) == {'A*02:02'}  # 'A*02:0' must not match A*02:02:01


def test_file_without_allele_rows_raises(tmp_path):
    (tmp_path / 'A_nuc.txt').write_text("# header only\n cDNA 1\n")
    with pytest.raises(ValueError):
        load_typed_exons(tmp_path, 'A', {'A*01:01'})


def test_phase_conjugate(ex):
    x = (ex['A*02:01:01:01'], ex['A*03:01:01:01'])
    y = (ex['A*02:02:01'], ex['A*03:02:01'])
    assert phase_relation(x, y) == ('phase', 0)


def test_identical_over_typed_exons_is_not_phase(ex):
    x = (ex['A*02:01:01:01'], ex['A*03:01:01:01'])
    y = (ex['A*03:01:01:01'], ex['A*02:99:01'])  # swapped slot order
    assert phase_relation(x, y) == ('identical', 0)


def test_different_genotype(ex):
    x = (ex['A*02:01:01:01'], ex['A*03:01:01:01'])
    y = (ex['A*02:02:01'], ex['A*03:01:01:01'])
    assert phase_relation(x, y)[0] == 'different'


def test_unknown_columns_are_counted(ex):
    x = (ex['A*02:01:01:01'], ex['A*03:01:01:01'])
    y = (ex['A*02:50'], ex['A*03:01:01:01'])
    assert phase_relation(x, y) == ('identical', 3)


def test_row_longer_than_reference_is_kept_not_an_error(tmp_path):
    # IMGT expression variants (e.g. C*04:09L) read past the reference stop codon.
    (tmp_path / 'A_nuc.txt').write_text(
        " A*01:01:01:01     ATG|ACG TAC|GGC C|TGA\n"
        " A*03:437Q         ---|--- ---|--- -|--C AAT\n")
    ex = load_typed_exons(tmp_path, 'A', {'A*03:437Q'})
    assert ex['A*03:437Q'] == 'ACGTACGGCC'


def test_row_shorter_than_reference_ends_unknown(tmp_path):
    (tmp_path / 'A_nuc.txt').write_text(
        " A*01:01:01:01     ATG|ACG TAC|GGC C|TT\n"
        " A*02:99           ---|--- ---|--\n")
    assert load_typed_exons(tmp_path, 'A', {'A*02:99'})['A*02:99'] == 'ACGTACGG??'
