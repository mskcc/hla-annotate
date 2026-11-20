"""Incomplete Allele: exon presence from HLA-HD's dictionary <GENE>/<GENE>.list."""
from pathlib import Path

import pytest

from hlahd_annotate.annotate import annotate_sample
from hlahd_annotate.dictionary import fmt_exons, load_exon_presence, missing_exons

DICT = Path(__file__).parent / "fixtures" / "dict"
COV = 'exon2:212.4:comp.0,exon3:198.7:comp.0'


@pytest.fixture
def table():
    return load_exon_presence(DICT, 'A')


def test_parses_exon_columns(table):
    assert table['A*02:01:01:01'] == []
    assert table['A*24:05:01'] == [1, 5, 6, 7]


def test_exact_name(table):
    assert missing_exons('HLA-A*24:05:01', table) == [1, 5, 6, 7]


def test_prefix_all_entries_incomplete_returns_most_complete(table):
    # A*68:99:01 covers :01 (missing 1,5-8) and :02 (missing 7,8)
    assert missing_exons('A*68:99:01', table) == [7, 8]


def test_prefix_with_one_complete_entry_is_complete(table):
    assert missing_exons('A*11:77:01', table) == []


def test_prefix_matches_field_boundary_only(table):
    assert missing_exons('A*24:0', table) is None


def test_not_in_dictionary(table):
    assert missing_exons('A*99:01', table) is None


@pytest.mark.parametrize('exons,text', [([1, 5, 6, 7], '1, 5-7'), ([7, 8], '7-8'), ([4], '4'),
                                        ([1, 4, 5, 6, 7, 8], '1, 4-8')])
def test_fmt_exons(exons, text):
    assert fmt_exons(exons) == text


def test_list_without_header_raises(tmp_path):
    (tmp_path / 'A').mkdir()
    (tmp_path / 'A' / 'A.list').write_text("HLA00005\tHLA-A*02:01:01:01\t+\n")
    with pytest.raises(ValueError):
        load_exon_presence(tmp_path, 'A')


def _sample(tmp_path, a1, a2):
    """A-locus-only sample (B and C missing -> Failed; irrelevant here)."""
    (tmp_path / 'S_final.result.txt').write_text(f"A\tHLA-{a1}\tHLA-{a2}\n")
    (tmp_path / 'S_A.est.txt').write_text(
        f"#Pair count\t1\n#Best allele pair\t1\nHLA-{a1}\tHLA-{a2}\t{COV}\t{COV}\n")
    return tmp_path


def _locus_a(df):
    return df[df['locus'] == 'A'].iloc[0]


def test_incomplete_reported_allele_is_flagged(tmp_path, pgroup_file):
    df = annotate_sample(_sample(tmp_path, 'A*02:01:01:01', 'A*24:05:01'), 'S', pgroup_file,
                         dictionary_dir=DICT)
    a = _locus_a(df)
    assert 'Incomplete Allele' in a['flags'].split(';')
    assert a['incomplete_allele_detail'] == 'A*24:05:01 missing exons 1, 5-7'
    assert 'QC review' not in a['call_status']


def test_complete_alleles_not_flagged(tmp_path, pgroup_file):
    a = _locus_a(annotate_sample(_sample(tmp_path, 'A*02:01:01:01', 'A*03:01:01:01'), 'S',
                                 pgroup_file, dictionary_dir=DICT))
    assert 'Incomplete Allele' not in a['flags'] and a['incomplete_allele_detail'] == ''


def test_allele_missing_from_dictionary_is_named_not_passed(tmp_path, pgroup_file):
    a = _locus_a(annotate_sample(_sample(tmp_path, 'A*02:01:01:01', 'A*99:01:01'), 'S',
                                 pgroup_file, dictionary_dir=DICT))
    assert 'Incomplete Allele' not in a['flags']
    assert a['incomplete_allele_detail'] == 'A*99:01:01 not in HLA-HD dictionary'


def test_no_dictionary_is_not_evaluated(tmp_path, pgroup_file):
    a = _locus_a(annotate_sample(_sample(tmp_path, 'A*02:01:01:01', 'A*24:05:01'), 'S', pgroup_file))
    assert a['incomplete_allele_detail'] == 'not evaluated' and 'Incomplete Allele' not in a['flags']


def test_missing_locus_list_warns_and_is_not_evaluated(tmp_path, pgroup_file, capsys):
    df = annotate_sample(_sample(tmp_path, 'A*02:01:01:01', 'A*24:05:01'), 'S', pgroup_file,
                         dictionary_dir=tmp_path / 'nodict')
    assert _locus_a(df)['incomplete_allele_detail'] == 'not evaluated'
    assert 'no HLA-HD dictionary list A/A.list' in capsys.readouterr().err
