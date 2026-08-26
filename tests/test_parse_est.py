from hlahd_annotate.parse_est import parse_est, cross_check_allele, _parse_coverage


# --- parse_est tests ---

def test_sample1_has_ambiguous_pair(sample1_est):
    result = parse_est(sample1_est)
    assert result['multiple_best_pairs'] is False
    assert result['has_ambiguous_pair'] is True


def test_sample1_ambiguous_alleles_recorded(sample1_est):
    result = parse_est(sample1_est)
    assert result['ambiguous_alleles'] is not None
    assert 'A*02:614:01' in result['ambiguous_alleles']
    assert 'A*68:164:02' in result['ambiguous_alleles']


def test_sample1_best_pair_alleles_collected(sample1_est):
    result = parse_est(sample1_est)
    assert len(result['best_pair_alleles']) > 0
    # All should have HLA- prefix stripped
    assert all(not a.startswith('HLA-') for a in result['best_pair_alleles'])
    # Spot-check a known allele from sample1_A.est.txt
    assert 'A*02:01:01:01' in result['best_pair_alleles']
    assert 'A*68:01:01:01' in result['best_pair_alleles']


def test_sample2_multiple_best_pairs(sample2_est):
    result = parse_est(sample2_est)
    assert result['multiple_best_pairs'] is True
    assert result['has_ambiguous_pair'] is False


def test_sample3_normal(sample3_est):
    result = parse_est(sample3_est)
    assert result['multiple_best_pairs'] is False
    assert result['has_ambiguous_pair'] is False
    assert result['ambiguous_alleles'] is None


def test_sample4_dash_allele2_not_collected(sample4_est):
    result = parse_est(sample4_est)
    assert '-' not in result['best_pair_alleles']
    assert len(result['best_pair_alleles']) > 0


# --- cross_check_allele tests ---

def test_cross_check_match_3field():
    best_alleles = ['A*02:01:01:01', 'A*02:01:01:198', 'A*68:01:01:01']
    match, detail = cross_check_allele('HLA-A*02:01:01', best_alleles)
    assert match is True
    assert detail is None


def test_cross_check_match_2field():
    best_alleles = ['A*02:01:01:01', 'A*68:01:01:01']
    match, detail = cross_check_allele('HLA-A*02:01', best_alleles)
    assert match is True


def test_cross_check_mismatch():
    best_alleles = ['A*03:01:01:01', 'A*24:02:01:01']
    match, detail = cross_check_allele('HLA-A*02:01:01', best_alleles)
    assert match is False
    assert detail is not None
    assert 'HLA-A*02:01:01' in detail


def test_cross_check_null_allele():
    match, detail = cross_check_allele(None, ['A*02:01:01:01'])
    assert match is True


def test_cross_check_empty_best_alleles():
    match, detail = cross_check_allele('HLA-A*02:01:01', [])
    assert match is True  # no est data to contradict


# --- _parse_coverage tests ---

def test_parse_coverage_complete():
    cov = _parse_coverage('exon2:262.426:comp.0,exon3:303.069:comp.0')
    assert cov['exon2_depth'] == 262.426
    assert cov['exon2_incomp'] == 0
    assert cov['exon3_depth'] == 303.069
    assert cov['exon3_incomp'] == 0


def test_parse_coverage_incomplete():
    cov = _parse_coverage('exon2:45.231:comp.0,exon3:51.002:incomp.62')
    assert cov['exon3_depth'] == 51.002
    assert cov['exon3_incomp'] == 62
    assert cov['exon2_incomp'] == 0


def test_parse_coverage_null_field():
    cov = _parse_coverage('-')
    assert cov['exon2_depth'] is None
    assert cov['exon2_incomp'] is None
    assert cov['exon3_depth'] is None
    assert cov['exon3_incomp'] is None


# --- best_pairs tests ---

def test_best_pairs_complete(sample3_est):
    result = parse_est(sample3_est)
    assert len(result['best_pairs']) == 1
    a1 = result['best_pairs'][0]['allele1_coverage']
    assert a1['exon2_depth'] == 262.426
    assert a1['exon2_incomp'] == 0


def test_best_pairs_incomplete(sample5_est):
    result = parse_est(sample5_est)
    bp = result['best_pairs'][0]
    assert bp['allele1_coverage']['exon3_incomp'] == 7
    assert bp['allele2_coverage']['exon2_incomp'] == 4


def test_best_pairs_multiple(sample2_est):
    result = parse_est(sample2_est)
    assert len(result['best_pairs']) == 2
    assert result['best_pairs'][1]['allele1_coverage']['exon2_depth'] is not None


def test_best_pairs_dash_side(sample4_est):
    # sample4 allele2 column is '-' (null), but its coverage column is present
    result = parse_est(sample4_est)
    assert len(result['best_pairs']) == 1
    assert result['best_pairs'][0]['allele1_coverage']['exon2_depth'] == 262.426
    assert result['best_pairs'][0]['allele2_coverage']['exon2_depth'] == 257.163


def test_best_pairs_representative_names(sample3_est):
    bp = parse_est(sample3_est)['best_pairs'][0]
    assert bp['allele1'] == 'A*02:01:01:01'
    assert bp['allele2'] == 'A*03:01:01:01'


def test_best_pairs_representative_name_dash_is_none(sample4_est):
    bp = parse_est(sample4_est)['best_pairs'][0]
    assert bp['allele1'] == 'A*02:01:01:01'
    assert bp['allele2'] is None
