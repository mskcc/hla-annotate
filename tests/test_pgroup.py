from hlahd_annotate.pgroup import load_pgroup_table, lookup_pgroup


def test_load_pgroup_table_returns_dict(pgroup_file):
    result = load_pgroup_table(pgroup_file)
    assert isinstance(result, dict)
    assert len(result) > 0


def test_load_pgroup_skips_comment_lines(pgroup_file):
    result = load_pgroup_table(pgroup_file)
    for key in result:
        assert not key.startswith('#')


def test_load_pgroup_multi_field_keys(pgroup_file):
    # A*02:01:01:01 is in fixture; lookup table should have 2-field, 3-field, 4-field keys
    result = load_pgroup_table(pgroup_file)
    assert 'A*02:01' in result
    assert 'A*02:01:01' in result
    assert 'A*02:01:01:01' in result


def test_load_pgroup_all_map_to_correct_pgroup(pgroup_file):
    result = load_pgroup_table(pgroup_file)
    # All three field-depth variants should map to same P-group
    assert result['A*02:01'] == 'A*02:01P'
    assert result['A*02:01:01'] == 'A*02:01P'
    assert result['A*02:01:01:01'] == 'A*02:01P'


def test_lookup_pgroup_3field_allele(pgroup_file):
    lookup = load_pgroup_table(pgroup_file)
    p_group, found = lookup_pgroup('HLA-A*02:01:01', lookup)
    assert found is True
    assert p_group == 'A*02:01P'


def test_lookup_pgroup_2field_allele(pgroup_file):
    lookup = load_pgroup_table(pgroup_file)
    p_group, found = lookup_pgroup('HLA-A*02:01', lookup)
    assert found is True
    assert p_group == 'A*02:01P'


def test_lookup_pgroup_strips_hla_prefix(pgroup_file):
    lookup = load_pgroup_table(pgroup_file)
    with_prefix = lookup_pgroup('HLA-A*03:01:01', lookup)
    without_prefix = lookup_pgroup('A*03:01:01', lookup)
    assert with_prefix == without_prefix


def test_lookup_pgroup_not_found(pgroup_file):
    lookup = load_pgroup_table(pgroup_file)
    p_group, found = lookup_pgroup('HLA-A*99:99:99', lookup)
    assert found is False
    assert p_group is None


def test_lookup_pgroup_not_typed():
    p_group, found = lookup_pgroup('Not typed', {})
    assert found is False
    assert p_group is None


def test_lookup_pgroup_dash():
    p_group, found = lookup_pgroup('-', {})
    assert found is False
    assert p_group is None


def test_lookup_pgroup_none():
    p_group, found = lookup_pgroup(None, {})
    assert found is False
    assert p_group is None


def test_load_pgroup_excludes_no_pgroup_alleles(pgroup_file):
    # A*01:06 has no P-group in fixture (empty third field)
    # It must NOT appear in the lookup table at all
    result = load_pgroup_table(pgroup_file)
    assert 'A*01:06' not in result
    # Also verify no garbage 'A*' value was inserted
    for key, value in result.items():
        assert value != 'A*', f"Garbage p_group 'A*' found for key {key!r}"


def test_load_pgroup_strips_allele_modifier(pgroup_file):
    # A*01:01:38L in fixture — the 'L' modifier must be stripped
    # so A*01:01:38 (not A*01:01:38L) is the key in the lookup table
    result = load_pgroup_table(pgroup_file)
    assert 'A*01:01:38' in result
    assert result['A*01:01:38'] == 'A*01:01P'
    # The unstripped key must NOT be in the table
    assert 'A*01:01:38L' not in result
