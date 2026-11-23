from pathlib import Path

import pytest

from hlahd_annotate.alignment import load_typed_exons
from hlahd_annotate.reporting import classify, two_field

ALIGN = Path(__file__).parent / "fixtures" / "align"
PG = {'A*02:01': 'A*02:01P', 'A*02:02': 'A*02:02P', 'A*02:99': 'A*02:01P',
      'A*03:01': 'A*03:01P', 'A*03:02': 'A*03:02P'}
GOOD = {'exon2_depth': 200.0, 'exon3_depth': 210.0,
        'exon2_incomplete': False, 'exon3_incomplete': False}
P1 = ('A*02:01:01:01', 'A*03:01:01:01')
P2 = ('A*02:02:01', 'A*03:02:01')


def est(slots, other=()):
    """Best pairs as (slot1, slot2); a slot is one allele name or a list of them."""
    as_list = lambda s: [s] if isinstance(s, str) else list(s)  # noqa: E731
    return {'best_pair_slots': [(as_list(a), as_list(b)) for a, b in slots],
            'other_pairs': list(other)}


def run(final, e, typed=None, cov=(GOOD, GOOD), depth=75):
    return classify(final, e, PG, typed, list(cov), depth)


@pytest.fixture
def typed():
    return load_typed_exons(ALIGN, 'A', set(P1) | set(P2) | {'A*02:99:01'})


@pytest.mark.parametrize('allele,expected', [
    ('HLA-A*01:01:01:02N', 'A*01:01'), ('A*24:09N', 'A*24:09N'),
    ('A*02:01:01', 'A*02:01'), ('A*02:01', 'A*02:01')])
def test_two_field(allele, expected):
    assert two_field(allele) == expected


def test_expression_variants_in_a_slot_stay_unique():
    # Real est.txt shape: one slot lists A*02:01:01:01, ...:02L, ...:Q variants.
    slot = ['A*02:01:01:01', 'A*02:01:01:02L', 'A*02:01:01:03Q']
    r = run([('A*02:01:01', 'A*03:01:01')], est([(slot, ['A*03:01:01:01'])]))
    assert r['call_status'] == ['Unique']


def test_unique():
    r = run([P1], est([P1]))
    assert r['call_status'] == ['Unique'] and r['flags'] == []
    assert r['two_field_genotype'] == 'A*02:01/A*03:01' and r['provenance'] == 'est'


def test_best_pairs_differing_only_past_two_fields_are_unique():
    # '#Best allele pair 2' but both collapse to A*02:01/A*03:01: raw flag would say multiple
    r = run([P1], est([P1, ('A*02:01:01:02', 'A*03:01:01:102')]))
    assert r['call_status'] == ['Unique']
    assert r['frequency_resolution_details'] is None


def test_frequency_resolved_lists_tied_pairs():
    r = run([P1], est([P1, P2]))
    assert r['call_status'] == ['Frequency-resolved']
    assert r['frequency_resolution_details'] == 'A*02:01/A*03:01; A*02:02/A*03:02'
    assert 'Frequency tie' not in r['flags']


def test_multiple_best_pairs_with_frequency_tie():
    r = run([P1, P2], est([P1, P2]))
    assert r['call_status'] == ['Multiple best pairs']
    assert 'Frequency tie' in r['flags']
    assert r['frequency_resolution_details'] is None  # nothing was dropped


def test_partial_frequency_resolution():
    third = ('A*02:01:01:01', 'A*03:02:01')
    r = run([P1, P2], est([P1, P2, third]))
    assert r['call_status'] == ['Multiple best pairs'] and 'Frequency tie' in r['flags']
    assert r['frequency_resolution_details'].count('/') == 3


def test_failed_is_exclusive():
    r = run([(None, None)], est([]), cov=())
    assert r['call_status'] == ['Failed'] and r['flags'] == []


def test_homozygous_slot_pairs_with_itself():
    r = run([('A*02:01:01:01', 'A*02:01:01:01')], est([(['A*02:01:01:01', 'A*02:01:01:02'], [])]))
    assert r['call_status'] == ['Unique']
    assert r['two_field_genotype'] == 'A*02:01/A*02:01'


def test_homozygous_slot_never_invents_heterozygous_candidates():
    # Slot lists two 2-field types; homozygous means each pairs only with itself.
    r = run([('A*02:01:01:01', 'A*02:01:01:01')],
            est([(['A*02:01:01:01', 'A*02:99:01'], [])]))
    assert r['call_status'] == ['Frequency-resolved']
    assert r['frequency_resolution_details'] == 'A*02:01/A*02:01; A*02:99/A*02:99'


def test_phase_ambiguity_from_other_pair(typed):
    # HLA-HD's '#Other ambiguous pair' is the phase alternative; sequence test confirms.
    r = run([P1], est([P1], other=[P2]), typed=typed)
    assert r['call_status'] == ['Unique']
    assert r['flags'] == ['Phase ambiguity'] and r['phase_check'] == 'evaluated'


def test_other_pair_identical_over_typed_exons_is_not_flagged(typed):
    r = run([P1], est([P1], other=[('A*02:99:01', 'A*03:01:01:01')]), typed=typed)
    assert r['flags'] == [] and r['phase_check'] == 'evaluated'


def test_other_pair_with_real_base_differences_is_not_phase(typed):
    r = run([P1], est([P1], other=[('A*02:02:01', 'A*03:01:01:01')]), typed=typed)
    assert 'Phase ambiguity' not in r['flags']
    assert r['other_pair_relation'] == 'different'  # reported, so the HTML never says 'identical'


def test_other_pair_relations_are_recorded_per_row(typed):
    r = run([P1], est([P1], other=[P2, ('A*02:99:01', 'A*03:01:01:01')]), typed=typed)
    assert r['other_pair_relation'] == 'phase;identical'


def test_unparsed_other_rows_are_flagged_incomplete_not_dropped(typed):
    e = {**est([P1]), 'other_rows_unparsed': 1}
    r = run([P1], e, typed=typed)
    assert r['flags'] == ['Phase ambiguity'] and r['phase_check'] == 'incomplete'
    assert r['other_pair_relation'] == 'untestable'


def test_unparsed_other_rows_without_alignment_are_not_evaluated():
    r = run([P1], {**est([P1]), 'other_rows_unparsed': 2})
    assert r['flags'] == ['Phase ambiguity'] and r['phase_check'] == 'not_evaluated'


def test_other_pair_never_sets_shared_pgroup(typed):
    r = run([P1], est([P1], other=[('A*02:99:01', 'A*03:01:01:01')]), typed=typed)
    assert 'Shared P-group' not in r['flags']


def test_shared_pgroup_between_read_tied_best_pairs(typed):
    r = run([P1], est([P1, ('A*02:99:01', 'A*03:01:01:01')]), typed=typed)
    assert r['call_status'] == ['Frequency-resolved'] and r['flags'] == ['Shared P-group']


def test_read_tied_best_pairs_never_set_phase_ambiguity(typed):
    r = run([P1], est([P1, P2]), typed=typed)
    assert r['call_status'] == ['Frequency-resolved'] and 'Phase ambiguity' not in r['flags']
    assert r['phase_check'] == 'no_other_pair'


def test_other_pair_without_alignment_is_flagged_unverified():
    r = run([P1], est([P1], other=[P2]))
    assert r['flags'] == ['Phase ambiguity'] and r['phase_check'] == 'not_evaluated'


def test_other_pair_allele_missing_from_alignment_is_flagged_incomplete(typed):
    r = run([P1], est([P1], other=[('A*02:77:01', 'A*03:01:01:01')]), typed=typed)
    assert r['flags'] == ['Phase ambiguity'] and r['phase_check'] == 'incomplete'


def test_other_pair_every_column_unknown_is_flagged_incomplete():
    unknown = {a: '??????????' for a in (*P1, *P2)}
    r = run([P1], est([P1], other=[P2]), typed=unknown)
    assert r['flags'] == ['Phase ambiguity'] and r['phase_check'] == 'incomplete'


def test_no_other_pair(typed):
    assert run([P1], est([P1]), typed=typed)['phase_check'] == 'no_other_pair'


def test_low_exon_depth_triggers_qc_review():
    low = {**GOOD, 'exon3_depth': 74.9}
    r = run([P1], est([P1]), cov=(GOOD, low))
    assert r['call_status'] == ['Unique', 'QC review'] and r['flags'] == ['Low exon depth']


def test_depth_at_threshold_passes():
    r = run([P1], est([P1]), cov=(GOOD, {**GOOD, 'exon2_depth': 75.0}))
    assert r['call_status'] == ['Unique']


def test_incomplete_coverage_triggers_qc_review():
    r = run([P1], est([P1]), cov=(GOOD, {**GOOD, 'exon2_incomplete': True}))
    assert r['call_status'] == ['Unique', 'QC review']
    assert r['flags'] == ['Incomplete exon coverage']


def test_unknown_depth_is_qc_review_not_a_pass():
    r = run([P1], est([P1]), cov=(GOOD, {**GOOD, 'exon2_depth': None}))
    assert r['call_status'] == ['Unique', 'QC review'] and r['qc_evaluable'] is False


def test_missing_est_uses_final_only():
    r = run([P1], {}, cov=())
    assert r['provenance'] == 'final_only'
    assert r['call_status'] == ['Unique', 'QC review'] and r['qc_evaluable'] is False


def _table(missing: dict, exons=range(1, 9)):
    from hlahd_annotate.dictionary import ExonTable
    t = ExonTable(missing)
    t.exons = tuple(exons)
    return t


EXON_2_3_ONLY = [1, 4, 5, 6, 7, 8]


def test_phase_alternative_with_other_exons_is_ruled_out_not_flagged(typed):
    # A*02:02:01 has exon 4: HLA-HD scored it there and it lost -> grey note, no flag.
    table = _table({'A*02:02:01': [1, 5, 6, 7, 8], 'A*03:02:01': EXON_2_3_ONLY})
    r = classify([P1], est([P1], other=[P2]), PG, typed, [GOOD, GOOD], 75, table)
    assert 'Phase ambiguity' not in r['flags']
    assert r['other_pair_relation'] == 'phase' and r['other_pair_evidence'] == 'ruled_out'


def test_phase_alternative_known_only_over_exons_2_3_stays_flagged(typed):
    table = _table({'A*02:02:01': EXON_2_3_ONLY, 'A*03:02:01': EXON_2_3_ONLY})
    r = classify([P1], est([P1], other=[P2]), PG, typed, [GOOD, GOOD], 75, table)
    assert r['flags'] == ['Phase ambiguity'] and r['other_pair_evidence'] == 'unresolved'


def test_phase_alternative_missing_from_dictionary_stays_flagged(typed):
    table = _table({'A*03:02:01': EXON_2_3_ONLY})  # A*02:02:01 absent
    r = classify([P1], est([P1], other=[P2]), PG, typed, [GOOD, GOOD], 75, table)
    assert r['flags'] == ['Phase ambiguity'] and r['other_pair_evidence'] == 'not_evaluated'


def test_no_dictionary_keeps_phase_flag(typed):
    r = run([P1], est([P1], other=[P2]), typed=typed)
    assert r['flags'] == ['Phase ambiguity'] and r['other_pair_evidence'] == 'not_evaluated'


def test_one_unresolved_row_keeps_flag_even_if_another_is_ruled_out(typed):
    other2 = ('A*02:02:01', 'A*03:02:01')
    table = _table({'A*02:02:01': EXON_2_3_ONLY, 'A*03:02:01': EXON_2_3_ONLY})
    e = {**est([P1], other=[P2, other2]), 'other_rows_unparsed': 0}
    table2 = _table({**table, 'A*02:02:01': [1, 5, 6, 7, 8]})
    assert classify([P1], e, PG, typed, [GOOD, GOOD], 75, table2)['flags'] == []
    assert 'Phase ambiguity' in classify([P1], e, PG, typed, [GOOD, GOOD], 75, table)['flags']


def test_unparsed_row_is_never_ruled_out(typed):
    table = _table({'A*02:02:01': [1, 5, 6, 7, 8], 'A*03:02:01': EXON_2_3_ONLY})
    e = {**est([P1], other=[P2]), 'other_rows_unparsed': 1}
    r = classify([P1], e, PG, typed, [GOOD, GOOD], 75, table)
    assert r['flags'] == ['Phase ambiguity']
    assert r['other_pair_evidence'] == 'not_evaluated;ruled_out'


def test_locus_exons_come_from_the_header():
    # HLA-B has exons 1-7: an allele missing 1,4-7 is exon 2-3 only.
    from hlahd_annotate.dictionary import known_beyond_exons_2_3
    t = _table({'B*40:21': [1, 4, 5, 6, 7], 'B*15:01:01:01': []}, exons=range(1, 8))
    assert known_beyond_exons_2_3('HLA-B*40:21', t) is False
    assert known_beyond_exons_2_3('B*15:01:01:01', t) is True
    assert known_beyond_exons_2_3('B*99:99', t) is None
