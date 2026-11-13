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
    r = run([P1], est([P1], other=[P2]), typed=typed)
    assert r['call_status'] == ['Unique']
    assert r['flags'] == ['Phase ambiguity'] and r['phase_check'] == 'evaluated'


def test_identical_other_pair_is_shared_pgroup_not_phase(typed):
    r = run([P1], est([P1], other=[('A*02:99:01', 'A*03:01:01:01')]), typed=typed)
    assert r['flags'] == ['Shared P-group']


def test_phase_not_evaluated_without_alignment():
    r = run([P1], est([P1], other=[P2]))
    assert r['phase_check'] == 'not_evaluated' and 'Phase ambiguity' not in r['flags']


def test_phase_incomplete_when_allele_missing_from_alignment(typed):
    r = run([P1], est([P1], other=[('A*02:77:01', 'A*03:01:01:01')]), typed=typed)
    assert r['phase_check'] == 'incomplete'


def test_no_competing_pairs(typed):
    assert run([P1], est([P1]), typed=typed)['phase_check'] == 'no_competing'


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
