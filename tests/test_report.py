from pathlib import Path
import pandas as pd
from hlahd_annotate.annotate import annotate_sample
from hlahd_annotate.report import render_report

TEMPLATE_DIR = Path(__file__).parents[1] / "src" / "hlahd_annotate" / "templates"
QC = dict(call_status='Unique;QC review', flags='Incomplete exon coverage')


def _make_row(
    sample='S1', locus='A', allele='HLA-A*02:01:01', allele_position='allele1',
    resolution=3, p_group='A*02:01P', p_group_found=True,
    multiple_best_pairs=False, has_ambiguous_pair=False, ambiguous_alleles=None,
    est_mismatch=False, est_mismatch_detail=None,
    exon2_depth=262.4, exon2_incomp=0, exon3_depth=303.0, exon3_incomp=0,
    incomplete_coverage=False, exon2_incomplete=None, exon3_incomplete=None,
    incomplete_bp=None,
    call_status='Unique', flags='', two_field_genotype='A*02:01/A*03:01',
    frequency_resolution_details=None, phase_check='no_competing',
    other_pair_relation='', other_pair_evidence='', incomplete_allele_detail='',
):
    if exon2_incomplete is None:
        exon2_incomplete = bool(exon2_incomp)
    if exon3_incomplete is None:
        exon3_incomplete = bool(exon3_incomp)
    if incomplete_bp is None:
        incomplete_bp = exon2_incomp + exon3_incomp
    return {
        'sample': sample, 'locus': locus, 'allele': allele,
        'allele_position': allele_position, 'resolution': resolution,
        'p_group': p_group, 'p_group_found': p_group_found,
        'multiple_best_pairs': multiple_best_pairs,
        'has_ambiguous_pair': has_ambiguous_pair,
        'ambiguous_alleles': ambiguous_alleles,
        'est_mismatch': est_mismatch,
        'est_mismatch_detail': est_mismatch_detail,
        'exon2_depth': exon2_depth, 'exon2_incomp': exon2_incomp,
        'exon3_depth': exon3_depth, 'exon3_incomp': exon3_incomp,
        'exon2_incomplete': exon2_incomplete, 'exon3_incomplete': exon3_incomplete,
        'incomplete_bp': incomplete_bp,
        'incomplete_coverage': incomplete_coverage,
        'call_status': call_status, 'flags': flags,
        'two_field_genotype': two_field_genotype,
        'frequency_resolution_details': frequency_resolution_details,
        'phase_check': phase_check,
        'other_pair_relation': other_pair_relation,
        'other_pair_evidence': other_pair_evidence,
        'incomplete_allele_detail': incomplete_allele_detail,
    }


def test_render_report_creates_html_file(fixtures_dir, pgroup_file, tmp_path):
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    out = render_report(df, 'sample3', str(pgroup_file), tmp_path, TEMPLATE_DIR)
    assert out.exists()
    assert out.suffix == '.html'


def test_render_report_contains_sample_id(fixtures_dir, pgroup_file, tmp_path):
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    out = render_report(df, 'sample3', str(pgroup_file), tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    assert 'sample3' in html


def test_render_report_contains_alleles(fixtures_dir, pgroup_file, tmp_path):
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    out = render_report(df, 'sample3', str(pgroup_file), tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    # Alleles in df have HLA- prefix (e.g. HLA-A*02:01:01); the substring A*02:01:01 is present
    assert 'A*02:01:01' in html
    assert 'A*03:01:01' in html


def test_render_report_contains_pgroups(fixtures_dir, pgroup_file, tmp_path):
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    out = render_report(df, 'sample3', str(pgroup_file), tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    assert 'A*02:01P' in html


def test_render_report_self_contained(fixtures_dir, pgroup_file, tmp_path):
    df = annotate_sample(fixtures_dir, 'sample3', pgroup_file)
    out = render_report(df, 'sample3', str(pgroup_file), tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    # No external stylesheet or script links (self-contained)
    assert 'href="http' not in html
    assert 'src="http' not in html


def test_render_report_no_flag_badges_for_clean_sample(tmp_path):
    rows = [
        _make_row(allele_position='allele1'),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P'),
    ]
    df = pd.DataFrame(rows)
    out = render_report(df, 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    # The note-box and CSS <style> block always contain these class names as
    # definitions; verify no badge <span> elements are actually rendered.
    assert '<span class="badge badge-yellow">' not in html
    assert '<span class="badge badge-orange">' not in html
    assert '<span class="badge badge-red">' not in html


def _pair(**model):
    return [_make_row(allele_position='allele1', **model),
            _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P', **model)]


def _html(rows, tmp_path):
    return render_report(pd.DataFrame(rows), 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR).read_text()


def test_render_report_unique_status_and_two_field(tmp_path):
    html = _html(_pair(), tmp_path)
    assert '<span class="badge badge-green">Unique</span>' in html
    assert 'A*02:01/A*03:01' in html


def test_render_report_shared_pgroup_shows_alternative_pair(tmp_path):
    html = _html(_pair(flags='Shared P-group', has_ambiguous_pair=True,
                       ambiguous_alleles='A*02:614:01 / A*68:164:02'), tmp_path)
    assert '<span class="badge badge-grey">Shared P-group</span>' in html
    assert 'A*02:614:01 / A*68:164:02' in html
    assert 'AMBIGUOUS PAIR' not in html


def test_render_report_phase_ambiguity_badge(tmp_path):
    html = _html(_pair(flags='Phase ambiguity', phase_check='evaluated'), tmp_path)
    assert '<span class="badge badge-darkred">Phase ambiguity</span>' in html


def test_render_report_phase_not_evaluated_is_stated(tmp_path):
    html = _html(_pair(phase_check='not_evaluated', flags='Phase ambiguity'), tmp_path)
    assert 'Phase ambiguity check:</strong> not evaluated' in html


def test_render_report_ruled_out_alternative_is_a_grey_note_not_a_flag(tmp_path):
    html = _html(_pair(ambiguous_alleles='HLA-C*07:19 / HLA-C*07:27:01', phase_check='evaluated',
                       other_pair_relation='phase', other_pair_evidence='ruled_out'), tmp_path)
    assert '<p class="note-grey"><strong>Alternative ruled out:</strong>' in html
    assert '[Phase ambiguity]' not in html
    assert 'identical to the reported pair' not in html  # a phase row is never called identical
    assert 'Phase ambiguity check:' not in html


def test_render_report_unresolved_alternative_keeps_flag_and_says_why(tmp_path):
    html = _html(_pair(ambiguous_alleles='HLA-C*12:98 / HLA-C*16:15:01', phase_check='evaluated',
                       flags='Phase ambiguity', other_pair_relation='phase',
                       other_pair_evidence='unresolved'), tmp_path)
    assert '[Phase ambiguity]' in html
    assert 'Not tested outside exons 2&#8211;3:</strong>' in html
    assert 'note-grey"><strong>Alternative ruled out' not in html


def test_render_report_frequency_resolved_details(tmp_path):
    html = _html(_pair(call_status='Frequency-resolved',
                       frequency_resolution_details='A*02:01/A*03:01; A*02:02/A*03:02'), tmp_path)
    assert '<span class="badge badge-yellow">Frequency-resolved</span>' in html
    assert 'A*02:02/A*03:02' in html

def test_render_report_multiple_best_pairs_badge(tmp_path):
    model = dict(call_status='Multiple best pairs', flags='Frequency tie',
                 two_field_genotype='A*02:01/A*03:01; A*24:02/A*68:01')
    rows = _pair(**model) + [
        _make_row(allele='HLA-A*68:01:01', allele_position='allele1_pair2', **model),
        _make_row(allele='HLA-A*24:02:01', allele_position='allele2_pair2', **model),
    ]
    html = _html(rows, tmp_path)
    assert '<span class="badge badge-orange">Multiple best pairs</span>' in html
    assert '<span class="badge badge-orange">Frequency tie</span>' in html
    assert 'A*24:02/A*68:01' in html

def test_render_report_est_mismatch_is_internal_not_a_flag(tmp_path):
    rows = [
        _make_row(allele_position='allele1', est_mismatch=True,
                  est_mismatch_detail='final=HLA-A*02:01:01; est_best_alleles=A*03:01,A*24:02'),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P'),
    ]
    html = _html(rows, tmp_path)
    assert 'EST MISMATCH</span>' not in html
    assert 'Internal check' in html and 'final=HLA-A*02:01:01' in html

def test_render_report_custom_hlahd_version(tmp_path):
    rows = [
        _make_row(allele_position='allele1'),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P'),
    ]
    df = pd.DataFrame(rows)
    out = render_report(df, 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR, hlahd_version='v1.9.0')
    html = out.read_text()
    assert 'v1.9.0' in html


def test_render_report_has_coverage_columns(tmp_path):
    rows = [
        _make_row(allele_position='allele1'),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P'),
    ]
    df = pd.DataFrame(rows)
    out = render_report(df, 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    assert 'Exon 2 cov' in html
    assert 'Exon 3 cov' in html
    assert '262x' in html  # depth rendered, rounded


def test_render_report_incomplete_coverage_badge(tmp_path):
    rows = [
        _make_row(allele_position='allele1', exon3_incomp=7, incomplete_coverage=True, **QC),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P',
                  exon2_incomp=4, incomplete_coverage=True, **QC),
    ]
    df = pd.DataFrame(rows)
    out = render_report(df, 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    assert '<span class="badge badge-red">Incomplete exon coverage</span>' in html
    assert '<span class="badge badge-red">QC review</span>' in html
    assert '(!)' in html  # marker on the incomplete allele's cell


def test_render_report_no_incomplete_badge_when_complete(tmp_path):
    rows = [
        _make_row(allele_position='allele1'),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P'),
    ]
    df = pd.DataFrame(rows)
    out = render_report(df, 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    assert 'Incomplete exon coverage</span>' not in html


def test_render_report_missing_pgroup_renders_dash_not_nan(tmp_path):
    # Regression: a NaN p_group must render as — not literal 'nan'
    rows = [
        _make_row(allele_position='allele1', p_group=float('nan'), p_group_found=False),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2',
                  p_group=float('nan'), p_group_found=False),
    ]
    df = pd.DataFrame(rows)
    out = render_report(df, 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    # The P-group cells must not show 'nan'
    assert '>nan<' not in html
    assert 'nan' not in html.replace('lang=', '')  # no stray 'nan' token


def test_render_report_pileup_block_present(tmp_path):
    rows = [
        _make_row(allele_position='allele1'),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P'),
    ]
    df = pd.DataFrame(rows)
    pileup = {
        "A": {
            "chart_png": "iVBORw0KGgo=",  # dummy base64
            "igv_url": "http://localhost:60151/load?file=/d/A.bam",
            "caption": "2 discriminating sites.",
        }
    }
    out = render_report(df, 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR, pileup=pileup)
    html = out.read_text()
    assert 'Read Evidence' in html
    assert 'data:image/png;base64,iVBORw0KGgo=' in html
    assert 'Open in IGV' in html
    assert 'http://localhost:60151/load?file=/d/A.bam' in html


def test_render_report_no_pileup_block_when_absent(tmp_path):
    rows = [_make_row(allele_position='allele1'),
            _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P')]
    df = pd.DataFrame(rows)
    out = render_report(df, 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    assert 'Read Evidence' not in html


def test_render_report_shows_incomplete_bp(tmp_path):
    rows = [
        _make_row(allele_position='allele1', exon3_incomp=7, incomplete_coverage=True, **QC),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P',
                  exon2_incomp=4, incomplete_coverage=True, **QC),
    ]
    html = render_report(pd.DataFrame(rows), 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR).read_text()
    assert '7 bp incomplete' in html
    assert '4 bp incomplete' in html


def test_render_report_flags_incomp_zero(tmp_path):
    # incomp.0: flagged and marked even though no bp count is positive
    rows = [
        _make_row(allele_position='allele1', exon2_incomplete=True, incomplete_coverage=True, **QC),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P'),
    ]
    html = render_report(pd.DataFrame(rows), 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR).read_text()
    assert '<span class="badge badge-red">Incomplete exon coverage</span>' in html
    assert '(!)' in html
    assert '0 bp incomplete' in html


OTHER = dict(has_ambiguous_pair=True, ambiguous_alleles='A*11:01:01 / A*24:02:01',
             phase_check='evaluated')


def test_render_report_other_pair_that_differs_is_never_called_identical(tmp_path):
    html = _html(_pair(other_pair_relation='different', **OTHER), tmp_path)
    assert 'differs from the reported pair' in html
    assert 'identical to the reported pair' not in html


def test_render_report_other_pair_identical_is_stated(tmp_path):
    html = _html(_pair(other_pair_relation='identical', **OTHER), tmp_path)
    assert 'identical to the reported pair' in html
    assert 'differs from the reported pair' not in html


def test_render_report_incomplete_allele_not_evaluated_is_stated(tmp_path):
    html = _html(_pair(incomplete_allele_detail='not evaluated'), tmp_path)
    assert 'Incomplete allele check:</strong> not evaluated' in html


def test_render_report_allele_absent_from_dictionary_is_stated(tmp_path):
    html = _html(_pair(incomplete_allele_detail='A*02:01:01 not in HLA-HD dictionary'), tmp_path)
    assert 'A*02:01:01 not in HLA-HD dictionary' in html


def test_render_report_complete_allele_has_no_incomplete_allele_line(tmp_path):
    assert 'Incomplete allele' not in _html(_pair(), tmp_path)


def test_render_report_allele_cells_drop_hla_prefix(tmp_path):
    html = _html(_pair(), tmp_path)
    assert '<td>A*02:01:01</td>' in html and '<td>HLA-A*02:01:01</td>' not in html


def test_render_report_details_are_not_collapsible(tmp_path):
    html = _html(_pair(flags='Phase ambiguity', phase_check='evaluated'), tmp_path)
    assert '<div class="locus-detail">' in html and '<details' not in html
