from pathlib import Path
import pandas as pd
from hlahd_annotate.annotate import annotate_sample
from hlahd_annotate.report import render_report

TEMPLATE_DIR = Path(__file__).parents[1] / "src" / "hlahd_annotate" / "templates"


def _make_row(
    sample='S1', locus='A', allele='HLA-A*02:01:01', allele_position='allele1',
    resolution=3, p_group='A*02:01P', p_group_found=True,
    multiple_best_pairs=False, has_ambiguous_pair=False, ambiguous_alleles=None,
    est_mismatch=False, est_mismatch_detail=None,
    exon2_depth=262.4, exon2_incomp=0, exon3_depth=303.0, exon3_incomp=0,
    incomplete_coverage=False,
):
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
        'incomplete_coverage': incomplete_coverage,
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


def test_render_report_ambiguous_pair_badge(tmp_path):
    rows = [
        _make_row(allele_position='allele1', has_ambiguous_pair=True,
                  ambiguous_alleles='A*68:01 / A*02:614'),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P',
                  has_ambiguous_pair=True),
    ]
    df = pd.DataFrame(rows)
    out = render_report(df, 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    assert '<span class="badge badge-yellow">AMBIGUOUS PAIR</span>' in html
    assert 'A*68:01' in html  # ambiguous_alleles in flag detail


def test_render_report_multiple_best_pairs_badge(tmp_path):
    rows = [
        _make_row(allele_position='allele1', multiple_best_pairs=True),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P',
                  multiple_best_pairs=True),
        _make_row(allele='HLA-A*68:01:01', allele_position='allele1_pair2',
                  p_group=None, p_group_found=False, multiple_best_pairs=True),
        _make_row(allele='HLA-A*24:02:01', allele_position='allele2_pair2',
                  p_group=None, p_group_found=False, multiple_best_pairs=True),
    ]
    df = pd.DataFrame(rows)
    out = render_report(df, 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    assert '<span class="badge badge-orange">MULTIPLE BEST PAIRS</span>' in html


def test_render_report_est_mismatch_badge(tmp_path):
    rows = [
        _make_row(allele_position='allele1', est_mismatch=True,
                  est_mismatch_detail='final=HLA-A*02:01:01; est_best_alleles=A*03:01,A*24:02'),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P'),
    ]
    df = pd.DataFrame(rows)
    out = render_report(df, 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    assert '<span class="badge badge-red">EST MISMATCH</span>' in html
    assert 'final=HLA-A*02:01:01' in html  # detail in flag section


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
        _make_row(allele_position='allele1', exon3_incomp=7, incomplete_coverage=True),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P',
                  exon2_incomp=4, incomplete_coverage=True),
    ]
    df = pd.DataFrame(rows)
    out = render_report(df, 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    assert '<span class="badge badge-darkred">INCOMPLETE COVERAGE</span>' in html
    assert '(!)' in html  # marker on the incomplete allele's cell


def test_render_report_no_incomplete_badge_when_complete(tmp_path):
    rows = [
        _make_row(allele_position='allele1'),
        _make_row(allele='HLA-A*03:01:01', allele_position='allele2', p_group='A*03:01P'),
    ]
    df = pd.DataFrame(rows)
    out = render_report(df, 'S1', 'pgroup.txt', tmp_path, TEMPLATE_DIR)
    html = out.read_text()
    assert '<span class="badge badge-darkred">INCOMPLETE COVERAGE</span>' not in html


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
