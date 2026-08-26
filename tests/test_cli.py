import subprocess
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"

def test_cli_produces_tsv(pgroup_file, tmp_path):
    result = subprocess.run(
        ['annotate_hlahd',
         '--result_dir', str(FIXTURES),
         '--sample', 'sample3',
         '--pgroup_file', str(pgroup_file),
         '--outdir', str(tmp_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / 'sample3_annotated.tsv').exists()

def test_cli_produces_html(pgroup_file, tmp_path):
    result = subprocess.run(
        ['annotate_hlahd',
         '--result_dir', str(FIXTURES),
         '--sample', 'sample3',
         '--pgroup_file', str(pgroup_file),
         '--outdir', str(tmp_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    html = (tmp_path / 'sample3_report.html').read_text()
    assert 'sample3' in html
    assert 'HLA-A' in html

def test_cli_missing_required_arg():
    result = subprocess.run(
        ['annotate_hlahd', '--result_dir', '/tmp'],
        capture_output=True, text=True,
    )
    assert result.returncode != 0

def test_cli_outdir_created_if_missing(pgroup_file, tmp_path):
    new_outdir = tmp_path / 'new_subdir'
    result = subprocess.run(
        ['annotate_hlahd',
         '--result_dir', str(FIXTURES),
         '--sample', 'sample3',
         '--pgroup_file', str(pgroup_file),
         '--outdir', str(new_outdir)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert new_outdir.exists()


def test_cli_custom_hlahd_version(pgroup_file, tmp_path):
    result = subprocess.run(
        ['annotate_hlahd',
         '--result_dir', str(FIXTURES),
         '--sample', 'sample3',
         '--pgroup_file', str(pgroup_file),
         '--outdir', str(tmp_path),
         '--hlahd_version', 'v2.0.0'],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    html = (tmp_path / 'sample3_report.html').read_text()
    assert 'v2.0.0' in html


def test_cli_skip_html(pgroup_file, tmp_path):
    result = subprocess.run(
        ['annotate_hlahd',
         '--result_dir', str(FIXTURES),
         '--sample', 'sample3',
         '--pgroup_file', str(pgroup_file),
         '--outdir', str(tmp_path),
         '--skip_html'],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / 'sample3_annotated.tsv').exists()
    assert not (tmp_path / 'sample3_report.html').exists()
