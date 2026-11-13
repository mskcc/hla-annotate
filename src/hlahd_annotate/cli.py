#!/usr/bin/env python3
"""
annotate_hlahd — Post-process HLA-HD class I output.

Produces per-sample annotated TSV and (optionally) a self-contained HTML report.

Usage:
    annotate_hlahd \
        --result_dir results/SAMPLE1/ \
        --sample SAMPLE1 \
        --pgroup_file /path/to/hla_nom_p.txt \
        --outdir output/
"""
import argparse
import sys
from importlib.metadata import version
from pathlib import Path

import pandas as pd

from hlahd_annotate.annotate import DEFAULT_MIN_EXON_DEPTH, annotate_sample
from hlahd_annotate.report import render_report

TEMPLATE_DIR = Path(__file__).parent / "templates"


def parse_args():
    p = argparse.ArgumentParser(
        description='Annotate HLA-HD class I output with P-groups and quality flags.'
    )
    p.add_argument('--version', action='version',
                   version=f"%(prog)s {version('hla-annotate')}")
    p.add_argument('--result_dir', required=True, type=Path,
                   help='Directory containing <sample>_final.result.txt and <sample>_{A,B,C}.est.txt')
    p.add_argument('--sample', required=True,
                   help='Sample ID (used as filename prefix for output files)')
    p.add_argument('--pgroup_file', required=True, type=Path,
                   help='IMGT wmda/hla_nom_p.txt P-group reference table')
    p.add_argument('--outdir', required=True, type=Path,
                   help='Output directory for TSV and HTML report (created if absent)')
    p.add_argument('--hlahd_version', default='v1.7.1',
                   help='HLA-HD version string for the report footer (default: v1.7.1)')
    p.add_argument('--skip_html', action='store_true',
                   help='Skip generating the HTML report; write only the annotated TSV')
    p.add_argument('--imgt_alignment_dir', type=Path, default=None,
                   help='Directory with IPD-IMGT/HLA alignments/{A,B,C}_nuc.txt (same release as '
                        '--pgroup_file); enables the Phase ambiguity test')
    p.add_argument('--min_exon_depth', type=float, default=DEFAULT_MIN_EXON_DEPTH,
                   help=f'Low exon depth threshold: exon 2/3 mean depth (default: {DEFAULT_MIN_EXON_DEPTH:g})')
    return p.parse_args()


def _summarize(df: pd.DataFrame) -> None:
    """Counts line per sample; names loci whose phase test could not run."""
    loci = df.drop_duplicates('locus')
    counts = loci['call_status'].str.split(';').explode().value_counts()
    print("CALL STATUS: " + ', '.join(f"{k}={v}" for k, v in counts.items()), file=sys.stderr)
    skipped = loci[loci['phase_check'].isin(['not_evaluated', 'incomplete'])]
    for _, r in skipped.iterrows():
        print(f"WARNING: Phase ambiguity not fully evaluated for locus {r['locus']} "
              f"({r['phase_check']})", file=sys.stderr)


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = annotate_sample(args.result_dir, args.sample, args.pgroup_file,
                         alignment_dir=args.imgt_alignment_dir,
                         min_exon_depth=args.min_exon_depth)
    if not df.empty:
        _summarize(df)

    incomplete = df[df['incomplete_coverage'].astype(bool)]
    if not incomplete.empty:
        known = incomplete['incomplete_bp'].dropna()
        total = f"{int(known.sum())} bp" + (" (some counts unknown)" if len(known) < len(incomplete) else "")
        print(f"INCOMPLETE COVERAGE: {len(incomplete)} allele(s), {total}", file=sys.stderr)
        for _, r in incomplete.iterrows():
            bp = f"{int(r['incomplete_bp'])} bp" if pd.notna(r['incomplete_bp']) else "unknown bp"
            print(f"  {r['locus']} {r['allele_position']} {r['allele']}: {bp}", file=sys.stderr)

    tsv_path = args.outdir / f"{args.sample}_annotated.tsv"
    df.to_csv(tsv_path, sep='\t', index=False)
    print(f"TSV written: {tsv_path}")

    if not args.skip_html:
        html_path = render_report(df, args.sample, str(args.pgroup_file), args.outdir, TEMPLATE_DIR,
                                  hlahd_version=args.hlahd_version,
                                  min_exon_depth=args.min_exon_depth)
        print(f"HTML written: {html_path}")


if __name__ == '__main__':
    main()
