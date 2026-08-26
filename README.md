# hla-annotate

Post-processes [HLA-HD](https://w3.genome.med.kyoto-u.ac.jp/HLA-HD/) class I output
(`<sample>_final.result.txt` + per-locus `<sample>_{A,B,C}.est.txt`) into a
per-allele TSV annotated with IMGT P-groups and quality flags, plus an optional
self-contained HTML report.

Extracted from [mskcc/HLA_HD_workflow](https://github.com/mskcc/HLA_HD_workflow) so it can be
shared between that pipeline and the [`annotate_hlahd`](https://github.com/mskcc-omics-workflows/modules)
Nextflow module without vendoring copies of the code.

## Install

```bash
pip install .
# or, for development:
pip install -e '.[dev]'
```

## Usage

```bash
annotate_hlahd \
    --result_dir results/SAMPLE1/ \
    --sample SAMPLE1 \
    --pgroup_file /path/to/hla_nom_p.txt \
    --outdir output/ \
    [--hlahd_version v1.7.1] \
    [--skip_html]
```

- `--skip_html` writes only the annotated TSV, skipping the HTML report.

## Tests

```bash
pytest
```
