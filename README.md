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
    [--imgt_alignment_dir /path/to/IMGTHLA/alignments/] \
    [--min_exon_depth 75] \
    [--hlahd_version v1.7.1] \
    [--skip_html]
```

- `--skip_html` writes only the annotated TSV, skipping the HTML report.
- `--imgt_alignment_dir` holds IPD-IMGT/HLA `alignments/{A,B,C}_nuc.txt` from the same
  release as `--pgroup_file`. It is needed for the Phase ambiguity test; without it,
  `phase_check` is `not_evaluated` and a warning is printed.
- `--min_exon_depth` is the Low exon depth threshold (HLA-HD mean depth, exon 2 or 3).

## Reporting model

Each locus gets a two-field result (`two_field_genotype`), P-groups, one or more
Call Status values (`call_status`, `;`-joined) and zero or more Flags (`flags`).
These locus-level columns repeat on each allele row of the TSV.

| Call Status | When |
|---|---|
| Unique | One genotype remains after truncating every candidate to two fields |
| Frequency-resolved | Several two-field genotypes tied on reads (`#Best allele pair`); HLA-HD's frequency tiebreak kept one. `frequency_resolution_details` lists the tied genotypes |
| Multiple best pairs | Several two-field genotypes remain in `final.result.txt` (always with the Frequency tie flag) |
| QC review | Low exon depth or Incomplete exon coverage, or coverage unknown (`qc_evaluable` false). Combines with the above |
| Failed | HLA-HD reported `Not typed`; no other status applies |

| Flag | When |
|---|---|
| Phase ambiguity | A competing genotype has the same unordered bases at every exon 2+3 position but on the opposite chromosomes (sequence test on the IMGT alignment) |
| Shared P-group | A competing two-field genotype maps to the same P-groups |
| Frequency tie | The frequency tiebreak left several genotypes |
| Low exon depth | Exon 2 or 3 mean depth below `--min_exon_depth` for either allele |
| Incomplete exon coverage | HLA-HD reports positions without coverage in exon 2 or 3 for either allele |

Competing genotypes are the read-tied candidates, the reported pairs and every
`#Other ambiguous pair` row of `est.txt`. `phase_check` records whether the phase test
ran (`evaluated`, `no_competing`, `not_evaluated`, `incomplete` when an allele is missing
from the alignment). `provenance` is `final_only` when `est.txt` is missing.
`multiple_best_pairs` and `has_ambiguous_pair` are kept as raw `est.txt` indicators
for backward compatibility; report `call_status`/`flags` instead.

## Tests

```bash
pytest
```
