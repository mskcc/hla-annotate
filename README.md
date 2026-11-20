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
    [--hlahd_dictionary_dir /path/to/hlahd/dictionary/] \
    [--min_exon_depth 75] \
    [--hlahd_version v1.7.1] \
    [--skip_html]
```

- `--skip_html` writes only the annotated TSV, skipping the HTML report.
- `--imgt_alignment_dir` holds IPD-IMGT/HLA `alignments/{A,B,C}_nuc.txt` from the same
  release as `--pgroup_file`. It is needed for the Phase ambiguity test; without it,
  `phase_check` is `not_evaluated` and a warning is printed.
- `--hlahd_dictionary_dir` is HLA-HD's `dictionary/` directory, the one the run used; needed for
  the Incomplete Allele flag.
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
| QC review | Low exon depth, Incomplete exon coverage, or coverage unknown (`qc_evaluable` false). Combines with the above |
| Failed | HLA-HD reported `Not typed`; no other status applies |

| Flag | When |
|---|---|
| Phase ambiguity | HLA-HD's `#Other ambiguous pair` is a phase conjugate of the reported pair: same unordered bases at every exon 2+3 position, opposite chromosomes (sequence test on the IMGT alignment), **and** HLA-HD could not test it outside exons 2+3: every alternative allele is known only for exons 2+3 in HLA-HD's dictionary (`other_pair_evidence` `unresolved`), or the dictionary is missing or lacks an allele (`not_evaluated`). If any alternative allele has sequence on other exons, HLA-HD scored it there and it lost: not flagged, `other_pair_evidence` `ruled_out`, grey note in the HTML. Untestable over exons 2+3 → treated like a phase conjugate, `phase_check` says so; identical over exons 2+3 → not flagged |
| Shared P-group | Another read-tied best-pair genotype (`#Best allele pair`) maps to the same P-groups |
| Frequency tie | The frequency tiebreak left several genotypes |
| Low exon depth | Exon 2 or 3 mean depth below `--min_exon_depth` for either allele |
| Incomplete exon coverage | HLA-HD reports positions without coverage in exon 2 or 3 for either allele |
| Incomplete Allele | A reported allele lacks one or more exons in HLA-HD's dictionary (`<GENE>/<GENE>.list`, passed with `--hlahd_dictionary_dir`); `incomplete_allele_detail` lists them. Without the dictionary: `not evaluated` |

Shared P-group compares HLA-HD's read-tied best pairs (`#Best allele pair`); Phase
ambiguity comes only from `#Other ambiguous pair` rows (in a 504-sample survey, 36 of
the 37 testable rows were phase conjugates). The two never overlap.
A class I locus missing from
`final.result.txt` is reported as Failed with a warning. `phase_check` records whether the phase test
ran (`evaluated`, `no_other_pair`, `not_evaluated`, `incomplete` when an allele is missing
from the alignment or an `#Other` row could not be read; unreadable rows are flagged, never
dropped). `other_pair_relation` gives each `#Other` row's result (`phase`, `identical`,
`different`, `untestable`): a `different` row is a real alternative genotype, not a phase
swap, and is not flagged, but the HTML says so. `other_pair_evidence` gives each row's
dictionary verdict (`ruled_out`, `unresolved`, `not_evaluated`; same order). In the
504-sample survey (2024 HLA-HD dictionary) 36 of 40 alternatives were `ruled_out` and 4
`unresolved`; HLA-HD's `hla_est` logs confirmed one of each on Iris. The HTML also shows
`incomplete_allele_detail` whenever it is not empty (`not evaluated`, or an allele missing
from the dictionary). `provenance` is `final_only` when `est.txt` is missing.
`multiple_best_pairs` and `has_ambiguous_pair` are kept as raw `est.txt` indicators
for backward compatibility; report `call_status`/`flags` instead.

## Tests

```bash
pytest
```
