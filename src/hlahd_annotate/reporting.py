"""Reporting model: per-locus Call Status, Flags and Frequency Resolution Details.

Call Status (one or more):
  Unique               one genotype remains after normalizing to two fields
  Frequency-resolved   several two-field genotypes tied on reads (est.txt
                       '#Best allele pair'); the frequency tiebreak kept one
  Multiple best pairs  several two-field genotypes remain in final.result.txt
  QC review            Low exon depth or Incomplete exon coverage (or coverage
                       unknown); combines with any of the above
  Failed               HLA-HD reported 'Not typed'; nothing else applies
Flags:
  Phase ambiguity      a competing genotype carries the same unordered bases at
                       every exon 2+3 position but on the other chromosome
  Shared P-group       a competing two-field genotype has the same P-groups
  Frequency tie        the frequency tiebreak left several genotypes
  Low exon depth       exon 2 or 3 mean depth below threshold, either allele
  Incomplete exon coverage
Competing genotypes are every other two-field genotype among HLA-HD's read-tied
best-pair candidates ('#Best allele pair'). '#Other ambiguous pair' rows are NOT
competing: HLA-HD ranks them below the best pair and prints no coverage for them
(504-sample survey 2026-11-13: always one row, always a different two-field
genotype). They stay in the TSV (ambiguous_alleles) as information only.
"""
from .alignment import phase_relation
from .pgroup import lookup_pgroup

def two_field(allele: str) -> str:
    """Truncate to two fields, as HLA-HD's pick_up_allele truncates to three:
    'HLA-A*02:01:01:02L' -> 'A*02:01'. A suffix that belongs to the second
    field stays ('A*24:09N'); one on a later field is dropped with it, so
    expression variants of a 3/4-field group never count as separate genotypes."""
    return ':'.join(allele.removeprefix('HLA-').split(':')[:2])


def genotype(pair: tuple[str, str]) -> tuple[str, str]:
    return tuple(sorted(two_field(a) for a in pair))


def fmt_genotype(g: tuple[str, str]) -> str:
    return f"{g[0]}/{g[1]}"


def candidate_genotypes(slots: list[tuple[list, list]]) -> dict:
    """Read-tied two-field genotypes from est.txt best pairs -> representative
    full-name pair (first seen). An empty second slot ('-') means homozygous:
    each slot-1 allele pairs only with itself."""
    out: dict = {}
    for s1, s2 in slots:
        pairs = [(a, b) for a in s1 for b in s2] if s2 else [(a, a) for a in s1]
        for a, b in pairs:
            out.setdefault(genotype((a, b)), (a, b))
    return out


def classify(final_pairs, est, pgroups, typed_exons, coverage, min_depth):
    """Reporting fields for one locus.

    final_pairs: [(a1, a2)] from final.result.txt, None for 'Not typed'
    est: parse_est() result; pgroups: P-group lookup
    typed_exons: {allele: exon 2+3 columns} or None when no alignment was given
    coverage: coverage dicts of every reported allele
    """
    pairs = [p for p in final_pairs if p[0] is not None and p[1] is not None]
    if not pairs:
        return _failed()

    reported = {}
    for p in pairs:
        reported.setdefault(genotype(p), p)
    est_cands = candidate_genotypes(est.get('best_pair_slots', []))
    provenance = 'est' if est_cands else 'final_only'
    cands = est_cands or dict(reported)
    for g in reported:  # prefer est.txt names: they are what HLA-HD scored
        reported[g] = cands.get(g, reported[g])

    status, flags, details = _resolution(reported, cands)
    primary = next(iter(reported))
    competing = {g: r for g, r in cands.items() if g != primary}
    if _shares_pgroup(primary, competing, pgroups):
        flags.append('Shared P-group')
    phase_check, phased = _phase(reported[primary], competing, typed_exons)
    if phased:
        flags.insert(0, 'Phase ambiguity')

    qc_flags, qc_evaluable = _qc(coverage, min_depth)
    flags += qc_flags
    if qc_flags or not qc_evaluable:
        status.append('QC review')

    return {
        'two_field_genotype': '; '.join(fmt_genotype(g) for g in reported),
        'call_status': status, 'flags': flags,
        'frequency_resolution_details': details,
        'provenance': provenance, 'phase_check': phase_check,
        'qc_evaluable': qc_evaluable,
    }


def _failed():
    return {'two_field_genotype': None, 'call_status': ['Failed'], 'flags': [],
            'frequency_resolution_details': None, 'provenance': 'final_only',
            'phase_check': 'not_applicable', 'qc_evaluable': False}


def _resolution(reported: dict, cands: dict):
    """Call status from reported vs read-tied two-field genotypes."""
    details = None
    if len(cands) > 1 and len(cands) > len(reported):
        details = '; '.join(fmt_genotype(g) for g in cands)
    if len(reported) > 1:
        return ['Multiple best pairs'], ['Frequency tie'], details
    if len(cands) > 1:
        return ['Frequency-resolved'], [], details
    return ['Unique'], [], None


def _pgroup_pair(g, pgroups):
    found = [lookup_pgroup(a, pgroups) for a in g]
    if not all(ok for _, ok in found):
        return None
    return tuple(sorted(p for p, _ in found))


def _shares_pgroup(primary, competing: dict, pgroups) -> bool:
    mine = _pgroup_pair(primary, pgroups)
    return mine is not None and any(_pgroup_pair(g, pgroups) == mine for g in competing)


def _phase(primary_rep, competing: dict, typed_exons):
    """('no_competing'|'not_evaluated'|'incomplete'|'evaluated', any phase conjugate)."""
    if not competing:
        return 'no_competing', False
    if typed_exons is None:
        return 'not_evaluated', False
    complete, phased = True, False
    for rep in competing.values():
        alleles = (*primary_rep, *rep)
        if not all(a in typed_exons for a in alleles):
            complete = False
            continue
        x = (typed_exons[alleles[0]], typed_exons[alleles[1]])
        y = (typed_exons[alleles[2]], typed_exons[alleles[3]])
        relation, unknown = phase_relation(x, y)
        if unknown >= len(x[0]):  # nothing compared: not evidence of anything
            complete = False
        phased |= relation == 'phase'
    return ('evaluated' if complete else 'incomplete'), phased


def _qc(coverage: list[dict], min_depth: float):
    """(QC flags, evaluable). Unknown depth is not a pass: evaluable=False."""
    low = incomplete = False
    evaluable = bool(coverage)
    for cov in coverage:
        for exon in ('exon2', 'exon3'):
            depth = cov.get(f'{exon}_depth')
            if depth is None:
                evaluable = False
            elif depth < min_depth:
                low = True
        incomplete |= bool(cov.get('exon2_incomplete') or cov.get('exon3_incomplete'))
    flags = (['Low exon depth'] if low else []) + (['Incomplete exon coverage'] if incomplete else [])
    return flags, evaluable
