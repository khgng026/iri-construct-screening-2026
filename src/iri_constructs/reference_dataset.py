"""Prepare source-free, unscored reference evidence; never grant T02/T07 approval.

The catalogue is a prospective proposal. Family guards check externally curated
evidence, while exact-copy/parent/lineage exclusions remain fail-closed. No model
or alignment scores are calculated here.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

from .runtime import InputError

AA = frozenset('ACDEFGHIKLMNPQRSTVWY')
FAMILY_RULES = {
    'gfp_like_homologs': {
        'pfam_all': ['PF01353'],
        'similarity_contains': 'Belongs to the GFP family',
        'scope': 'GFP molecular homolog family, including fluorescent and nonfluorescent chromoproteins; color is not a class.',
    },
    'coleopteran_luciferase_homologs': {
        'pfam_all': ['PF00501', 'PF13193'],
        'similarity_contains': 'ATP-dependent AMP-binding enzyme family',
        'lineage_contains': 'Coleoptera',
        'recommended_name': 'Luciferin 4-monooxygenase',
        'scope': 'Curated coleopteran luciferase homolog proposal, narrower than the AMP-binding superfamily; domain matches alone do not define luciferases.',
    },
    'renilla_luciferase_lineage': {
        'pfam_all': ['PF00561'],
        'interpro_all': ['IPR000639'],
        'organism': 'Renilla reniformis',
        'recommended_name': 'Coelenterazine h 2-monooxygenase',
        'scope': 'Renilla lineage anchor only; not an independent positive for Renilla-derived query constructs. Not the firefly homolog family.',
    },
    'omp_decarboxylase_homologs': {
        'pfam_all': ['PF00215'],
        'similarity_contains': 'Belongs to the OMP decarboxylase family',
        'scope': 'Standalone OMP-decarboxylase molecular homolog family, used here for a source-documented URA3 component; no inference from growth selection alone.',
    },
}


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _require(condition, message):
    if not condition:
        raise InputError(message)


def _references(document):
    result = []
    for ref in document.get('references', []):
        citation = ref.get('citation', {})
        xrefs = {r['database']: r['id'] for r in citation.get('citationCrossReferences', [])}
        result.append({
            'title': citation.get('title'), 'publication_date': citation.get('publicationDate'),
            'pubmed_id': xrefs.get('PubMed'), 'doi': xrefs.get('DOI'),
            'citation_type': citation.get('citationType'),
            'reference_positions': ref.get('referencePositions', []),
        })
    return result


def summarize_reference(document, specification, *, raw_sha256, retrieved_at_utc):
    """Validate exact reviewed record and catalogue-specific family evidence.

    Returns no sequence. This structural/evidence check is not independent
    scientific clearance, even for externally reviewed Swiss-Prot records.
    """
    accession = specification.get('accession')
    _require(bool(re.fullmatch(r'[A-Z0-9]{6,10}', str(accession))), 'Invalid UniProt accession specification')
    _require(document.get('primaryAccession') == accession, 'Reference accession changed or redirected')
    _require(document.get('entryType') == 'UniProtKB reviewed (Swiss-Prot)',
             'Reference is not an exact reviewed Swiss-Prot entry')
    sequence = document.get('sequence', {}).get('value', '')
    _require(isinstance(sequence, str) and sequence and not set(sequence) - AA,
             'Reference has an empty or unsupported protein alphabet')
    _require(document['sequence'].get('length') == len(sequence), 'Reference sequence length is inconsistent')
    _require(len(sequence) <= 1022, 'Reference exceeds the frozen shared 1022-residue cap')
    audit = document.get('entryAudit', {})
    _require(type(audit.get('entryVersion')) is int and type(audit.get('sequenceVersion')) is int
             and audit['entryVersion'] > 0 and audit['sequenceVersion'] > 0,
             'Reference lacks immutable entry and sequence versions')
    organism = document.get('organism', {})
    _require(organism.get('taxonId') == specification.get('taxon_id'), 'Reference organism differs from proposal')
    family = specification.get('proposed_family')
    _require(family in FAMILY_RULES, 'Unsupported reference family proposal')
    rule = FAMILY_RULES[family]
    crossrefs = document.get('uniProtKBCrossReferences', [])
    pfam = sorted({r['id'] for r in crossrefs if r.get('database') == 'Pfam'})
    interpro = sorted({r['id'] for r in crossrefs if r.get('database') == 'InterPro'})
    similarity = [t for c in document.get('comments', []) if c.get('commentType') == 'SIMILARITY'
                  for t in c.get('texts', [])]
    _require(set(rule.get('pfam_all', [])) <= set(pfam), 'Required curated Pfam evidence absent')
    _require(set(rule.get('interpro_all', [])) <= set(interpro), 'Required curated InterPro evidence absent')
    if rule.get('similarity_contains'):
        _require(any(rule['similarity_contains'] in t.get('value', '') for t in similarity),
                 'Required reviewed family annotation absent')
    if rule.get('lineage_contains'):
        _require(rule['lineage_contains'] in organism.get('lineage', []), 'Curated reference taxonomic lineage mismatch')
    name = document.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value')
    if rule.get('recommended_name'):
        _require(name == rule['recommended_name'], 'Curated reference protein identity mismatch')
    if rule.get('organism'):
        _require(organism.get('scientificName') == rule['organism'], 'Curated reference source organism mismatch')
    citations = _references(document)
    _require(any(c['pubmed_id'] or c['doi'] for c in citations), 'No independently citable source in reference entry')
    _require(bool(re.fullmatch(r'[0-9a-f]{64}', raw_sha256)), 'Missing original reference response fingerprint')
    _require(bool(retrieved_at_utc), 'Missing reference retrieval date')
    return {
        'reference_id': 'uniprot-' + accession, 'accession': accession,
        'entry_id': document.get('uniProtkbId'), 'entry_type': document['entryType'],
        'entry_version': audit['entryVersion'], 'sequence_version': audit['sequenceVersion'],
        'last_sequence_update': audit.get('lastSequenceUpdateDate'),
        'last_annotation_update': audit.get('lastAnnotationUpdateDate'),
        'source_url': f'https://rest.uniprot.org/uniprotkb/{accession}.json',
        'display_url': f'https://www.uniprot.org/uniprotkb/{accession}/entry',
        'retrieved_at_utc': retrieved_at_utc, 'original_json_sha256': raw_sha256,
        'protein_sequence_sha256': hashlib.sha256(sequence.encode('ascii')).hexdigest(),
        'length_aa': len(sequence), 'recommended_name': name,
        'organism': organism.get('scientificName'), 'taxon_id': organism['taxonId'],
        'proposed_family': family, 'family_scope': rule['scope'],
        'family_evidence': {'pfam': pfam, 'interpro': interpro,
                            'reviewed_similarity_annotations': similarity, 'citations': citations},
        'lineage_group': specification['lineage_group'],
        'ancestor_anchor': specification.get('ancestor_anchor', False),
        'independence_caveat': specification['independence_caveat'],
        'benign_basis': 'Naturally occurring reporter/chromoprotein or ordinary pyrimidine-biosynthesis homolog used for protein-family retrieval evidence; no sequence generation, hazardous source panel or function engineering.',
        'candidate_selected': True, 'source_rights_proposal': 'CC_BY_4_0_with_UniProt_attribution',
        'local_analysis_proposed_permitted': True,
        'private_cloud_copy_proposed_permitted': True,
        'source_admitted_for_scoring': False, 'family_panel_approved': False,
        'run_authorized': False, 'review_status': 'READY_FOR_DELEGATED_REVIEW',
    }


def proposed_eligibility(query, references):
    """Exclude copies of ANY documented parent component before any scoring.

    Unknown parent components still enter the copy-exclusion set. Provenance and
    family decisions are proposals only; missing lineage/family evidence blocks
    eligibility rather than being interpreted as independence.
    """
    required = ('query_id', 'parent_id', 'protein_sequence_sha256', 'proposed_family',
                'lineage_group', 'family_evidence')
    _require(all(query.get(k) for k in required), 'Query lacks family or lineage evidence')
    _require(query['proposed_family'] in FAMILY_RULES, 'Query family is outside the proposed ontology')
    components = query.get('all_parent_component_sha256')
    _require(isinstance(components, list) and components
             and query['protein_sequence_sha256'] in components,
             'Query must bind all documented parent component fingerprints including itself')
    component_hashes = set(components)
    _require(all(re.fullmatch(r'[0-9a-f]{64}', h) for h in component_hashes), 'Invalid parent component fingerprint')
    known_ancestors = set(query.get('known_ancestor_accessions', []))
    allowed, removals = [], {}
    seen = set()
    for reference in sorted(references, key=lambda r: r['reference_id']):
        rid = reference['reference_id']
        _require(rid not in seen, 'Duplicate reference identity')
        seen.add(rid)
        _require(reference.get('lineage_group') and reference.get('proposed_family') in FAMILY_RULES,
                 'Reference lacks proposed family or source lineage')
        reasons = []
        if reference['protein_sequence_sha256'] in component_hashes:
            reasons.append('exact_copy_of_any_documented_parent_component')
        if reference.get('parent_id') == query['parent_id']:
            reasons.append('same_source_parent')
        if reference['lineage_group'] == query['lineage_group']:
            reasons.append('shared_origin_or_vector_lineage')
        if reference['accession'] in known_ancestors:
            reasons.append('documented_query_ancestor')
        if reasons:
            removals[rid] = reasons
        else:
            allowed.append(reference)
    positives = [r['reference_id'] for r in allowed if r['proposed_family'] == query['proposed_family']]
    competitors = [r['reference_id'] for r in allowed if r['proposed_family'] != query['proposed_family']]
    return {
        'query_id': query['query_id'], 'parent_id': query['parent_id'],
        'proposed_family': query['proposed_family'],
        'eligible_reference_ids_proposed': [r['reference_id'] for r in allowed],
        'positive_reference_ids_proposed': positives,
        'competing_reference_ids_proposed': competitors, 'removed_references': removals,
        'eligible_count': len(allowed), 'positive_count': len(positives), 'competing_count': len(competitors),
        'top_k': 5, 'budget_saturated': len(allowed) <= 5,
        'support_status': ('PROPOSED_SUPPORTED_PENDING_REVIEW' if positives and competitors
                           else 'HOLD_NO_INDEPENDENT_POSITIVE_OR_COMPETITOR'),
        'scientific_clearance': False, 'scoring_performed': False,
    }


def freeze_proposal(references, queries, *, source_manifest_sha256, catalogue_spec_sha256):
    """Bind evidence/panel identities prospectively without making it approved."""
    refs = sorted(references, key=lambda r: r['reference_id'])
    _require(len({r['reference_id'] for r in refs}) == len(refs), 'Duplicate reference IDs')
    _require(len({r['protein_sequence_sha256'] for r in refs}) == len(refs),
             'Exact duplicate reference proteins must be resolved before panel proposal')
    qrows = sorted(queries, key=lambda q: q['query_id'])
    _require(len({q['query_id'] for q in qrows}) == len(qrows), 'Duplicate query IDs')
    eligibility = [proposed_eligibility(q, refs) for q in qrows]
    result = {
        'schema_version': 'week2-reference-proposal-1',
        'review_status': ('READY_FOR_DELEGATED_REVIEW' if refs and qrows
                          else 'CATALOGUE_ONLY_PENDING_QUERY_EVIDENCE' if refs
                          else 'HOLD_NO_REFERENCE_EVIDENCE'),
        'scientific_clearance': False, 'inference_executed': False, 'scoring_performed': False,
        'source_manifest_sha256': source_manifest_sha256,
        'catalogue_spec_sha256': catalogue_spec_sha256,
        'reference_ids': [r['reference_id'] for r in refs],
        'references': refs, 'query_proposals': qrows, 'eligibility_proposals': eligibility,
        'top_k': 5, 'model_residue_cap': 1022,
        'strict_relatedness_policy': 'HOLD_missing_independently_reviewed_fingerprint_bound_policy',
        'limitations': [
            'This is an evidence proposal, not a scoring panel accepted by the experimental runners.',
            'Species/source differences do not prove evolutionary, pretraining or statistical independence.',
            'Families denote proposed homology groups, not equivalent biochemical activity or fluorescence color.',
            'Natural ancestor proteins are excluded from their engineered derivative query lineages.',
            'T02 rights/source admission and ChatGPT delegated scientific clearance remain separate gates.',
        ],
    }
    result['proposal_sha256'] = canonical_hash(result)
    return result


def build_source_queries(source_review, unit_provenance, specification):
    """Bind explicit evidence proposals to exact source CDSs and normalized units.

    No family is guessed from a product string. The human-readable specification
    must name the source protein, proposed homolog family and supporting facts;
    these remain proposals requiring separate delegated scientific review.
    """
    sources = {r['source_id']: r for r in source_review['records']}
    units = unit_provenance['units']
    queries, parent_coverage = [], {}
    seen = set()
    for spec in specification['queries']:
        source_id, protein_id = spec['source_id'], spec['protein_accession_version']
        _require((source_id, protein_id) not in seen, 'Duplicate source protein query specification')
        seen.add((source_id, protein_id))
        source = sources.get(source_id)
        _require(source is not None and source.get('eligible_for_preparatory_QA') is True
                 and source.get('status') == 'READY_FOR_REVIEW',
                 'Query source is absent or held from preparatory QA')
        matches = [c for c in source.get('cds_evidence', [])
                   if c.get('protein_accession_version') == protein_id]
        _require(len(matches) == 1 and matches[0].get('translation_state') == 'known',
                 'Query must name one validated source-provided CDS translation')
        cds = matches[0]
        candidates = [u for u in units if u['source_id'] == source_id
                      and u['sequence_sha256'] == cds['translation_sha256']
                      and u['source_segments'] == cds['segments']]
        _require(len(candidates) == 1, 'Query source CDS must map to one deterministic normalized unit')
        unit = candidates[0]
        _require(unit['source_sha256'] == source['source_retrieval']['sha256']
                 and unit['length_aa'] == cds['translation_length_aa'],
                 'Query unit differs from the original source evidence')
        _require(unit['length_aa'] <= 1022, 'Query exceeds the frozen model residue cap')
        _require(spec.get('proposed_family') in FAMILY_RULES and spec.get('family_evidence')
                 and spec.get('lineage_group') and spec.get('known_ancestor_accessions'),
                 'Query specification lacks explicit homolog-family and ancestry evidence')
        components = source.get('cds_evidence', [])
        component_hashes = sorted({c['translation_sha256'] for c in components
                                   if c.get('translation_sha256')})
        queries.append({
            **spec, 'query_id': unit['input_id'], 'parent_id': source_id,
            'source_url': source['source_url'],
            'source_sha256': source['source_retrieval']['sha256'],
            'protein_sequence_sha256': unit['sequence_sha256'], 'length_aa': unit['length_aa'],
            'feature_ids': unit['feature_ids'], 'source_segments': unit['source_segments'],
            'all_parent_component_sha256': component_hashes,
            'review_status': 'PROPOSED_NOT_APPROVED', 'source_admitted_for_scoring': False,
            'family_panel_approved': False, 'run_authorized': False,
        })
        supported_ids = {s['protein_accession_version'] for s in specification['queries']
                         if s['source_id'] == source_id}
        eligible = [c for c in components if c.get('protein_accession_version') in supported_ids
                    and c.get('translation_state') == 'known']
        parent_coverage[source_id] = {
            'source_documented_cds_count': len(components),
            'validated_source_translation_count': sum(c.get('translation_state') == 'known' for c in components),
            'proposed_supported_component_count': len(eligible),
            'proposed_supported_distinct_families': sorted({s['proposed_family'] for s in specification['queries']
                                                          if s['source_id'] == source_id}),
            'full_documented_family_set_status': 'HOLD_unreviewed_component_families',
            'unsupported_components': [
                {'protein_accession_version': c.get('protein_accession_version'), 'product': c.get('product'),
                 'length_aa': c.get('translation_length_aa'), 'translation_sha256': c.get('translation_sha256'),
                 'reason': 'No supported family/reference proposal; preserved outside the supported-subset denominator'}
                for c in components if c.get('protein_accession_version') not in supported_ids],
            'source_longest_translation_aa': max((c.get('translation_length_aa') or 0 for c in components), default=None),
            'proposed_eligible_longest_translation_aa': max((c['translation_length_aa'] for c in eligible), default=None),
            'scientific_clearance': False,
        }
    return {'schema_version': 'week2-query-proposals-1', 'queries': queries,
            'parent_component_coverage': dict(sorted(parent_coverage.items())),
            'scientific_clearance': False, 'scoring_performed': False}
