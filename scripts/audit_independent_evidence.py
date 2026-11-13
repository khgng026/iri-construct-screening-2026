"""Reproduce bounded source/reference audit from privately cached provider originals.

Acquisition URLs and hashes are retained in the generated source-free proposal;
this command performs no downloads, scoring, translation invention or admission.
"""
import hashlib
import json
import sys
from pathlib import Path
from urllib.parse import urlencode
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.runtime import InputError, digest_file, write_json


def require(condition, reason):
    if not condition:
        raise InputError(reason)


def main():
    cache = ROOT / 'results/independence-2026-11-12'
    output = ROOT / 'docs/results/analysis-2026-11-12'
    source = json.loads((ROOT / 'data/manifests/week2/source_review.json').read_text())
    original = json.loads((ROOT / 'data/manifests/week2/references/reference_proposal.json').read_text())
    groups = {}
    candidates = []
    for r in source['records']:
        require(digest_file(ROOT / r['original_path']) == r['source_retrieval']['sha256'], 'Changed acquired source: ' + r['source_id'])
        for cds in r.get('cds_evidence', []):
            if cds.get('translation_state') == 'known':
                groups.setdefault(cds['translation_sha256'], []).append({'source_id': r['source_id'], 'product': cds['product'], 'protein_id': cds.get('protein_accession_version'), 'covered_by_current_family_proposal': cds.get('protein_accession_version') in ['AAA89082.1','AAA89084.1','AAV52869.1','AAW66991.1','AAB02574.1','AAG34528.1','AAG34529.1']})
        candidates.append({'source_id':r['source_id'], 'provided_validated_proteins':r.get('validated_source_translation_count',0), 'archetype':r.get('proposed_archetype'), 'preparatory_eligible':r.get('eligible_for_preparatory_QA',False), 'scored_expansion':False})
    renilla = cache / 'renilla-reviewed.json'
    accessions = [r['primaryAccession'] for r in json.loads(renilla.read_text())['results']]
    require(accessions == ['P27652'], 'Changed Renilla evidence: review search outcome')
    audit = {'status':'PROSPECTIVE_EXPANSION_REQUIRES_SEPARATE_REVIEW', 'existing_candidates':candidates,
             'existing_unique_translation_groups':list(groups.values()),
             'renilla_reviewed_search':{'url':'https://rest.uniprot.org/uniprotkb/search?' + urlencode({'query':'(protein_name:luciferase) AND (taxonomy_name:Renilla) AND (reviewed:true)','format':'json','size':50}), 'original_response_sha256':digest_file(renilla),'accessions':accessions,'conclusion':'Only the already excluded ancestor was found. No added independent Renilla positive; exclusion retained. Search scope is reviewed UniProt, not all possible literature.'},
             'limitations':['Ten iGEM candidates lack provided translations; native exports inspected previously did not repair this.', 'Four original GenBank candidates remain held for annotation limitations.', 'Three firefly constructs encode one identical reporter; marker/backbone sharing also connects parents.', 'Distinct source accessions or translations alone do not establish independent clusters.'], 'new_scored_panel':False}
    new_sources = []
    for accession, path in [('AF234299.1', cache/'pCAMBIA1303.gb'), ('MK448012.1', ROOT/'results/week2/supplemental-fusion/raw_sources/MK448012.1.gb')]:
        expected = {'AF234299.1':'2a56e2ae5eaaffcddb1e7bd7266047fe60b778253ede1fcacae6251ccda1c480', 'MK448012.1':'038fe3d5fb47562aca2003bc20a06a91809ccab8380221d4573976214b952c07'}
        require(digest_file(path) == expected[accession], 'Changed prospective source bytes')
        record = SeqIO.read(path, 'genbank')
        require(record.id == accession, 'Wrong supplemental accession')
        evidence = []
        for i, feature in enumerate(record.features):
            if feature.type != 'CDS':
                continue
            provided = ''.join(feature.qualifiers.get('translation',[''])[0].split()).rstrip('*')
            translated = str(feature.extract(record.seq)[int(feature.qualifiers.get('codon_start',['1'])[0])-1:].translate(table=int(feature.qualifiers.get('transl_table',['1'])[0]), to_stop=True))
            require(provided and translated == provided, 'Supplemental source translation mismatch')
            evidence.append({'source_feature_ordinal_1based':i+1, 'location':str(feature.location), 'protein_id':feature.qualifiers.get('protein_id',[None])[0], 'product':feature.qualifiers.get('product',[None])[0], 'length_aa':len(provided), 'sequence_sha256':hashlib.sha256(provided.encode()).hexdigest()})
        new_sources.append({'accession':accession, 'title':record.description, 'archetype':'tagged_fusion', 'private_cache_path':str(path.relative_to(ROOT)), 'original_sha256':digest_file(path), 'sequence_sha256':hashlib.sha256(str(record.seq).encode()).hexdigest(), 'validated_cds':evidence, 'proposal':'E1 annotation preservation only; no invented fusion component boundaries or new family truth', 'rights':'NCBI molecular-data policy; no conflicting notice found; submitter IP caveat; private analysis only, no public source release', 'rights_policy':'https://www.ncbi.nlm.nih.gov/home/about/policies/', 'benign_context':'Annotated research reporter vector; source regulatory and resistance-marker annotations retained, no engineering or functional validation', 'relatedness_limit':'Shared reporter ancestry/backbone/markers preclude treating accession count as independent clusters', 'source_url':'https://www.ncbi.nlm.nih.gov/nuccore/'+accession})
    ref_path = cache/'O58462.json'
    require(digest_file(ref_path) == 'a31843af945087f53c4717fbe4cba7a93e4151c12de51534b557c5a4061e60a4', 'Changed prospective reference bytes')
    ref = json.loads(ref_path.read_text())
    require(ref['primaryAccession'] == 'O58462' and 'reviewed' in ref['entryType'].lower() and 'unreviewed' not in ref['entryType'].lower(), 'Wrong reference status')
    seq = ref['sequence']['value']; fp = hashlib.sha256(seq.encode()).hexdigest()
    require(ref['sequence']['length'] == 208, 'Changed reference length')
    require(fp not in {r['protein_sequence_sha256'] for r in original['references'] + original['query_proposals']}, 'Exact copy cannot expand panel')
    new_ref = {'reference_id':'uniprot-O58462', 'accession':'O58462', 'family_proposal':'omp_decarboxylase_homologs', 'original_sha256':digest_file(ref_path), 'sequence_sha256':fp, 'length_aa':208, 'organism':ref['organism']['scientificName'], 'lineage_group':'natural-origin-taxon-70601', 'evidence':'Reviewed Swiss-Prot monofunctional OMP decarboxylase, EC 4.1.1.23; natural archaeal origin rather than a query/vector derivative. Curated functional assertion includes homology evidence, not independent wet-lab validation.', 'source_url':'https://rest.uniprot.org/uniprotkb/O58462.json', 'rights':'UniProt CC BY 4.0; attribution/accession retained; existing license evidence data/manifests/week2/references/reference_proposal.json', 'role':'One prospective reference only; does not add a query parent or independent cluster'}
    elig = {r['query_id']:sorted(r['eligible_reference_ids_proposed'] + [new_ref['reference_id']]) for r in original['eligibility_proposals'] if r['query_id'] in {q['query_id'] for q in original['query_proposals'] if q['proposal_role']=='PROPOSED_SUPPORTED_PENDING_REVIEW'}}
    proposal = {'schema_version':'prospective-panel-2026-11-12-v2', 'status':'HOLD_PENDING_ACTUAL_DELEGATED_REVIEW', 'scientific_clearance':False, 'scored':False, 'original_panel_sha256':digest_file(ROOT/'data/manifests/week2/references/reference_proposal.json'), 'new_E1_only_sources':new_sources, 'new_reference':new_ref, 'reference_ids':sorted(original['reference_ids']+[new_ref['reference_id']]), 'queries':[q for q in original['query_proposals'] if q['proposal_role']=='PROPOSED_SUPPORTED_PENDING_REVIEW'], 'eligible_reference_ids':elig, 'exclusions':'All original exact-copy/ancestor/lineage exclusions retained; Renilla query remains excluded; no stricter tier inferred', 'E2_scope':'Same AF298787 supported subset, unchanged max-per-reference primary. Separately proposed stable round-robin K=5 and component-attributed secondary endpoint. Full-family endpoint remains null.', 'models':{'esm2_revision':'c731040fcd8d73dceaa04b0a8e6329b345b0f5df','weight_sha256':'24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189','alignment':'Original frozen BLOSUM62 configuration'}, 'criteria':{'K':5,'seed':20261105,'precision':'float32','max_residues':1022,'min_independent_clusters_for_CI':8,'primary':'Original documented-family retrieval endpoints; normalization sensitivity remains secondary; no post-score selection'}, 'bounds':{'memory_GiB':6,'wall_seconds_per_model':1800,'unique_proteins':64,'processed_residues':65536,'paid_services':False,'model_downloads':False}, 'limitations':['This is reference-composition sensitivity and E1 architecture expansion, not an independent-query population benchmark.', 'Existing four queries still represent three parents and two connected clusters; CI and general decomposition claims remain unavailable.', 'New fusion sources lack complete independently reviewed family/boundary truth; E2/E3 scored use is not requested for them.']}
    proposal['proposal_sha256'] = hashlib.sha256(json.dumps(proposal,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    audit['supplemental_sources'] = [{k:v for k,v in r.items() if k not in ('validated_cds',)} for r in new_sources]
    audit['additional_independent_reference'] = new_ref
    write_json(output/'independence-audit.json',audit)
    write_json(ROOT/'data/manifests/prospective-2026-11-12/panel.v2.json',proposal)
    print(proposal['proposal_sha256'])

if __name__ == '__main__':
    main()
