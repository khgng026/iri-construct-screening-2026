#!/usr/bin/env python3
"""Governed metadata-first catalogue acquisition; no models or alignments.

Metadata fetching never requests sequences. Sequence fetching is a separately
gated stage that requires specific recorded approval for every accession.
"""
from __future__ import annotations
import argparse
import csv
import datetime as dt
import io
import json
from pathlib import Path
import re
import sys
import time
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode,urlparse,parse_qs
from urllib.request import Request,urlopen

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.session_store import atomic_bytes,atomic_json,append_event,canonical,digest,file_hash


def require(value,message):
    if not value:raise ValueError(message)


def read(path):return json.loads(Path(path).read_text())


def utcnow():return dt.datetime.now(dt.timezone.utc).isoformat()


def frozen_json(path,value):
    """Idempotent freeze, never silently replace a different admitted panel."""
    if path.exists():require(read(path)==value,'Frozen artifact already exists with different contents: '+str(path))
    else:atomic_json(path,value)


def selection_key(accession,design):
    return digest((design['selection_seed']+'|'+accession).encode()),accession


def metadata_eligibility_errors(record):
    """Mechanical identity gates supplement, and cannot replace, source review."""
    errors=[]
    if record.get('reviewed','').lower()!='reviewed':errors.append('not reviewed Swiss-Prot')
    if record.get('organism_id')!='559292':errors.append('not S288c taxon 559292')
    if '-' in record.get('accession',''):errors.append('isoform accession')
    if record.get('fragment','').strip():errors.append('fragment annotation')
    if not record.get('version','').isdigit() or not record.get('sequence_version','').isdigit():errors.append('missing entry/sequence version')
    if not record.get('length','').isdigit() or not 1<=int(record.get('length','0'))<=1022:errors.append('length outside cap')
    if record.get('protein_existence','').lower()=='uncertain':errors.append('uncertain protein existence')
    return errors


def source_family_statements(text):
    """Explicit annotation only; same first-sentence rule as original labels."""
    return sorted(set(re.findall(r'Belongs to the [^.{}]*? family(?=\.|\s*\{|$)',text)))


def nested_prefixes(base_order,new_order,targets):
    require(len(set(base_order+new_order))==len(base_order)+len(new_order),'Representative accession repeated')
    full=base_order+new_order
    return {str(n):{'requested_unique_sequences':n,'achieved_unique_sequences':min(n,len(full)),
                    'endpoint_reached':len(full)>=n,'accessions':full[:n]}
            for n in targets}


def retrieve(url,path,design,ledger):
    """One bounded official public request; original bytes/provenance retained."""
    require(urlparse(url).scheme=='https' and urlparse(url).hostname=='rest.uniprot.org','Unofficial acquisition endpoint')
    provenance=path.with_suffix(path.suffix+'.provenance.json')
    if path.exists() or provenance.exists():
        require(path.exists() and provenance.exists(),'Incomplete previous download; preserve and inspect')
        saved=read(provenance);require(saved['url']==url and file_hash(path)==saved['sha256'],'Download cache identity/corruption mismatch')
        return path.read_bytes(),saved
    attempts=[];network=design['network_policy']
    for attempt in range(network['maximum_retries']+1):
        try:
            request=Request(url,headers={'Accept':'text/tab-separated-values, application/json;q=0.9',
                                         'User-Agent':'IRI-benign-catalogue-review/1.0'})
            with urlopen(request,timeout=network['timeout_seconds']) as response:
                limit=design['metadata_pagination']['maximum_response_bytes'];body=response.read(limit+1)
                require(len(body)<=limit,'Response exceeds bounded acquisition size')
                headers=dict(response.headers);status=response.status
            result={'url':url,'retrieved_utc':utcnow(),'status':status,'headers':headers,
                    'sha256':digest(body),'bytes':len(body),'prior_failed_attempts':attempts}
            atomic_bytes(path,body);atomic_json(provenance,result)
            append_event(ledger,{'event':'download_complete','path':str(path),**result})
            return body,result
        except (HTTPError,URLError,TimeoutError) as exc:
            transient=not isinstance(exc,HTTPError) or exc.code==429 or 500<=exc.code<600
            failure={'attempt':attempt+1,'type':type(exc).__name__,'message':str(exc),'transient':transient}
            attempts.append(failure);append_event(ledger,{'event':'download_failed','url':url,**failure})
            if not transient or attempt==network['maximum_retries']:
                atomic_json(path.with_suffix(path.suffix+'.blocked.json'),{'url':url,'attempts':attempts,'status':'blocked'});raise
            time.sleep(network['retry_backoff_seconds'][attempt])


def field_labels(config):
    labels={}
    def walk(value):
        if isinstance(value,dict):
            if 'name' in value and 'label' in value:labels[value['name']]=value['label']
            for child in value.values():walk(child)
        elif isinstance(value,list):
            for child in value:walk(child)
    walk(config);return labels


def next_link(headers):
    value=next((v for k,v in headers.items() if k.lower()=='link'),'')
    # Commas are legal inside the fields= query value; do not split blindly.
    value=value.replace('\\"','"')
    for match in re.finditer(r'<([^>]+)>;\s*rel="?([^";,\s]+)"?',value):
        if match.group(2)=='next':return match.group(1)
    return None


def setup(args):
    session=args.session.resolve();design=read(args.design);root=session/'expansion';root.mkdir(exist_ok=True)
    require(design['version']=='iri-expanded-catalogue-v1','Unknown expansion design')
    require(file_hash(root/'evidence/result_fields.json')==design['official_field_configuration_sha256'],'Official field configuration changed')
    require(file_hash(root/'evidence/license.json')==design['license_evidence_sha256'],'Source terms evidence changed')
    license_document=read(root/'evidence/license.json')
    require('Creative Commons Attribution 4.0' in license_document['content'],'Missing affirmative reuse licence')
    audit=read(session/'audit/independent_esm2-650m-cpu-full-recovery01.json')
    require(audit['status']=='passed' and audit['query_count']==1146,'Shared larger-model comparison gate not passed')
    require(file_hash(ROOT/'data/manifests/t05_expanded_catalogue.json')==design['source_catalogue_manifest_sha256'],
            'Original catalogue changed')
    return session,root,design


def metadata(args):
    session,root,design=setup(args);destination=root/'metadata';destination.mkdir(exist_ok=True)
    fields=design['metadata_fields'];require('sequence' not in fields,'Metadata stage cannot request sequence')
    labels=field_labels(read(root/'evidence/result_fields.json'))
    require(all(field in labels for field in fields),'Unknown current official metadata field')
    expected_headers=[labels[field] for field in fields]
    pagination=design['metadata_pagination']
    params={'query':design['metadata_query'],'format':'tsv','fields':','.join(fields),'size':pagination['size'],
            'sort':pagination['sort'],'includeIsoform':'false'}
    url=design['primary_source']+'?'+urlencode(params);records={};pages=[];releases=set();complete=False;total_results=None
    for page in range(pagination['maximum_pages']):
        parsed=parse_qs(urlparse(url).query)
        require(parsed.get('fields')==[','.join(fields)] and parsed.get('format')==['tsv'],'Pagination altered metadata field policy')
        body,provenance=retrieve(url,destination/f'page-{page:03d}.tsv',design,root/'acquisition_ledger.jsonl')
        reader=csv.DictReader(io.StringIO(body.decode('utf-8')),delimiter='\t')
        require(reader.fieldnames==expected_headers and 'Sequence' not in reader.fieldnames,
                'Metadata response fields differ or expose sequence column')
        page_records=list(reader)
        for index,row in enumerate(page_records):
            require(set(row)==set(expected_headers) and None not in row,'Malformed metadata row')
            record={field:row[labels[field]] for field in fields};accession=record['accession']
            require(accession not in records,'Duplicate accession across paginated metadata; release may have changed')
            records[accession]={'metadata':record,'metadata_row_sha256':digest(canonical(record)),
                                'source_page':str((destination/f'page-{page:03d}.tsv').relative_to(session)),
                                'source_page_sha256':provenance['sha256'],'source_row_number_1based':index+2,
                                'source_url':provenance['url'],'source_retrieved_utc':provenance['retrieved_utc'],
                                'review_status':'pending','safety_review':'pending','licensing_status':'pending'}
        release=next((v for k,v in provenance['headers'].items() if k.lower()=='x-uniprot-release'),None)
        declared=next((v for k,v in provenance['headers'].items() if k.lower()=='x-total-results'),None)
        require(declared is not None,'Metadata pagination lacks total-result count')
        if total_results is None:total_results=int(declared)
        require(int(declared)==total_results,'Source total changed during pagination')
        if release:releases.add(release)
        require(len(releases)<=1,'Release changed during pagination; freeze cannot mix releases')
        pages.append({'page':page,'records':len(page_records),'response_sha256':provenance['sha256'],'release':release})
        print(json.dumps({'stage':'metadata','page':page+1,'records':len(records),'release':release}),flush=True)
        url=next_link(provenance['headers'])
        if url is None:
            require(len(records)==total_results,'Pagination ended before declared total; refuse partial completeness')
            complete=True;break
    status='complete_metadata_pending_review' if complete else 'blocked_page_budget'
    result={'status':status,'design_sha256':file_hash(args.design),'records':records,'record_count':len(records),
            'pages':pages,'source_releases':sorted(releases),'source_total_results':total_results,'metadata_only':True,'sequences_requested':False,
            'remaining_url':url,'review_instruction':'Every record remains pending until its own source annotations and provenance receive a recorded review.'}
    path=root/'metadata_index.json'
    if path.exists():require(read(path)==result,'Refusing to overwrite differing metadata freeze')
    else:atomic_json(path,result)
    print(json.dumps({'status':status,'record_count':len(records),'path':str(path),'sha256':file_hash(path)}),flush=True)


def sequence_candidates(metadata_index,review,design,base_accessions):
    require(review['metadata_index_sha256']==metadata_index['file_sha256'],'Review refers to different source metadata')
    approved=[]
    for accession,decision in review['records'].items():
        require(accession in metadata_index['records'],'Review accession absent from frozen metadata')
        source=metadata_index['records'][accession]
        require(decision['metadata_row_sha256']==source['metadata_row_sha256'],'Reviewed metadata row changed')
        if all(decision.get(k)=='approved' for k in ('review_status','safety_review','licensing_status')):
            require(decision.get('rationale') and decision.get('source_annotation_basis'),'Approved record lacks individual evidence/rationale')
            require(decision['license_evidence_sha256']==design['license_evidence_sha256'],'Individual licence evidence mismatch')
            require(not metadata_eligibility_errors(source['metadata']),'Approved metadata fails mechanical source gate: '+accession)
            if accession not in base_accessions:approved.append(accession)
    return sorted(approved,key=lambda a:selection_key(a,design))


def verify_metadata_sources(session,index,design):
    """Rebuild every indexed source row from its exact downloaded page bytes."""
    labels=field_labels(read(session/'expansion/evidence/result_fields.json'));fields=design['metadata_fields']
    cached={}
    for accession,record in index['records'].items():
        source=session/record['source_page']
        if source not in cached:
            provenance=read(source.with_suffix(source.suffix+'.provenance.json'))
            require(file_hash(source)==provenance['sha256'],'Metadata response bytes changed')
            reader=csv.DictReader(io.StringIO(source.read_text()),delimiter='\t')
            require(reader.fieldnames==[labels[f] for f in fields],'Frozen metadata response header changed')
            cached[source]=(provenance,list(reader))
        provenance,rows=cached[source]
        require(provenance['sha256']==record['source_page_sha256'] and provenance['url']==record['source_url'],
                'Metadata provenance identity changed')
        raw=rows[record['source_row_number_1based']-2];parsed={field:raw[labels[field]] for field in fields}
        require(parsed==record['metadata'] and parsed['accession']==accession and
                digest(canonical(parsed))==record['metadata_row_sha256'],'Indexed metadata row no longer matches original bytes')
    require(len(index['records'])==index['record_count']==index['source_total_results'],'Frozen metadata count mismatch')


def compile_review(args):
    """Bind explicit per-accession reviewer decisions to immutable source rows."""
    session,root,design=setup(args)
    require(args.decision_file and args.review is not None,'Review compilation requires explicit decision files and destination')
    index_path=root/'metadata_index.json';index=read(index_path)
    require(index['status']=='complete_metadata_pending_review' and index['design_sha256']==file_hash(args.design),
            'Source metadata is incomplete or from another design')
    verify_metadata_sources(session,index,design);records={};inputs=[]
    snapshot_dir=root/'review_inputs';snapshot_dir.mkdir(exist_ok=True)
    for path in args.decision_file:
        reader=csv.DictReader(io.StringIO(path.read_text()),delimiter='\t')
        require(reader.fieldnames==['accession','decision','ordinary_role_basis','review_note'],'Manual decision schema differs')
        sha=file_hash(path);snapshot=snapshot_dir/(sha+'.tsv')
        if snapshot.exists():require(file_hash(snapshot)==sha,'Manual review snapshot changed')
        else:atomic_bytes(snapshot,path.read_bytes())
        inputs.append({'path':str(snapshot.relative_to(session)),'sha256':sha})
        for row in reader:
            accession=row['accession'];status=row['decision']
            require(accession in index['records'] and accession not in records,'Unknown or multiply reviewed accession')
            require(status in ('approved','pending','uncertain','rejected') and row['ordinary_role_basis'] and row['review_note'],
                    'Missing explicit individual review decision/rationale')
            source=index['records'][accession];metadata=source['metadata']
            if status=='approved':require(not metadata_eligibility_errors(metadata),'Approved row violates source identity gate')
            records[accession]={'accession':accession,'metadata_row_sha256':source['metadata_row_sha256'],
                 'source_page':source['source_page'],'source_page_sha256':source['source_page_sha256'],
                 'source_url':source['source_url'],'source_retrieved_utc':source['source_retrieved_utc'],
                 'review_status':status,'safety_review':status,'licensing_status':'approved',
                 'annotation_status':'known' if status=='approved' else 'uncertain',
                 'rationale':row['ordinary_role_basis'],'review_note':row['review_note'],
                 'source_annotation_basis':{k:metadata[k] for k in ('protein_name','gene_names','cc_function','cc_similarity',
                                            'protein_families','cc_caution','cc_sequence_caution','keyword')},
                 'source_identity_basis':{k:metadata[k] for k in ('reviewed','organism_id','length','fragment','protein_existence',
                                         'version','sequence_version','date_modified','date_sequence_modified')},
                 'license_evidence_sha256':design['license_evidence_sha256'],'license_url':design['license_url'],
                 'reuse_scope':design['license'],'reviewer':'Codex source-annotation review under authorized benign-only scope; not institutional approval',
                 'manual_decision_source_sha256':inputs[-1]['sha256']}
    require(records,'No individual records reviewed')
    prior_manifests=[]
    for prior_path in args.prior_review or []:
        prior=read(prior_path)
        require(prior['metadata_index_sha256']==file_hash(index_path) and prior['design_sha256']==file_hash(args.design),
                'Prior milestone reviewed different metadata/design')
        for accession,previous in prior['records'].items():
            require(accession in records,'Previously reviewed record missing from extended review')
            for key in ('review_status','safety_review','licensing_status','metadata_row_sha256','rationale','review_note',
                        'source_annotation_basis','source_identity_basis','license_evidence_sha256'):
                require(records[accession][key]==previous[key],'Earlier frozen review decision changed')
        previous_max=max(selection_key(a,design) for a in prior['records'])
        require(all(selection_key(a,design)>previous_max for a in set(records)-set(prior['records'])),
                'New review would precede an earlier frozen hash-order milestone')
        prior_manifests.append({'path':str(prior_path.resolve()),'sha256':file_hash(prior_path)})
    result={'version':'iri-expanded-admission-review-v1','status':'frozen_before_new_sequence_acquisition',
        'design_sha256':file_hash(args.design),'metadata_index_sha256':file_hash(index_path),'decision_files':inputs,
        'records':records,'prior_frozen_reviews':prior_manifests,
        'counts':{s:sum(r['review_status']==s for r in records.values()) for s in ('approved','pending','uncertain','rejected')},
        'unreviewed_metadata_records':len(index['records'])-len(records),
        'unreviewed_policy':'Every unlisted metadata accession remains pending and is ineligible for sequence acquisition.',
        'limitations':'Source annotation review establishes benign inclusion only, not validated protein function, screening ability or institutional approval.'}
    frozen_json(args.review,result)
    print(json.dumps({'status':result['status'],'counts':result['counts'],'sha256':file_hash(args.review),'path':str(args.review)}),flush=True)


def sequences(args):
    # Deliberately separate from metadata and review. This stage cannot create
    # approvals and never fills pending or rejected records automatically.
    session,root,design=setup(args)
    require(args.review is not None,'Sequence acquisition requires a specific approved review manifest')
    metadata_path=root/'metadata_index.json';index=read(metadata_path)
    require(index['status']=='complete_metadata_pending_review','Metadata catalogue is incomplete')
    require(index['design_sha256']==file_hash(args.design),'Metadata was acquired under another design')
    verify_metadata_sources(session,index,design)
    index['file_sha256']=file_hash(metadata_path);review=read(args.review)
    require(review.get('status')=='frozen_before_new_sequence_acquisition' and review.get('design_sha256')==file_hash(args.design),
            'Individual admission review is not independently frozen under this design')
    base=read(ROOT/'data/manifests/t05_expanded_catalogue.json');base_ids={r['accession'] for r in base['records']}
    ordered=sequence_candidates(index,review,design,base_ids)[:design['maximum_new_sequence_records']]
    require(ordered,'No specifically approved new records to acquire')
    if args.limit is not None:ordered=ordered[:args.limit]
    output=root/'sequences';output.mkdir(exist_ok=True);accepted={};rejected={};seen={r['sequence_sha256'] for r in base['records']}
    target=args.target or max(design['nested_unique_sequence_targets'])
    require(target in design['nested_unique_sequence_targets'],'Unplanned target endpoint')
    milestone=root/('milestone-'+str(target));milestone.mkdir(exist_ok=True);batchsize=design['sequence_batch_size']
    require(not (milestone/'admitted/expanded_catalogue_v1.json').exists(),'Milestone already frozen; cannot rewrite its acquisition progress')
    for offset in range(0,len(ordered),batchsize):
        batch=ordered[offset:offset+batchsize]
        query=' OR '.join('accession:'+a for a in batch)
        fields=design['metadata_fields']+['sequence']
        url=design['primary_source']+'?'+urlencode({'query':'('+query+')','format':'tsv','fields':','.join(fields),
                                                  'size':batchsize,'includeIsoform':'false','sort':'accession asc'})
        identifier=digest(canonical(batch))[:16]
        body,provenance=retrieve(url,output/('batch-'+identifier+'.tsv'),design,root/'acquisition_ledger.jsonl')
        labels=field_labels(read(root/'evidence/result_fields.json'));reader=csv.DictReader(io.StringIO(body.decode()),delimiter='\t')
        require(reader.fieldnames==[labels[f] for f in fields],'Sequence response field policy changed')
        rows={row[labels['accession']]:row for row in reader};require(set(rows)==set(batch),'Sequence batch omitted/added accessions')
        for accession in batch:
            original=index['records'][accession]['metadata'];row=rows[accession]
            current={field:row[labels[field]] for field in design['metadata_fields']}
            sequence=row[labels['sequence']]
            reason=None
            if current!=original:reason='Source metadata/version differs from approved review; re-review required'
            elif not sequence or not set(sequence)<=set('ACDEFGHIKLMNPQRSTVWY'):reason='Noncanonical/unsupported sequence alphabet; unchanged and excluded'
            elif len(sequence)!=int(original['length']) or len(sequence)>1022:reason='Sequence length mismatch/over cap; no truncation'
            sha=digest(sequence.encode())
            if reason:rejected[accession]={'reason':reason,'response_sha256':provenance['sha256']};continue
            duplicate=sha in seen
            accepted[accession]={'accession':accession,'sequence_sha256':sha,'length':len(sequence),
                                'duplicate_of_previously_admitted_sequence':duplicate,'metadata':original,
                                'metadata_row_sha256':index['records'][accession]['metadata_row_sha256'],
                                'review':review['records'][accession],'source_response':str((output/('batch-'+identifier+'.tsv')).relative_to(session)),
                                'source_response_sha256':provenance['sha256'],'sequence':sequence}
            seen.add(sha)
        atomic_json(milestone/'sequence_acquisition_progress.json',{'design_sha256':file_hash(args.design),
             'review_sha256':file_hash(args.review),'accepted':accepted,'rejected':rejected,
             'requested_milestone':target,'unique_sequences_including_base':len(seen),'processed_approved_records':offset+len(batch)})
        print(json.dumps({'stage':'sequences','processed':offset+len(batch),'unique_including_base':len(seen),'rejected':len(rejected)}),flush=True)
        if len(seen)>=target:break
    # The sequence-bearing progress artifact remains local/ignored. Final nested
    # catalogue construction is a separate frozen manifest step after review.


def finalize(args):
    """Freeze achieved, nested unique references after exact saved-byte checks."""
    session,root,design=setup(args)
    require(args.review is not None,'Finalization requires the exact source review manifest')
    index_path=root/'metadata_index.json';index=read(index_path)
    require(index['status']=='complete_metadata_pending_review' and index['design_sha256']==file_hash(args.design),
            'Metadata completeness/design mismatch')
    verify_metadata_sources(session,index,design);index['file_sha256']=file_hash(index_path)
    review=read(args.review)
    require(review.get('status')=='frozen_before_new_sequence_acquisition' and review.get('design_sha256')==file_hash(args.design),
            'Individual admission review is not frozen under this design')
    base=read(ROOT/'data/manifests/t05_expanded_catalogue.json')
    base_ids={r['accession'] for r in base['records']}
    ordered=sequence_candidates(index,review,design,base_ids)[:design['maximum_new_sequence_records']]
    target=args.target or max(design['nested_unique_sequence_targets'])
    require(target in design['nested_unique_sequence_targets'],'Unplanned target endpoint')
    milestone=root/('milestone-'+str(target))
    progress_path=milestone/'sequence_acquisition_progress.json';progress=read(progress_path)
    require(progress['design_sha256']==file_hash(args.design) and progress['review_sha256']==file_hash(args.review),
            'Acquisition used a different design/review')
    require(progress['requested_milestone']==target,'Acquisition milestone differs')
    processed=progress['processed_approved_records'];expected=ordered[:processed]
    require(set(progress['accepted']).isdisjoint(progress['rejected']) and
            set(progress['accepted'])|set(progress['rejected'])==set(expected),'Acquisition does not cover a deterministic approval prefix')
    require(file_hash(session/'inputs/challenge_manifest.json')==design['source_challenge_manifest_sha256'],'Frozen queries changed')
    challenge=read(session/'inputs/challenge_manifest.json')
    targets=sorted({a for q in challenge['queries'] for a in q['components']})
    require(len(base_ids)==105 and len(targets)==16 and set(targets)<=base_ids,'Original catalogue/target membership changed')
    families=read(session/'inputs/family_labels_v1.json');catalogue={};by_hash={};aliases={};sequences_by_rep={}
    for record in sorted(base['records'],key=lambda r:r['accession']):
        accession=record['accession'];sha=record['sequence_sha256']
        require(all(record[k]=='approved' for k in ('review_status','safety_review','licensing_status')),'Base source lacks approval')
        require(file_hash(ROOT/record['raw_path'])==record['checksum_sha256'] and
                file_hash(ROOT/record['source_page_path'])==record['source_page_sha256'],'Base source evidence changed')
        fasta=ROOT/record['sequence_ref']['path'];lines=fasta.read_text().splitlines()
        require(lines and lines[0].startswith('>') and sum(line.startswith('>') for line in lines)==1,'Malformed base FASTA')
        sequence=''.join(lines[1:]);require(digest(sequence.encode())==sha and len(sequence)==record['sequence_length_aa'],
                                          'Base canonical sequence changed')
        representative=by_hash.setdefault(sha,accession);aliases[accession]=representative
        if representative not in catalogue:
            catalogue[representative]={'accession':representative,'sequence_sha256':sha,'length':len(sequence),
                 'origin':'original_reviewed_105','aliases':[],'original_source':record,
                 'explicit_family_statements':sorted({e['family_statement'] for e in families['source_evidence'][accession]['annotation_evidence']})}
            sequences_by_rep[representative]=sequence
        catalogue[representative]['aliases'].append(accession)
    base_order=sorted(catalogue);require(len(base_order)==101,'Original unique sequence count changed')
    labels=field_labels(read(root/'evidence/result_fields.json'));response_rows={};new_order=[];new_duplicates={}
    for accession in expected:
        if accession in progress['rejected']:continue
        accepted=progress['accepted'][accession];source=session/accepted['source_response']
        if source not in response_rows:
            provenance=read(source.with_suffix(source.suffix+'.provenance.json'))
            require(file_hash(source)==provenance['sha256'],'Saved sequence response changed')
            reader=csv.DictReader(io.StringIO(source.read_text()),delimiter='\t')
            require(reader.fieldnames==[labels[f] for f in design['metadata_fields']+['sequence']],'Sequence response fields changed')
            rows=list(reader);require(len({r[labels['accession']] for r in rows})==len(rows),'Duplicate response accessions')
            response_rows[source]=(provenance,{r[labels['accession']]:r for r in rows})
        provenance,rows=response_rows[source];raw=rows[accession]
        metadata={field:raw[labels[field]] for field in design['metadata_fields']};sequence=raw[labels['sequence']]
        require(provenance['sha256']==accepted['source_response_sha256'] and metadata==index['records'][accession]['metadata'] and
                metadata==accepted['metadata'] and sequence==accepted['sequence'] and digest(sequence.encode())==accepted['sequence_sha256'],
                'Accepted sequence/provenance differs from the approved source')
        require(set(sequence)<=set('ACDEFGHIKLMNPQRSTVWY') and len(sequence)==int(metadata['length'])<=1022,
                'Accepted sequence violates canonical alphabet/length gate')
        sha=accepted['sequence_sha256'];representative=by_hash.setdefault(sha,accession);aliases[accession]=representative
        if representative!=accession:
            new_duplicates[accession]=representative;catalogue[representative]['aliases'].append(accession);continue
        catalogue[accession]={'accession':accession,'sequence_sha256':sha,'length':len(sequence),'origin':'new_reviewed_metadata',
             'aliases':[accession],'metadata':metadata,'review':review['records'][accession],
             'metadata_row_sha256':accepted['metadata_row_sha256'],'source_response':accepted['source_response'],
             'source_response_sha256':accepted['source_response_sha256'],'source_provenance':provenance,
             'explicit_family_statements':source_family_statements(metadata['cc_similarity'])}
        new_order.append(accession);sequences_by_rep[accession]=sequence
    require(len(by_hash)==progress['unique_sequences_including_base'],'Final unique count differs from acquisition')
    full_order=base_order+new_order
    subsets=nested_prefixes(base_order,new_order,[n for n in design['nested_unique_sequence_targets'] if n<=target])
    earlier_milestones=[]
    for earlier in design['nested_unique_sequence_targets']:
        prior_path=root/('milestone-'+str(earlier))/'admitted/expanded_catalogue_v1.json'
        if earlier>=target or not prior_path.exists():continue
        prior=read(prior_path)
        require(prior['design_sha256']==file_hash(args.design) and prior['metadata_index_sha256']==file_hash(index_path),
                'Earlier admitted milestone source/design differs')
        for size,panel in prior['nested_catalogues'].items():
            previous=panel['accessions'];current=subsets[size]['accessions']
            require(current[:len(previous)]==previous and (not panel['endpoint_reached'] or current==previous),
                    'Previously frozen nested endpoint changed')
            require(all(catalogue[a]['sequence_sha256']==prior['catalogue'][a]['sequence_sha256'] for a in previous),
                    'Previously frozen representative sequence changed')
        earlier_milestones.append({'path':str(prior_path.relative_to(session)),'sha256':file_hash(prior_path)})
    target_families={a:sorted({e['family_statement'] for e in families['source_evidence'][a]['annotation_evidence']}) for a in targets}
    family_composition={size:{target:[a for a in panel['accessions'] if
                          set(catalogue[a]['explicit_family_statements'])&set(statements)]
                          for target,statements in target_families.items()} for size,panel in subsets.items()}
    output=milestone/'admitted';output.mkdir(exist_ok=True)
    fasta_bytes=''.join('>'+a+'\n'+sequences_by_rep[a]+'\n' for a in full_order).encode()
    fasta_path=output/'canonical_references.fasta'
    if fasta_path.exists():require(fasta_path.read_bytes()==fasta_bytes,'Frozen canonical FASTA changed')
    else:atomic_bytes(fasta_path,fasta_bytes)
    manifest={'version':'iri-expanded-admitted-v1','status':'frozen_before_expanded_scoring',
       'design_sha256':file_hash(args.design),'metadata_index_sha256':file_hash(index_path),'review_sha256':file_hash(args.review),
       'sequence_acquisition_progress':str(progress_path.relative_to(session)),
       'sequence_acquisition_progress_sha256':file_hash(progress_path),'acquisition_script_sha256':file_hash(Path(__file__)),
       'source_catalogue_manifest_sha256':design['source_catalogue_manifest_sha256'],
       'original_family_labels_sha256':file_hash(session/'inputs/family_labels_v1.json'),
       'source_releases':index['source_releases'],'source_base_note':'Original reviewed105 source versions remain unchanged; new records use the separately recorded current release.',
       'catalogue':catalogue,'base_representatives':base_order,'new_representatives_hash_order':new_order,
       'aliases_to_representative':aliases,'target_aliases_to_representative':{a:aliases[a] for a in targets},
       'nested_catalogues':subsets,'new_exact_sequence_duplicates':new_duplicates,
       'requested_endpoints':design['nested_unique_sequence_targets'],'current_milestone':target,
       'earlier_frozen_milestones':earlier_milestones,
       'later_endpoints_pending_review':[n for n in design['nested_unique_sequence_targets'] if n>target],
       'total_unique_sequences':len(catalogue),
       'reviewed_approved_new_records':len(ordered),'acquired_new_records':processed,'rejected_after_acquisition':progress['rejected'],
       'unacquired_approved_records':len(ordered)-processed,'unreviewed_metadata_records':len(index['records'])-len(review['records']),
       'review_decision_counts':{status:sum(r.get('review_status')==status for r in review['records'].values())
                                for status in ('approved','pending','uncertain','rejected')},
       'target_explicit_families':target_families,'target_family_composition_by_subset':family_composition,
       'relatedness_interpretation':'Descriptive exact source-family overlap, not functional validation or a pure catalogue-size effect; no scores informed selection.',
       'canonical_fasta':str(fasta_path.relative_to(session)),'canonical_fasta_sha256':file_hash(fasta_path),
       'sequence_transformations':'none; canonical strings preserved; exact SHA256 duplicates represented once with aliases',
       'attribution':design['attribution'],'local_ignored_only':True}
    path=output/'expanded_catalogue_v1.json';frozen_json(path,manifest)
    print(json.dumps({'status':manifest['status'],'unique_sequences':len(catalogue),'nested_counts':
                      {k:v['achieved_unique_sequences'] for k,v in subsets.items()},'path':str(path),'sha256':file_hash(path)}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('stage',choices=['metadata','review','sequences','finalize'])
    parser.add_argument('--session',type=Path,required=True);parser.add_argument('--design',type=Path,required=True)
    parser.add_argument('--review',type=Path);parser.add_argument('--limit',type=int);parser.add_argument('--target',type=int)
    parser.add_argument('--decision-file',type=Path,action='append')
    parser.add_argument('--prior-review',type=Path,action='append')
    args=parser.parse_args();{'metadata':metadata,'review':compile_review,'sequences':sequences,'finalize':finalize}[args.stage](args)
