#!/usr/bin/env python3
"""Audit up to fourteen publicly published, benign iGEM candidates without credentials.

This script does NOT grant scientific admission. It stores downloaded JSON only
under ignored results, and publishes a source-free accession/annotation audit.
"""
import datetime
import hashlib
import json
from pathlib import Path
import re
import time
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

ROOT=Path(__file__).resolve().parents[1]
INVENTORY=ROOT/"data/manifests/igem_candidates_2026_11_10.json"
OUT=ROOT/"results/igem_audit_2026_11_10"
BASE="https://api.registry.igem.org/v1"

def get_json(path):
    u=BASE+path
    req=Request(u,headers={"User-Agent":"IRI-Construct-Screening/0.2 (public academic metadata review; no authentication)","Accept":"application/json"})
    raw=None
    for attempt in range(7):
        try:
            with urlopen(req,timeout=30) as h:
                raw=h.read(1_500_001)
            break
        except HTTPError as e:
            if e.code != 429 or attempt == 6:
                raise
            retry=e.headers.get("Retry-After", "")
            delay=min(120, max(40, int(retry) if retry.isdigit() else (attempt + 1)*40))
            print("Registry rate-limited (HTTP 429); respecting cooldown",delay,"seconds",flush=True)
            time.sleep(delay)
    if raw is None:
        raise ValueError("No source bytes received")
    if len(raw)>1_500_000:
        raise ValueError("Oversized response")
    data=json.loads(raw)
    return data,raw,u

def secret_free_error(e):
    if isinstance(e,HTTPError): return f"HTTP {e.code}"
    return f"{type(e).__name__}: {str(e)[:130]}"

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    raw=OUT/"raw"
    raw.mkdir(exist_ok=True)
    items=json.loads(INVENTORY.read_text())["candidates"]
    report=[]
    for x in items:
        slug=x["slug"]
        entry={"part_id":x["part_id"],"slug":slug,"source_url":x["source_url"],
               "archetype_proposed":x["proposed_archetype"],"status":"HOLD",
               "source_admission":False,"model_scoring_approved":False}
        try:
            info,b,u=get_json("/parts/slugs/"+quote(slug))
            (raw/(slug+".json")).write_bytes(b)
            seq=info.get("sequence")
            if isinstance(seq,dict):seq=seq.get("sequence") or seq.get("elements")
            if not isinstance(seq,str) or not re.fullmatch("[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+",seq):
                raise ValueError("No complete nucleotide sequence in anonymous public part response")
            if len(seq)<400: raise ValueError("Sequence too short for reviewed construct candidate")
            status=str(info.get("status","")).lower()
            uuid=info.get("uuid") or info.get("id")
            entry.update(title=info.get("title"),status="ACCESSIBLE_PUBLIC_SEQUENCE",
                         published_status=status,sequence_length_nt=len(seq),
                         source_json_sha256=hashlib.sha256(b).hexdigest(),
                         source_sequence_sha256=hashlib.sha256(seq.upper().encode()).hexdigest(),
                         source_api=u,
                         license_metadata=str(info.get("licenseUUID") or info.get("license") or "")[:160],
                         part_type=str(info.get("type") or info.get("typeUUID") or "")[:180],
                         uuid=uuid)
            if status!="published":
                entry["status"]="HOLD_NOT_PUBLISHED"
            if uuid:
                try:
                    ann,b2,annurl=get_json("/parts/"+quote(str(uuid))+"/sequence-features?page=1&pageSize=100")
                    (raw/(slug+".annotations.json")).write_bytes(b2)
                    feats=ann.get("data",[]) if isinstance(ann,dict) else ann if isinstance(ann,list) else []
                    entry.update(annotation_count=len(feats),annotation_source_url=annurl,
                                 annotation_response_sha256=hashlib.sha256(b2).hexdigest())
                    if isinstance(ann,dict) and ann.get("total",0)>len(feats):
                        entry["annotation_pagination_incomplete"]=True
                except (HTTPError,URLError,ValueError,TimeoutError) as e:
                    entry["annotation_error"]=secret_free_error(e)
            if entry["status"]=="ACCESSIBLE_PUBLIC_SEQUENCE":
                entry["status"]="DOWNLOADED_REVIEW_REQUIRED"
        except (HTTPError,URLError,ValueError,TimeoutError,KeyError,TypeError) as e:
            entry["error"]=secret_free_error(e)
        report.append(entry)
        print(slug,entry["status"],entry.get("sequence_length_nt","-"),entry.get("annotation_count","-"),flush=True)
        time.sleep(13.0)
    output={"audited_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "candidate_count":len(report),
            "public_sequence_download_count":sum(z["status"]=="DOWNLOADED_REVIEW_REQUIRED" for z in report),
            "no_login_access_verified":True,
            "contains_raw_sequences":False,
            "scientific_admission_count":0,
            "warning":"Not an independent safety, license, feature-level, family or model-scoring review. Public part!=complete plasmid; sequence-features and composite subparts require manual interpretation.",
            "records":report}
    (OUT/"report.json").write_text(json.dumps(output,indent=2)+"\n")
    print("PUBLIC SEQUENCE FILES",output["public_sequence_download_count"],"of",len(report),flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
