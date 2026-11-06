"""Versioned T07 v1.1 configuration; the frozen v1.0 contracts stay unchanged."""
from __future__ import annotations

from pathlib import Path

from .model_adapters import ESM2_REVISION
from .runtime import InputError, read_json, require_text

EXTENSIONS = {'E1', 'E2', 'E3', 'E4'}
MODELS = {'alignment', 'esm2', 'esmc', 'esm3'}


def validate_extension_config(config: dict, expected: str | None = None) -> None:
    if not isinstance(config, dict) or config.get('schema_version') != 't07-extensions-1':
        raise InputError('Expected versioned t07-extensions-1 configuration')
    if config.get('protocol_version') != '1.1':
        raise InputError('Extensions require T07 addendum v1.1, not an edited v1.0 design')
    require_text(config, 'protocol_id')
    extension = config.get('extension')
    if extension not in EXTENSIONS or expected is not None and extension != expected:
        raise InputError('Extension discriminator must match E1, E2, E3 or E4')
    if config.get('design_status') != 'design_frozen_pending_admission':
        raise InputError('Extension design must freeze before source admission/scoring')
    if type(config.get('top_k')) is not int or type(config.get('retrieval_budget')) is not int:
        raise InputError('Retrieval budgets must be integers')
    if config['top_k'] != 5 or config['retrieval_budget'] != 5:
        raise InputError('T07 v1.1 requires five TOTAL returned hits per parent/query')
    for key,low,high in [('max_residues',1,1022),('seed',0,2**32-1),
                         ('min_independent_clusters',2,10000),('bootstrap_replicates',100,10000)]:
        value = config.get(key)
        if type(value) is not int or not low <= value <= high:
            raise InputError('Invalid extension configuration field: ' + key)
    models = config.get('models')
    if not isinstance(models,list) or any(not isinstance(m,str) or m not in MODELS for m in models) or len(models) != len(set(models)):
        raise InputError('models must be unique supported method names')
    if extension in ('E1','E4') and models:
        raise InputError('E1/E4 do not initiate pretrained inference')
    if extension == 'E2' and (not models or 'alignment' in models):
        raise InputError('E2 requires at least one supported fixed protein model')
    if extension == 'E3' and ('alignment' not in models or not {'esm2','esmc'} <= set(models)):
        raise InputError('E3 predeclares alignment, ESM-2 and ESM-C on identical tier panels')
    revisions = config.get('model_revisions',{})
    if not isinstance(revisions,dict) or any(not isinstance(v,str) or not v.strip() for v in revisions.values()):
        raise InputError('Model revisions must be nonempty identifiers; resolve and freeze them before inference')


def load_extension_config(path: str | Path, expected: str | None = None) -> dict:
    config = read_json(Path(path))
    validate_extension_config(config,expected)
    return config


def base_scoring_config(config: dict, experiment: str = 'C') -> dict:
    """Adapt existing frozen-panel utilities without editing their v1.0 contracts.

    E2 uses its own parent panel; E3 uses existing protein-target scorer semantics.
    The extra v1.1 extension field remains part of the canonical protocol hash.
    """
    validate_extension_config(config)
    models = list(config['models'])
    if experiment == 'C' and 'alignment' not in models:
        models.insert(0,'alignment')
    return {'schema_version':'pilot-scoring-1','protocol_id':config['protocol_id'],
        'protocol_version':'1.1','design_status':config['design_status'],
        'experiment':experiment,'extension':config['extension'],
        **{key:config[key] for key in ('top_k','retrieval_budget','seed','min_independent_clusters',
                                       'bootstrap_replicates','max_residues')},
        'models':models,'model_revisions':config.get('model_revisions',{'esm2':ESM2_REVISION}),
        'alignment':{'mode':'local','substitution_matrix':'BLOSUM62','open_gap_score':-10.0,
                     'extend_gap_score':-0.5,'primary':'raw_score','sensitivity':'reference_self_normalized'}}
