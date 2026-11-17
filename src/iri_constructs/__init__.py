"""Annotation-guided preprocessing of annotated DNA constructs.

The public toolkit API is re-exported here. Historical experiment modules stay
importable as submodules but are not part of the supported interface.
"""
__version__ = '0.1.0'

_PUBLIC = ('BatchResult', 'FeatureResult', 'InputFile', 'Issue', 'OutputError', 'PreprocessError',
           'PreprocessOptions', 'RecordResult', 'UnitResult', 'preprocess_file', 'preprocess_files',
           'process_inputs', 'write_bundle')
__all__ = ['__version__', *_PUBLIC]


def __getattr__(name):
    # Lazy so that importing a historical submodule never pays for (or depends on) the toolkit.
    if name in _PUBLIC:
        from . import preprocess
        return getattr(preprocess, name)
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
