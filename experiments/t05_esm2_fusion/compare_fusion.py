"""Path-independent entry point; see README.md for reviewed-input configuration."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
from iri_constructs.experiment_cli import main
from iri_constructs.fusion_legacy import (
    AA_ALPHABET, REQUIRED, OPTIONAL, parse_fasta, segment_bounds, validate_pair,
    load_pairs, compare_pair, average, digest, cosine,
)

if __name__ == '__main__':
    raise SystemExit(main(default_experiment='fusion'))
