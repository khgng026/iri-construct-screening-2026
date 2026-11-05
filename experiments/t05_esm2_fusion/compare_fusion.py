"""Path-independent entry point; see README.md for reviewed-input configuration."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
from iri_constructs.experiment_cli import main

if __name__ == '__main__':
    raise SystemExit(main(default_experiment='fusion'))
