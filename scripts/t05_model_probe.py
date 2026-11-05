"""ESM-2 feasibility probe on explicitly reviewed local protein inputs."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.experiment_cli import main

if __name__ == '__main__':
    raise SystemExit(main(default_experiment='probe'))
