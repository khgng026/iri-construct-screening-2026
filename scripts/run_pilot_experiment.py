"""Run a frozen reviewed-input pilot experiment, or record explicit admission blockers."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.pilot_experiments import main

if __name__ == '__main__':
    raise SystemExit(main())
