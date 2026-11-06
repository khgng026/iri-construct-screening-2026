"""Execute T07 v1.1 E1-E4 with the unchanged source-review restrictions."""
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from iri_constructs.extension_experiments import main

if __name__=='__main__':
    raise SystemExit(main())
