"""Dependency-free smoke check using fictional metadata only."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from iri_constructs.config import configure_logging, load_config
from iri_constructs.schema_check import check_schema_example
from iri_constructs.validate_manifest import validate_manifest


def main() -> int:
    try:
        cfg = load_config(ROOT / "config" / "default.json")
        log = configure_logging(cfg["log_level"])
        log.info("event=smoke_start project=%s", cfg["project"])
        issues = validate_manifest(ROOT / cfg["manifest"])
        examples = sorted((ROOT / cfg["schema_examples"]).glob("*.json"))
        if len(examples) < 3:
            issues.append("Expected three fictional schema examples")
        for path in examples:
            issues.extend(f"{path.name}: {issue}" for issue in check_schema_example(path))
        if issues:
            for issue in issues:
                log.error("event=metadata_validation_failed detail=%s", issue)
            return 1
        log.info("event=smoke_pass examples=%d", len(examples))
        print("Week 1 smoke check passed (fictional metadata only).")
        return 0
    except (OSError, ValueError, TypeError) as exc:
        print(f"Configuration or smoke-check error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
