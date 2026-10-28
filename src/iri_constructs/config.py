"""Standard-library configuration and privacy-conscious logging for metadata checks."""

from __future__ import annotations

import json
import logging
from pathlib import Path

VALID_LOG_LEVELS = {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}
REQUIRED_KEYS = {"project", "manifest", "log_level", "schema_examples"}


def load_config(path: Path) -> dict[str, str]:
    """Load a small repo-relative configuration without shell/env interpretation."""
    path = Path(path)
    settings = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(settings, dict):
        raise ValueError("Configuration must be a JSON object")
    missing = REQUIRED_KEYS - settings.keys()
    if missing:
        raise ValueError(f"Missing configuration settings: {', '.join(sorted(missing))}")
    for key in REQUIRED_KEYS:
        if not isinstance(settings[key], str) or not settings[key].strip():
            raise ValueError(f"Invalid configuration setting: {key}")
    if settings["log_level"] not in VALID_LOG_LEVELS:
        raise ValueError("log_level must be DEBUG, INFO, WARNING, ERROR, or CRITICAL")
    return settings


def configure_logging(level: str) -> logging.Logger:
    """Log metadata and outcomes, never sequence content or customer data."""
    logging.basicConfig(
        level=level,
        format="%(levelname)s %(name)s %(message)s",
        force=True,
    )
    return logging.getLogger("iri_constructs")
