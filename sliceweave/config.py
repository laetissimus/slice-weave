"""Settings. Stored as JSON in the data folder; environment variables override.

Nothing here is secret and nothing identifies a patient.
"""
from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path

# Version = the GitHub release tag (v0.1.1 -> "0.1.1"). The build stamps it into _version.py: from the
# release tag when the release workflow builds, otherwise from `git describe` ("0.1.1+3.g1a2b3c4": 3
# commits after v0.1.1). Running from source without a stamp shows the last release + "+dev".
BASE_VERSION = "0.1.1"


def _version() -> str:
    try:
        from ._version import VERSION as stamped  # written by packaging/build_windows.py
        return str(stamped)
    except ImportError:
        return BASE_VERSION + "+dev"


VERSION = _version()
NAME = "SliceWeave"
RESEARCH_BANNER = "RESEARCH USE ONLY - NOT FOR CLINICAL DIAGNOSIS"


@dataclass
class Settings:
    data_dir: str = str(Path.home() / "SliceWeave-Data")
    port: int = 8766
    retention_days: int = 30
    max_upload_gb: float = 8.0
    max_voxels: int = 30_000_000

    @property
    def root(self) -> Path:
        return Path(self.data_dir).expanduser()

    @property
    def cases_dir(self) -> Path:
        return self.root / "cases"

    @property
    def tmp_dir(self) -> Path:
        return self.root / "tmp"

    def ensure(self) -> "Settings":
        for d in (self.root, self.cases_dir, self.tmp_dir):
            d.mkdir(parents=True, exist_ok=True)
            try:
                os.chmod(d, 0o700)
            except OSError:
                pass
        return self

    def save(self):
        self.ensure()
        (self.root / "settings.json").write_text(json.dumps(asdict(self), indent=2), encoding="utf-8")


def load(data_dir: str | None = None) -> Settings:
    d = data_dir or os.environ.get("SLICEWEAVE_DATA") or Settings().data_dir
    s = Settings(data_dir=str(Path(d).expanduser()))
    p = s.root / "settings.json"
    if p.exists():
        try:
            raw = json.loads(p.read_text(encoding="utf-8"))
            for k, v in raw.items():
                if hasattr(s, k) and k != "data_dir":
                    setattr(s, k, v)
        except (OSError, ValueError):
            pass
    return s
