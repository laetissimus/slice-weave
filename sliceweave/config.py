"""Settings. Stored as JSON in the data folder; environment variables override.

Nothing here is secret and nothing identifies a patient.
"""
from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path

VERSION = "1.1.0"
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
