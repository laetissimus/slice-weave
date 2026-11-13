"""De-identification at ingest.

Policy (see SECURITY.md):
* Original DICOM files are never copied into the data folder. Uploads are
  read from a temporary folder that is deleted when ingest finishes.
* Images are stored as freshly written NIfTI files numbered by stack; no header
  field of any input file is copied.
* Of the DICOM attributes, only a whitelist of technical acquisition
  parameters is kept. Free-text fields that are kept (series description) are
  scrubbed of the patient's own name/ID tokens, long digit runs and dates.
* Dates and times are not stored. Series are ordered by series number.
* The case is known only by a coded study ID entered by the user. The tool
  refuses IDs that look like an MRN or that contain an identifier found in
  the data. No hash or "fingerprint" of patient identifiers is stored.
* Series that may carry burned-in identifiers (secondary capture, colour,
  BurnedInAnnotation=YES) are refused. Overlay planes are never imported.
  Pixel text inside ordinary greyscale MR images is NOT detected (no OCR).
"""
from __future__ import annotations

import datetime as _dt
import re

from .series import SeriesInfo

CASE_ID_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_-]{2,39}$")
MRN_LIKE = re.compile(r"\d{6,}")
_DATE_PATTERNS = [
    re.compile(r"\b(19|20)\d{2}[-/._]?(0[1-9]|1[0-2])[-/._]?(0[1-9]|[12]\d|3[01])\b"),
    re.compile(r"\b(0?[1-9]|1[0-2])[-/.](0?[1-9]|[12]\d|3[01])[-/.](19|20)?\d{2}\b"),
]


class DeidError(ValueError):
    pass


def _phi_tokens(phi: dict) -> set[str]:
    toks = set()
    for k, v in phi.items():
        if not v:
            continue
        if k in ("PatientName", "OtherPatientNames", "ReferringPhysicianName", "PerformingPhysicianName",
                 "OperatorsName"):
            for t in re.split(r"[\^\s,\\]+", v):
                if len(t) >= 2:
                    toks.add(t.lower())
        else:
            for t in re.split(r"[\\]+", v):
                t = t.strip()
                if len(t) >= 3:
                    toks.add(t.lower())
    return toks


def scrub(text: str, phi_tokens: set[str]) -> str:
    if not text:
        return ""
    out = text
    for t in sorted(phi_tokens, key=len, reverse=True):
        out = re.sub(re.escape(t), "[x]", out, flags=re.I)
    for pat in _DATE_PATTERNS:
        out = pat.sub("[date]", out)
    out = re.sub(r"\d{5,}", "[n]", out)
    return out[:64]


def validate_case_id(case_id: str, series: list[SeriesInfo] | None = None, allow_digits: bool = False) -> list[str]:
    problems = []
    if not CASE_ID_RE.match(case_id or ""):
        problems.append("Case ID must start with a letter and use 3-40 letters, digits, '-' or '_'.")
    if not allow_digits and MRN_LIKE.search(case_id or ""):
        problems.append("Case ID contains 6 or more consecutive digits and may be an MRN. Use a coded study ID.")
    if series:
        low = (case_id or "").lower()
        for s in series:
            for k, v in s.phi.items():
                if k in ("PatientID", "AccessionNumber", "OtherPatientIDs", "StudyID") and v and len(v) >= 3 \
                        and v.lower() in low:
                    problems.append(f"Case ID contains the {k} found in the images.")
                if k == "PatientName" and v:
                    for t in re.split(r"[\^\s]+", v):
                        if len(t) >= 3 and t.lower() in low:
                            problems.append("Case ID contains part of the patient name found in the images.")
                            break
        problems = sorted(set(problems))
    return problems


def validate_site(site: str) -> list[str]:
    if not site:
        return []
    if len(site) > 32 or not re.match(r"^[A-Za-z0-9 _.-]+$", site):
        return ["Site label: up to 32 letters, digits, spaces, '.', '-' or '_'."]
    if MRN_LIKE.search(site):
        return ["Site label must not contain long digit runs."]
    return []


