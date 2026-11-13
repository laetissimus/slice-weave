"""Case IDs and de-identification."""
from types import SimpleNamespace

from sliceweave.deid import scrub, validate_case_id


def test_case_id_rules():
    assert validate_case_id("FET-012") == []
    assert validate_case_id("MRN1234567")          # long digit run: may be an MRN
    assert validate_case_id("")
    s = SimpleNamespace(phi={"PatientName": "DOE^JANE", "PatientID": "X99Y"})
    assert validate_case_id("DOE-01", [s])          # contains a name found in the images
    assert validate_case_id("CASE-X99Y", [s])       # contains the patient ID
    assert validate_case_id("FET-013", [s]) == []


def test_scrub_removes_names_dates_and_numbers():
    out = scrub("HASTE DOE JANE 2026-03-04 acc 1234567", {"DOE", "JANE"})
    assert "DOE" not in out and "JANE" not in out and "2026" not in out and "1234567" not in out
