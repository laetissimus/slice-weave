# Security, privacy and PHI

SliceWeave is research software for de-identified data on an encrypted workstation. It is
**not** a medical device, it is **not** HIPAA-certified and it has **not** had an
information-security review. Obtain IRB approval and your institution's security review of
this exact build before identifiable data are used.

## Data flow

```
DICOM / NIfTI on the workstation
   |  (browser upload to 127.0.0.1, or the command line reading a local folder)
   v
<data>/tmp/staging/<case>/   random file names only; deleted when the import ends
   |  headers parsed in memory -> stacks built -> descriptions scrubbed
   v
<data>/cases/<coded ID>/stacks/*.npz   pixel data, slice positions, slice thickness,
                                       acquisition order (rank only), optional mask
<data>/cases/<coded ID>/case.json      whitelisted technical fields (plane, pixel size,
                                       thickness, TE/TR, field strength, scrubbed description)
   v
local processing (numpy / scipy)
   v
<data>/cases/<coded ID>/output/        srr.nii.gz, average.nii.gz, report.json, slices.csv
   v
viewer in the local browser; exports downloaded from the local server
```

Not stored: patient name, ID, birth date, sex, accession number, study/series UIDs, dates,
times (acquisition times are reduced to a rank order and discarded), institution, physicians,
device serial numbers, or any original file. Series descriptions are kept only after removing
the patient's own name/ID tokens, dates and runs of 5 or more digits. The case ID must be a coded
study ID: IDs with 6 or more consecutive digits are refused, and an import is refused if the
case ID contains the patient name or ID found in the images.

### Exports

* **NIfTI (SRR / average), research bundle, slice table**: no patient attributes.
* **DICOM (de-identified)**: PatientID = coded case ID, empty PatientName, a new study UID,
  `ImageType DERIVED\SECONDARY`, "RESEARCH" in the series description and image comments.
* **DICOM (to original study)**: the user chooses one image file of the original study at export
  time. Its patient and study attributes (name, ID, birth date, sex, age, study UID, date, time,
  accession number, study ID and description, referring physician, frame of reference,
  institution, manufacturer) are read **in memory**, copied into the exported series only, and
  never written to the case folder. The zip is identifiable data: handle it as PHI. The series
  still carries "RESEARCH" in its description and comments.

## Network: offline by design and enforced

SliceWeave opens one listening socket, on 127.0.0.1, for the viewer. It makes no outbound
connections and no DNS lookups. This is enforced in the process (`sliceweave/netguard.py`): a
Python audit hook, installed before the server starts, refuses socket connect, bind, sendto and
DNS resolution for anything other than a loopback address, and the program's own browser launch
is limited to the local login link. The guard
is tested (`tests/test_network.py`), the program code is scanned for network, subprocess and URL
use (`tools/audit_network.py`), and the Windows build workflow runs the executables with all
outbound traffic blocked by Windows Firewall while tracing every network event
(`SLICEWEAVE_NETTRACE`).

Limits: programs started by SliceWeave run in their own processes and are not covered by the
in-process guard (only the read-only disk-encryption checks: PowerShell on Windows, fdesetup on
macOS, lsblk on Linux), and native code that opened sockets without Python would not be seen;
none of the bundled libraries (numpy, scipy, pydicom, GDCM) does so in the code paths used.

The viewer page loads only files served by the local server (fonts are bundled; no CDN). The
`Content-Security-Policy: default-src 'self'` header enforces this in the browser.

## Local server

* Listens on 127.0.0.1 only; `Host` header must be 127.0.0.1 / localhost (DNS rebinding).
* Login with a one-time random token printed in the console; session cookie HttpOnly,
  SameSite=Strict.
* State-changing requests need a CSRF header and, when the browser sends one, a matching Origin.
* Security headers: CSP, X-Frame-Options DENY, no-referrer, nosniff, no-store.

## Storage

* Data folder `~/SliceWeave-Data` (Windows: `%USERPROFILE%\SliceWeave-Data`), permissions
  restricted to the user where the OS allows.
* `doctor` and the viewer report whether the disk is encrypted (BitLocker on Windows, FileVault
  on macOS, dm-crypt on Linux). SliceWeave does not encrypt files itself: use full-disk
  encryption.
* Cases older than the retention period (default 30 days) are flagged in the viewer and can be
  purged; `sliceweave wipe` deletes cases.

## Reporting a problem

Open a private issue in the repository or contact the maintainer directly. Do not attach images.
