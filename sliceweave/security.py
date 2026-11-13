"""Best-effort environment checks reported at start-up and by ``sliceweave doctor``.

These checks inform; they do not replace an institutional information-security review.
"""
from __future__ import annotations

import os
import platform
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path


def _run(cmd, timeout=10) -> str:
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return (r.stdout or "") + (r.stderr or "")
    except (OSError, subprocess.TimeoutExpired):
        return ""


@lru_cache(maxsize=8)
def disk_encryption_status(path) -> dict:
    sysname = platform.system()
    p = str(Path(path).expanduser())
    if sysname == "Darwin":
        out = _run(["fdesetup", "status"])
        if "FileVault is On" in out:
            return {"status": "encrypted", "detail": "FileVault is on"}
        if "FileVault is Off" in out:
            return {"status": "NOT ENCRYPTED", "detail": "FileVault is off"}
        return {"status": "unknown", "detail": "could not run fdesetup"}
    if sysname == "Windows":
        drive = os.path.splitdrive(p)[0] or "C:"
        # Shell property System.Volume.BitLockerProtection is readable without admin rights:
        # 1 = on, 2 = off, 3 = encrypting, 4 = decrypting, 5 = suspended, 6 = locked
        ps = ("(New-Object -ComObject Shell.Application).NameSpace('" + drive + "\\')"
              ".Self.ExtendedProperty('System.Volume.BitLockerProtection')")
        val = _run(["powershell", "-NoProfile", "-NonInteractive", "-Command", ps]).strip()
        states = {"1": ("encrypted", "BitLocker on"), "3": ("encrypted", "BitLocker encryption in progress"),
                  "6": ("encrypted", "BitLocker on (locked)"), "2": ("NOT ENCRYPTED", "BitLocker off"),
                  "4": ("NOT ENCRYPTED", "BitLocker decryption in progress"),
                  "5": ("NOT ENCRYPTED", "BitLocker protection suspended")}
        if val in states:
            st, det = states[val]
            return {"status": st, "detail": f"{det} for {drive}"}
        out = _run(["manage-bde", "-status", drive])
        if "Protection On" in out:
            return {"status": "encrypted", "detail": f"BitLocker protection on for {drive}"}
        if "Protection Off" in out:
            return {"status": "NOT ENCRYPTED", "detail": f"BitLocker protection off for {drive}"}
        return {"status": "unknown", "detail": "manage-bde needs administrator rights; ask IT to confirm BitLocker"}
    if sysname == "Linux" and shutil.which("lsblk") and shutil.which("df"):
        src = _run(["df", "--output=source", p]).splitlines()
        dev = src[-1].strip() if len(src) >= 2 else ""
        if dev.startswith("/dev/"):
            out = _run(["lsblk", "-s", "-n", "-o", "TYPE", dev])
            if "crypt" in out.split():
                return {"status": "encrypted", "detail": f"{dev} is on a dm-crypt/LUKS device"}
            return {"status": "NOT ENCRYPTED", "detail": f"no dm-crypt layer found under {dev}"}
        return {"status": "unknown", "detail": f"data folder is on {dev or 'an unknown device'}"}
    return {"status": "unknown", "detail": f"no check implemented for {sysname}"}


def folder_permissions(path) -> dict:
    p = Path(path).expanduser()
    if os.name == "nt" or not p.exists():
        return {"status": "unknown"}
    mode = p.stat().st_mode & 0o777
    return {"status": "private" if mode & 0o077 == 0 else "readable by others", "mode": oct(mode)}
