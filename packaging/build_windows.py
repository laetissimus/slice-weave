"""Build the standalone Windows executables with PyInstaller.

Produces in ./dist:
  SliceWeave.exe               one file; double-click to run (unpacks itself to a temp folder)
  SliceWeave-portable.zip      a folder version that starts faster (unzip, run SliceWeave.exe)
  SHA256SUMS.txt, build-report.json

Run on Windows 10/11 x64 with Python 3.13:
  python -m pip install -r requirements.txt "pyinstaller>=6.10,<7"
  python packaging/build_windows.py
(The GitHub workflow .github/workflows/windows-build.yml does exactly this, then tests the .exe.)
The same script builds a Linux or macOS binary when run there.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import platform
import shutil
import sys
import zipfile
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DIST = ROOT / "dist"
BUILD = ROOT / "build"
NAME = "SliceWeave"
EXCLUDES = ["tkinter", "_tkinter", "pytest", "IPython", "notebook", "pandas", "skimage", "torch", "matplotlib",
            "PyQt5", "PyQt6", "PySide2", "PySide6", "wx", "gi", "sphinx", "docutils", "PIL"]
sys.path.insert(0, str(ROOT))
from sliceweave.config import VERSION  # noqa: E402

README = f"""SliceWeave {VERSION} - standalone build
Super-resolution reconstruction of 3-plane 2-D HASTE stacks into one isotropic volume
(adult, paediatric and fetal brain MRI).
RESEARCH USE ONLY - NOT FOR CLINICAL DIAGNOSIS

Nothing to install. Double-click SliceWeave.exe:
  * a console window opens and shows a one-time login link;
  * your browser opens SliceWeave at http://127.0.0.1:8766 (this computer only);
  * keep the console window open while you work; close it to stop SliceWeave.

Works offline. Nothing is sent over the network: only 127.0.0.1 can be contacted.
Data are stored in %USERPROFILE%\\SliceWeave-Data (change with --data or SLICEWEAVE_DATA).
Use a computer with full-disk encryption (BitLocker); the console shows its status.

From a terminal the same commands as `sliceweave` are available, for example:
  SliceWeave.exe recon --mode fetal -o out.nii.gz <DICOM folders or NIfTI stacks>
  SliceWeave.exe doctor      checks the installation and disk encryption
  SliceWeave.exe selftest    end-to-end test on synthetic phantoms (no patient data)
  SliceWeave.exe -h          all commands

The executable is not code-signed. Windows SmartScreen may say "Windows protected your PC":
choose "More info" then "Run anyway", or ask IT to allow or sign it. Managed computers with
application control (AppLocker / WDAC) need IT approval before it can run.
"""


def _args(mode: str) -> list[str]:
    sep = os.pathsep
    args = ["--noconfirm", "--clean", "--console", "--name", NAME, f"--{mode}",
            "--paths", str(ROOT),
            "--add-data", f"{ROOT / 'sliceweave' / 'web'}{sep}sliceweave/web",
            "--collect-submodules", "uvicorn", "--collect-submodules", "sliceweave",
            "--hidden-import", "python_multipart",
            "--distpath", str(DIST if mode == "onefile" else DIST / "portable"),
            "--workpath", str(BUILD / mode), "--specpath", str(BUILD / mode)]
    for m in EXCLUDES:
        args += ["--exclude-module", m]
    # optional JPEG-family DICOM decoding: pydicom (MIT) + GDCM (BSD), bundled when installed
    if importlib.util.find_spec("pydicom"):
        args += ["--collect-submodules", "pydicom", "--collect-data", "pydicom", "--copy-metadata", "pydicom"]
    if importlib.util.find_spec("gdcm"):
        args += ["--collect-all", "gdcm"]
    if os.name == "nt":
        args += ["--version-file", str(ROOT / "packaging" / "version_info.txt")]
    return args + [str(ROOT / "packaging" / "launcher.py")]


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    import PyInstaller.__main__
    modes = ("onefile",) if "--onefile-only" in sys.argv else ("onefile", "onedir")
    shutil.rmtree(DIST, ignore_errors=True)
    shutil.rmtree(BUILD, ignore_errors=True)
    for mode in modes:
        print(f"== PyInstaller {mode}", flush=True)
        PyInstaller.__main__.run(_args(mode))
    exe = DIST / (NAME + (".exe" if os.name == "nt" else ""))
    (DIST / "README-WINDOWS.txt").write_text(README, encoding="utf-8")
    outputs = [exe]
    if "onedir" in modes:
        folder = DIST / "portable" / NAME
        (folder / "README-WINDOWS.txt").write_text(README, encoding="utf-8")
        zpath = DIST / f"{NAME}-portable.zip"
        with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
            for f in sorted(folder.rglob("*")):
                if f.is_file():
                    z.write(f, str(Path(NAME) / f.relative_to(folder)))
        outputs.append(zpath)
    sums = {p.name: _sha256(p) for p in outputs}
    (DIST / "SHA256SUMS.txt").write_text("".join(f"{v}  {k}\n" for k, v in sums.items()), encoding="utf-8")
    pkgs = {}
    for name in ("numpy", "scipy", "starlette", "uvicorn", "python-multipart", "pyinstaller", "pydicom",
                 "python-gdcm"):
        try:
            pkgs[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            pkgs[name] = None
    report = {"sliceweave_version": VERSION, "python": platform.python_version(), "platform": platform.platform(),
              "machine": platform.machine(), "packages": pkgs, "sha256": sums,
              "sizes_mb": {p.name: round(p.stat().st_size / 1e6, 1) for p in outputs}}
    (DIST / "build-report.json").write_text(json.dumps(report, indent=1), encoding="utf-8")
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
