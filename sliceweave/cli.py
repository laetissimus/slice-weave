"""Command line: ``sliceweave <command>`` (the Windows executable accepts the same commands).

  serve      start the local viewer (default when the executable is double-clicked)
  recon      reconstruct from DICOM folders / NIfTI files without the viewer
  simulate   write simulated stacks (NIfTI) of a phantom, with the ground truth
  selftest   end-to-end test on synthetic phantoms
  doctor     versions, optional decoders, data folder, disk encryption, network guard
  wipe       delete cases (all, one, or older than N days)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import netguard
from .config import NAME, RESEARCH_BANNER, VERSION, load


def main(argv=None):
    netguard.install()
    ap = argparse.ArgumentParser(prog="sliceweave", description=f"{NAME} {VERSION}: 3-plane 2-D stacks (single-shot T2 or FLAIR) -> "
                                 f"isotropic volume (super-resolution reconstruction). {RESEARCH_BANNER}.")
    ap.add_argument("--version", action="version", version=f"{NAME} {VERSION}")
    ap.add_argument("--data", help="data folder (default ~/SliceWeave-Data or SLICEWEAVE_DATA)")
    sub = ap.add_subparsers(dest="cmd")
    s = sub.add_parser("serve", help="start the local viewer")
    s.add_argument("--port", type=int)
    s.add_argument("--no-browser", action="store_true")
    r = sub.add_parser("recon", help="reconstruct without the viewer")
    r.add_argument("inputs", nargs="+", help="DICOM folders/zips and/or NIfTI stacks (masks: files named *mask*)")
    r.add_argument("-o", "--out", required=True, help="output NIfTI (.nii.gz)")
    r.add_argument("--mode", choices=["adult", "fetal", "flair"], default="adult",
                   help="adult: single-shot T2 (adult/paediatric); fetal; flair: 2-D multi-shot FLAIR (motion "
                        "between series only)")
    r.add_argument("--resolution", type=float)
    r.add_argument("--reference", type=int, help="reference stack number (1 = first); default: least motion")
    r.add_argument("--motion", choices=["none", "stack", "slice"])
    r.add_argument("--cycles", type=int)
    r.add_argument("--alpha", type=float, help="regularisation weight (default: automatic from the noise)")
    r.add_argument("--roi-center", help="x,y,z (mm, RAS) of the region centre (fetal: the brain)")
    r.add_argument("--roi-size", help="full sizes a,b,c (mm) of the region along the grid axes")
    r.add_argument("--dicom-out", help="also write a de-identified DICOM series to this folder")
    r.add_argument("--report", help="write the report (JSON)")
    m = sub.add_parser("simulate", help="write simulated stacks of a phantom")
    m.add_argument("--kind", choices=["adult", "fetal", "flair"], default="adult")
    m.add_argument("-o", "--out", required=True)
    t = sub.add_parser("selftest", help="end-to-end test on synthetic phantoms")
    t.add_argument("--full", action="store_true", help="also the fetal (slice motion) test; slower")
    t.add_argument("--out", help="write the results (JSON)")
    sub.add_parser("doctor", help="check the installation")
    di = sub.add_parser("dicom-info", help="decode DICOM files and print a checksum of their pixels")
    di.add_argument("files", nargs="+")
    w = sub.add_parser("wipe", help="delete cases")
    w.add_argument("case", nargs="?")
    w.add_argument("--older-than", type=int)
    w.add_argument("--all", action="store_true")
    a = ap.parse_args(argv)
    settings = load(a.data)
    cmd = a.cmd or "serve"
    if cmd == "serve":
        from .server import serve
        serve(settings, getattr(a, "port", None), open_browser=not getattr(a, "no_browser", False))
    elif cmd == "recon":
        _recon(a)
    elif cmd == "simulate":
        _simulate(a)
    elif cmd == "selftest":
        from . import selftest
        res = selftest.run(full=a.full, out=Path(a.out) if a.out else None)
        print(f"\nself-test {'PASSED' if res['passed'] else 'FAILED'} in {res['seconds']} s")
        sys.exit(0 if res["passed"] else 1)
    elif cmd == "doctor":
        _doctor(settings)
    elif cmd == "dicom-info":
        import hashlib

        import numpy as np

        from . import dicom
        for f in a.files:
            ds = dicom.read(f)
            px = np.ascontiguousarray(dicom.pixel_array(ds, f, rescale=False))
            print(f"{Path(f).name}: {px.shape} {px.dtype} pixels sha256 {hashlib.sha256(px.tobytes()).hexdigest()}")
    elif cmd == "wipe":
        from .store import CaseStore
        st = CaseStore(settings)
        if a.all:
            for c in st.list():
                st.delete(c["case_id"])
            print("all cases deleted")
        elif a.older_than is not None:
            print("deleted:", st.purge_older_than(a.older_than))
        elif a.case:
            st.delete(a.case)
            print("deleted", a.case)
        else:
            print("nothing to do: give a case ID, --older-than N or --all")


def _recon(a):
    import numpy as np

    from . import series
    from .dicom_write import new_uid, write_series
    from .ingest import suggest_reference
    from .nifti import read_nifti, write_nifti
    from .recon import FLAIR_SLICE_MSG, MODE_LABELS, Params, Reconstruction, default_roi
    from .stacks import Stack
    paths = [Path(p) for p in a.inputs]
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        found = series.scan(paths, tmpdir=Path(td), log=print)
        stacks, masks = [], []
        for s in sorted(found, key=lambda s: (s.series_number or 0, s.uid)):
            if s.reject:
                print(f"  skipped: {s.reject}")
                continue
            if s.kind == "nifti" and series.is_mask_name(s.description):
                masks.append(s)
                continue
            st = series.build_stack(s, "stack")
            st.name = f"{st.plane} {len(stacks) + 1}"
            side = st.meta.pop("sidecar_text", "")
            st.meta["flair"] = series.looks_like_flair(s, text=side, ti_ms=s.ti or st.meta.get("ti_ms"),
                                                       te_ms=s.te or st.meta.get("te_ms"))
            st.meta["single_shot"] = not st.meta["flair"] and (
                series.looks_like_haste(s) or any(h in side.lower() for h in series.HASTE_HINTS))
            st.meta["desc"] = s.description
            if st.nz >= 5:
                stacks.append(st)
        for mk in masks:
            v = read_nifti(mk.nifti_path)
            probe = Stack.from_volume("m", v.data, v.affine)
            for st in stacks:
                if st.shape == probe.shape and np.allclose(st.slice_aff[:, :3, :], probe.slice_aff[:, :3, :], atol=0.05):
                    st.mask = probe.data > 0.5
                    print(f"  mask attached to {st.name}")
    if not stacks:
        sys.exit("no usable stacks")
    # one contrast per reconstruction: FLAIR stacks in flair mode, single-shot T2 otherwise (when recognised)
    key = "flair" if a.mode == "flair" else "single_shot"
    if any(st.meta.get(key) for st in stacks):
        for st in stacks:
            if not st.meta.get(key):
                print(f"  not used in {a.mode} mode: {st.name} ({st.meta.get('desc') or 'series'})")
        stacks = [st for st in stacks if st.meta.get(key)]
        for i, st in enumerate(stacks):
            st.name = f"{st.plane} {i + 1}"
    elif a.mode == "flair":
        print("NOTE: no series looks like FLAIR (name or inversion time): all stacks are used")
    if a.mode == "flair" and a.motion == "slice":
        sys.exit(FLAIR_SLICE_MSG)
    entries = [{"id": str(i), "include": True, "plane": s.plane} for i, s in enumerate(stacks)]
    ref = (a.reference - 1) if a.reference else int(suggest_reference(stacks, entries, a.mode))
    p = Params.preset(a.mode, resolution=a.resolution, motion=a.motion, cycles=a.cycles, reference=ref)
    p.alpha = a.alpha
    roi = default_roi(stacks, ref, a.mode)
    if a.roi_center:
        roi.center = np.array([float(v) for v in a.roi_center.split(",")])
    if a.roi_size:
        roi.half = np.array([float(v) for v in a.roi_size.split(",")]) / 2
    if a.mode == "fetal" and not a.roi_center and stacks[ref].mask is None:
        print("NOTE: fetal mode without --roi-center or a brain mask: the region is the centre of the reference stack")
    print(f"reference stack: {stacks[ref].name}; {len(stacks)} stacks")
    res = Reconstruction(stacks, roi, p, log=print).run()
    write_nifti(a.out, res.volume, res.grid.affine, f"SliceWeave {VERSION} SRR {MODE_LABELS[a.mode]} RESEARCH")
    print(f"wrote {a.out}  ({res.grid.shape}, {res.grid.spacing:.2f} mm)")
    if a.dicom_out:
        v = res.volume
        hi = float(np.percentile(v[v > 0], 99.9)) if (v > 0).any() else 1.0
        write_series(v, res.grid.affine, Path(a.dicom_out),
                     {"SeriesDescription": f"SliceWeave SRR {MODE_LABELS[a.mode]} RESEARCH", "PatientID": "SLICEWEAVE",
                      "PatientName": "", "StudyInstanceUID": new_uid(), "SeriesNumber": 9100,
                      "ImageType": ["DERIVED", "SECONDARY", "SRR"]}, slope=max(hi / 4000, 1e-6))
        print(f"wrote DICOM series to {a.dicom_out}")
    if a.report:
        Path(a.report).write_text(json.dumps(res.report, indent=1), encoding="utf-8")


def _simulate(a):
    from . import simulate
    from .nifti import write_nifti
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    truth, aff, info = simulate.phantom_head(a.kind, spacing=0.7 if a.kind == "fetal" else 1.0)
    ax = info["brain_half_axes_mm"]
    kw = {"static": truth, "fetal_radius": tuple(v * 1.3 for v in ax)} if a.kind == "fetal" else {}
    stacks, _ = simulate.simulate_stacks(truth, aff, a.kind, log=print, **kw)
    write_nifti(out / "truth.nii.gz", truth, aff, "SliceWeave phantom SIMULATED")
    for n, st in enumerate(stacks):
        A = st.slice_aff[0].copy()
        A[:3, 2] = st.normal * st.spacing
        write_nifti(out / f"stack{n + 1}_{st.plane}.nii.gz", st.data, A, "SliceWeave SIMULATED stack")
        side = {"SliceThickness": st.thickness,
                "SeriesDescription": f"SIM {st.plane} " + ("2D FLAIR" if a.kind == "flair" else "HASTE")}
        if st.meta.get("acquired_matrix"):   # dcm2niix-style description of the zero-filled acquisition
            side.update({"PhaseEncodingDirection": "i", "AcquisitionMatrixPE": st.meta["acquired_matrix"][0],
                         "ReconMatrixPE": st.meta["recon_matrix"][0], "BaseResolution": st.meta["acquired_matrix"][1],
                         "InversionTime": 2.5, "EchoTime": 0.086, "RepetitionTime": 9.0})
        (out / f"stack{n + 1}_{st.plane}.json").write_text(json.dumps(side), encoding="utf-8")
    print(f"wrote {len(stacks)} simulated stacks and truth.nii.gz to {out}")


def _doctor(settings):
    import platform

    import numpy
    import scipy

    from . import hooks, security
    print(f"{NAME} {VERSION}   Python {platform.python_version()} on {platform.platform()}")
    h = hooks.status()
    print(f"numpy {numpy.__version__}, scipy {scipy.__version__}; JPEG-family DICOM: pydicom "
          f"{h['pydicom'] or 'not installed'}, GDCM {h['gdcm'] or 'not installed'}")
    settings.ensure()
    print(f"data folder: {settings.root}  ({security.folder_permissions(settings.root)['status']})")
    enc = security.disk_encryption_status(settings.root)
    print(f"disk encryption: {enc['status']} ({enc['detail']})")
    st = netguard.status()
    print(f"Network guard: {'ACTIVE' if st['active'] else 'OFF'} (only 127.0.0.1 can be contacted)")
