"""Reconstruction engine: geometry, the acquisition operator, registration and super-resolution."""
import numpy as np
import pytest

from sliceweave import metrics, operators as ops, simulate
from sliceweave.geometry import (about_center, apply, make_grid, matrix_to_rotvec, nearest_ras_axes, rigid_params,
                                 rotvec_to_matrix)
from sliceweave.recon import Params, Reconstruction, auto_alpha, default_roi
from sliceweave.register import Sampler, register
from sliceweave.stacks import Stack, packets


def test_rotations_round_trip():
    rng = np.random.default_rng(0)
    for _ in range(20):
        w = rng.normal(0, 1, 3)
        w *= rng.uniform(0, 3.0) / np.linalg.norm(w)
        R = rotvec_to_matrix(w)
        assert np.allclose(R @ R.T, np.eye(3), atol=1e-10) and np.isclose(np.linalg.det(R), 1)
        assert np.allclose(rotvec_to_matrix(matrix_to_rotvec(R)), R, atol=1e-8)


def test_nearest_ras_axes_is_proper_rotation():
    R = rotvec_to_matrix([0.2, -0.1, 0.3]) @ np.eye(3)[:, [2, 0, 1]]
    A = nearest_ras_axes(R)
    assert np.isclose(np.linalg.det(A), 1) and np.all(np.diag(A) > 0.5)


@pytest.fixture(scope="module")
def adult():
    truth, aff, info = simulate.phantom_head("adult", spacing=1.2)
    stacks, Ts = simulate.simulate_stacks(truth, aff, "adult", motion=False, noise=0.01)
    return truth, aff, info, stacks, Ts


def test_operator_adjoint(adult):
    _, _, _, stacks, Ts = adult
    g = make_grid(nearest_ras_axes(stacks[0].axes), [0, 0, 0], [40, 40, 40], 2.0, margin_mm=4)
    masks = [np.ones(s.shape, bool) for s in stacks]
    op = ops.Operator.build(g, stacks, Ts, masks)
    rng = np.random.default_rng(1)
    x = rng.random(g.shape).astype(np.float32)
    r = [rng.random(gr.n_pix).astype(np.float32) for gr in op.groups]
    lhs = sum(float(np.dot(a.astype(np.float64), b)) for a, b in zip(op.forward(x), r))
    rhs = float(np.vdot(x.astype(np.float64), op.adjoint(r)))
    assert abs(lhs - rhs) / abs(rhs) < 1e-4


def test_registration_recovers_known_motion(adult):
    truth, aff, _, stacks, _ = adult
    from sliceweave.geometry import Grid
    smp = Sampler(truth, Grid(truth.shape, aff))
    st = stacks[0]
    ii, jj, ss = np.nonzero(st.data > 0.2 * st.data.max())
    sel = np.random.default_rng(2).choice(ii.size, 20000, replace=False)
    A = st.slice_aff[ss[sel]]
    P = np.einsum("nij,nj->ni", A[:, :3, :3], np.stack([ii[sel], jj[sel], np.zeros(sel.size)], 1)) + A[:, :3, 3]
    y, _, _ = smp.sample(P, grad=False)
    D_true = about_center(rotvec_to_matrix(np.deg2rad([3.0, -2.0, 4.0])), [0, 0, 0], [2.0, -1.5, 1.0])
    # pixels displaced by the inverse motion: registration must find D_true
    r = register(smp, apply(np.linalg.inv(D_true), P), y, center=[0, 0, 0], iters=40)
    ang, disp = rigid_params(np.linalg.inv(D_true) @ r.D, [0, 0, 0])
    assert r.ok and ang < 0.2 and disp < 0.2


def test_srr_beats_average_and_recovers_series_motion(adult):
    truth, aff, info, stacks, Ts = adult
    p = Params.preset("adult", resolution=2.0)
    p.final_iterations, p.cycle_iterations = 15, 8
    res = Reconstruction(stacks, default_roi(stacks, 0, "adult"), p).run()
    g = res.grid
    from scipy import ndimage as ndi
    W = g.vox_to_world(np.indices(g.shape).reshape(3, -1).T)
    tg = ndi.map_coordinates(truth, apply(np.linalg.inv(aff), W).T, order=1).reshape(g.shape)
    a, b, c = info["brain_half_axes_mm"]
    m = (((W[:, 0] / a) ** 2 + (W[:, 1] / b) ** 2 + (W[:, 2] / c) ** 2) < 1).reshape(g.shape)
    srr, avg = metrics.compare(res.volume, tg, m), metrics.compare(res.average, tg, m)
    assert srr["psnr_db"] > avg["psnr_db"]
    for k in (1, 2):
        ang, disp = rigid_params(np.linalg.inv(Ts[k][0]) @ res.transforms[k][0], [0, 0, 0])
        assert ang < 0.5 and disp < 0.5


def test_stack_from_volume_puts_thick_axis_last():
    d = np.zeros((30, 8, 40), np.float32)
    A = np.diag([0.8, 4.0, 0.8, 1.0])
    st = Stack.from_volume("s", d, A)
    assert st.shape == (30, 40, 8) and np.isclose(st.spacing, 4.0) and st.plane == "COR"


def test_packets_cover_every_slice_once():
    st = Stack("s", np.zeros((4, 4, 23), np.float32), np.repeat(np.eye(4)[None], 23, 0), 3.0)
    pk = packets(st, 4)
    assert sorted(s for p in pk for s in p) == list(range(23)) and len(pk) == 4
    st.meta["acq_order"] = list(range(23))[::-1]
    assert packets(st, 4)[0][0] == 22


def test_auto_alpha_grows_with_noise():
    assert auto_alpha(0.01, "huber") < auto_alpha(0.05, "huber") < auto_alpha(0.1, "huber")
    # final volumes of FLAIR (small lesions) and adult / paediatric T2 (cortex, sulci): weaker weight
    # than the registration volumes; fetal unchanged
    for mode in ("flair", "adult"):
        assert auto_alpha(0.05, "huber", mode) < 0.5 * auto_alpha(0.05, "huber")
    assert auto_alpha(0.05, "huber", "fetal") == auto_alpha(0.05, "huber")


@pytest.mark.parametrize("motion,cycles,outliers,bias,expect_reg,expect_weights", [
    ("stack", 2, True, True, 2, 2),     # two cycles: register, then intensity/outlier check, twice
    ("stack", 0, True, True, 0, 1),     # 0 cycles: no refinement, intensity/outlier check still once
    ("none", 5, False, True, 0, 1),     # no motion correction: cycles ignored, intensity correction once
    ("none", 5, False, False, 0, 0),    # nothing to do between SRRs
])
def test_cycles_schedule(adult, monkeypatch, motion, cycles, outliers, bias, expect_reg, expect_weights):
    truth, aff, info, stacks, Ts = adult
    p = Params.preset("adult", resolution=3.0, motion=motion, cycles=cycles, outlier_rejection=outliers,
                      bias_correction=bias)
    p.final_iterations, p.cycle_iterations = 2, 2
    calls = {"reg": 0, "weights": 0}
    rec = Reconstruction(stacks, default_roi(stacks, 0, "adult"), p, log=lambda m: None)
    orig_reg, orig_w = rec.register_stacks, rec.update_weights

    def reg(*a, **k):
        calls["reg"] += 1
        return orig_reg(*a, **k)

    def weights(*a, **k):
        calls["weights"] += 1
        return orig_w(*a, **k)
    monkeypatch.setattr(rec, "register_stacks", reg)
    monkeypatch.setattr(rec, "update_weights", weights)
    rec.run()
    # register_stacks also runs once before the first SRR when motion correction is on
    assert calls["reg"] == expect_reg + (1 if motion != "none" else 0)
    assert calls["weights"] == expect_weights
