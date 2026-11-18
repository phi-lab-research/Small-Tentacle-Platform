import sys
from pathlib import Path

import equinox as eqx
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "model"))

from small_tentacle_platform_soromox import SoftArm

# (name, d_x along backbone from model base [m], [d_y, d_z] offset in local frame [m])
# from Small-Tentacle-Platform/data/soft_arm_parameters.yaml
MARKERS = [
    ("base", 0.013, (0.0, 0.030)),
    ("1", 0.1291, (0.0, -0.025)),
    ("2", 0.2195, (0.0, -0.023)),
    ("3", 0.280, (0.0, -0.021)),
    ("tip", 0.370, (0.0, 0.0)),
]

model = SoftArm().model
L_model = float(jnp.sum(model.segment_length))

S = jnp.array([min(d_x, L_model) for _, d_x, _ in MARKERS])
OFFSETS = jnp.array([[0.0, d_y, d_z, 1.0] for _, _, (d_y, d_z) in MARKERS])
TENDON_REST = model.actuators[0].path_lengths(model, jnp.zeros(model.num_internal_dofs))


def marker_positions(q):
    """FK positions of all markers for strain coefficients q, shape (len(MARKERS), 3)."""
    g = jax.vmap(lambda s: model.forward_kinematics(q, s))(S)
    return jnp.einsum("nij,nj->ni", g, OFFSETS)[:, :3]


def tendon_lengths(q):
    """Length change of the 4 tendons relative to rest [m], shape (4,)."""
    return model.actuators[0].path_lengths(model, q) - TENDON_REST


def residuals(q, p_meas, l_meas, w_tendon=1.0):
    """Marker position errors (15,) followed by sqrt(w_tendon) * tendon length errors (4,).

    p_meas: (len(MARKERS), 3) marker positions [m]; l_meas: (4,) tendon length changes [m].
    """
    return jnp.concatenate(
        [
            (marker_positions(q) - p_meas).ravel(),
            jnp.sqrt(w_tendon) * (tendon_lengths(q) - l_meas),
        ]
    )


def cost(q, p_meas, l_meas, w_tendon=1.0):
    return jnp.sum(residuals(q, p_meas, l_meas, w_tendon) ** 2)


cost_jit = eqx.filter_jit(cost)
grad_jit = eqx.filter_jit(jax.grad(cost))


@eqx.filter_jit
def solve(q0, p_meas, l_meas, w_tendon=1.0, iters=10, lam0=1e-3):
    """Levenberg-Marquardt: solve (J^T J + lam I) dq = -J^T r; accept dq and shrink lam
    if the cost drops, otherwise reject dq and grow lam (more gradient-descent-like)."""

    def step(carry, _):
        q, lam = carry
        r = residuals(q, p_meas, l_meas, w_tendon)
        J = jax.jacfwd(residuals)(q, p_meas, l_meas, w_tendon)
        dq = jnp.linalg.solve(J.T @ J + lam * jnp.eye(q.size), -J.T @ r)
        better = cost(q + dq, p_meas, l_meas, w_tendon) < r @ r
        return (jnp.where(better, q + dq, q), jnp.where(better, lam / 3, lam * 3)), None

    return jax.lax.scan(step, (q0, lam0), None, length=iters)[0][0]


def qd_finite_difference(q, t, N=5):
    """Causal least-squares slope of q over the last N samples, shape of q; first N-1 rows are zero."""
    idx = jnp.arange(N - 1, len(t))[:, None] - jnp.arange(N)[::-1]
    tc = t[idx] - jnp.mean(t[idx], axis=1, keepdims=True)
    qd = jnp.einsum("kw,kwd->kd", tc, q[idx]) / jnp.sum(tc**2, axis=1, keepdims=True)
    return jnp.concatenate([jnp.zeros((N - 1, q.shape[1])), qd])


def kalman_step(x, P, z, dt, sigma_a, sigma_z):
    """One predict+update of the per-DoF constant-velocity filter; x (n, 2) = [q, qdot], P (2, 2) shared by all DoFs."""
    H = jnp.array([1.0, 0.0])
    F = jnp.array([[1.0, dt], [0.0, 1.0]])
    Qm = sigma_a**2 * jnp.array([[dt**4 / 4, dt**3 / 2], [dt**3 / 2, dt**2]])
    x = x @ F.T
    P = F @ P @ F.T + Qm
    K = P @ H / (H @ P @ H + sigma_z**2)
    x = x + (z - x @ H)[:, None] * K
    P = P - jnp.outer(K, H @ P)
    return x, P


def qd_kalman(q, t, sigma_a=1.0, sigma_z=1e-2):
    """Constant-velocity Kalman filter per DoF, state [q, qdot], measurement q; returns qdot with shape of q, first row zero."""

    def step(carry, inp):
        x, P = kalman_step(*carry, *inp, sigma_a, sigma_z)
        return (x, P), x[:, 1]

    x0 = jnp.stack([q[0], jnp.zeros_like(q[0])], axis=1)
    P0 = jnp.diag(jnp.array([sigma_z**2, 1e2]))
    qd = jax.lax.scan(step, (x0, P0), (q[1:], jnp.diff(t)))[1]
    return jnp.concatenate([jnp.zeros((1, q.shape[1])), qd])


if __name__ == "__main__":
    robot = SoftArm()
    n = model.num_internal_dofs
    # tendon forces [N], positive = pulling
    U_CASES = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 3.0],
        [2.0, 0.0, 0.0, 2.0],
        [0.0, 5.0, 0.0, 0.0],
        [1.0, 2.0, 3.0, 4.0],
        [5.0, 0.0, 0.0, 5.0],
    ]
    for u in U_CASES:
        q_true = robot.equilibrium(jnp.asarray(u))
        p_meas = marker_positions(q_true)
        l_meas = tendon_lengths(q_true)
        q_hat = solve(jnp.zeros(n), p_meas, l_meas)
        print(
            f"u={u}  cost={float(cost(q_hat, p_meas, l_meas)):.2e}"
            f"  marker err [mm]={1e3 * float(jnp.abs(marker_positions(q_hat) - p_meas).max()):.1f}"
            f"  tendon err [mm]={1e3 * float(jnp.abs(tendon_lengths(q_hat) - l_meas).max()):.1f}"
            f"  |q_hat-q_true|={float(jnp.abs(q_hat - q_true).max()):.3f}"
        )
