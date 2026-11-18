from __future__ import annotations

from typing import Any

import equinox as eqx
import jax
import jax.numpy as jnp
import optimistix as optx
from jax import Array

from soromox.actuation import ThreadlikeActuator, ThreadlikeRouting
from soromox.actuation.threadlike import BaseThreadlikeRoutingParams
from soromox.systems import (
    GVS,
    GVSSegment,
    JointSpec,
    LinearProfile,
    LinkSpec,
    StrainBasisSpec,
    SystemState,
)

__all__ = [
    "SegmentedLinearThreadlikeRoutingParams",
    "segmented_linear_threadlike_routing",
    "segmented_linear_threadlike_routing_derivative",
    "linear_segment",
]


def linear_segment(
    R_start: Array, R_end: Array, angle_start: Array, angle_end: Array, length: float
):
    y_start, z_start = R_start * jnp.cos(angle_start), R_start * jnp.sin(angle_start)
    y_end, z_end = R_end * jnp.cos(angle_end), R_end * jnp.sin(angle_end)
    return (
        y_start,
        (y_end - y_start) / length,
        z_start,
        (z_end - z_start) / length,
    )


def _index_tuple(value: Any) -> tuple[int, ...]:
    if isinstance(value, int):
        return (int(value),)
    if isinstance(value, (tuple, list)):
        return tuple(int(index) for index in value)
    array = jnp.asarray(value)
    if array.ndim == 0:
        return (int(array.item()),)
    return tuple(int(index) for index in array.tolist())


def _float_tuple(value: Any) -> tuple[float, ...]:
    array = jnp.asarray(value)
    return tuple(float(length) for length in array.tolist())


def _segment_index_and_local_s(
    segment_lengths: tuple[float, ...], s: Array
) -> tuple[Array, Array]:
    boundaries = jnp.cumsum(jnp.asarray((0.0, *segment_lengths)))
    num_segments = len(segment_lengths)
    segment_index = jnp.clip(
        jnp.searchsorted(boundaries[1:], s, side="left"), 0, num_segments - 1
    )
    s_local = s - boundaries[segment_index]
    return segment_index, s_local


def segmented_linear_threadlike_routing(
    params: "SegmentedLinearThreadlikeRoutingParams", s: Array
) -> Array:
    segment_index, s_local = _segment_index_and_local_s(params.segment_lengths, s)
    y = params.y_intercept[segment_index] + params.y_slope[segment_index] * s_local
    z = params.z_intercept[segment_index] + params.z_slope[segment_index] * s_local
    return jnp.stack([jnp.zeros_like(y), y, z])


def segmented_linear_threadlike_routing_derivative(
    params: "SegmentedLinearThreadlikeRoutingParams", s: Array
) -> Array:
    segment_index, _ = _segment_index_and_local_s(params.segment_lengths, s)
    dy = params.y_slope[segment_index]
    dz = params.z_slope[segment_index]
    return jnp.stack([jnp.zeros_like(dy), dy, dz])


class SegmentedLinearThreadlikeRoutingParams(BaseThreadlikeRoutingParams):
    y_intercept: Array
    y_slope: Array
    z_intercept: Array
    z_slope: Array
    segment_lengths: tuple[float, ...] = eqx.field(static=True, converter=_float_tuple)
    start_segment_index: tuple[int, ...] = eqx.field(
        static=True, converter=_index_tuple
    )
    end_segment_index: tuple[int, ...] = eqx.field(static=True, converter=_index_tuple)

    def __check_init__(self) -> None:
        self.validate()

    @property
    def num_paths(self) -> int:
        return int(jnp.asarray(self.y_intercept).shape[0])

    def validate(self) -> None:
        y_intercept = jnp.asarray(self.y_intercept)
        if y_intercept.ndim != 2:
            raise ValueError(
                "SegmentedLinearThreadlikeRoutingParams fields must be "
                "two-dimensional with shape (num_paths, num_segments)."
            )
        num_segments = len(self.segment_lengths)
        expected_shape = (y_intercept.shape[0], num_segments)
        if y_intercept.shape != expected_shape:
            raise ValueError(
                f"y_intercept must have shape {expected_shape}, got {y_intercept.shape}."
            )
        for name in ("y_slope", "z_intercept", "z_slope"):
            value = jnp.asarray(getattr(self, name))
            if value.shape != expected_shape:
                raise ValueError(
                    f"{name} must have shape {expected_shape}, got {value.shape}."
                )
        num_paths = expected_shape[0]
        if len(self.start_segment_index) != num_paths:
            raise ValueError("start_segment_index must contain one entry per path.")
        if len(self.end_segment_index) != num_paths:
            raise ValueError("end_segment_index must contain one entry per path.")
        for start, end in zip(self.start_segment_index, self.end_segment_index):
            if start < 0 or start > end or end >= num_segments:
                raise ValueError(
                    "each routing span must satisfy 0 <= start <= end < "
                    f"num_segments; got ({start}, {end}) for {num_segments} segments."
                )


class SoftArm:
    def __init__(self, num_gauss_points=(8, 8, 8, 8), gravity=(9.81, 0.0, 0.0)):
        E_all = 3.02e5
        nu_all = 0.4
        rho_all = 1310.0
        eta_all = 5e3

        L0 = 0.02
        segment_lengths = jnp.array([0.025, 0.255, 0.025, 0.055])
        L_total = L0 + jnp.sum(segment_lengths)
        cum_L = L0 + jnp.concatenate([jnp.zeros(1), jnp.cumsum(segment_lengths)[:-1]])
        Rb0 = 0.016
        Rt4 = 0.0048
        Rt0123 = Rb0 + (Rt4 - Rb0) * cum_L / L_total
        Rb = Rt0123
        Rt = jnp.concatenate([Rt0123[1:], jnp.array([Rt4])])

        link1, link2, link3, link4 = [
            LinkSpec.circular(
                length=segment_lengths[i],
                radius=LinearProfile(base=Rb[i], tip=Rt[i]),
                density=rho_all,
                young_modulus=E_all,
                poisson_ratio=nu_all,
                material_damping_coefficient=eta_all,
                reference_strain=[0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            )
            for i in range(4)
        ]

        joint = JointSpec.fixed()

        basis1 = StrainBasisSpec(
            type="monomial",
            strain_selector=[0, 1, 1, 1, 0, 0],
            basis_order=[0, 0, 0, 0, 0, 0],
        )
        basis2 = StrainBasisSpec(
            type="monomial",
            strain_selector=[1, 1, 1, 1, 0, 0],
            basis_order=[1, 1, 1, 1, 0, 0],
        )
        basis3 = StrainBasisSpec(
            type="monomial",
            strain_selector=[0, 1, 1, 1, 0, 0],
            basis_order=[0, 0, 0, 0, 0, 0],
        )
        basis4 = StrainBasisSpec(
            type="monomial",
            strain_selector=[0, 1, 1, 0, 0, 0],
            basis_order=[0, 0, 0, 0, 0, 0],
        )

        clearances = jnp.array([0.004, 0.004, 0.003, 0.003])
        R_disc = Rb - clearances
        tendon_angles_rad = jnp.deg2rad(jnp.array([30.0, 60.0, 120.0, 150.0]))
        angles_crossed = tendon_angles_rad[jnp.array([0, 2, 1, 3])]

        seg1 = linear_segment(
            R_disc[0], R_disc[1], tendon_angles_rad, tendon_angles_rad, link1.length
        )
        seg2 = linear_segment(
            R_disc[1], R_disc[2], tendon_angles_rad, angles_crossed, link2.length
        )
        seg3 = linear_segment(
            R_disc[2], R_disc[3], angles_crossed, angles_crossed, link3.length
        )
        seg4 = jnp.zeros(4)

        segmented_tendon_routing = SegmentedLinearThreadlikeRoutingParams(
            y_intercept=jnp.stack([seg1[0], seg2[0], seg3[0], seg4], axis=1),
            y_slope=jnp.stack([seg1[1], seg2[1], seg3[1], seg4], axis=1),
            z_intercept=jnp.stack([seg1[2], seg2[2], seg3[2], seg4], axis=1),
            z_slope=jnp.stack([seg1[3], seg2[3], seg3[3], seg4], axis=1),
            segment_lengths=segment_lengths,
            start_segment_index=(0, 0, 0, 0),
            end_segment_index=(2, 2, 2, 2),
        )
        active_tendon_routing = ThreadlikeRouting(
            params=segmented_tendon_routing,
            offset_fn=segmented_linear_threadlike_routing,
            derivative_fn=segmented_linear_threadlike_routing_derivative,
        )

        p0 = jnp.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0])  # [qx,qy,qz,qw,x,y,z]

        segments = [
            GVSSegment(
                link=link1, joint=joint, basis=basis1, num_gauss_points=num_gauss_points[0]
            ),
            GVSSegment(
                link=link2, joint=joint, basis=basis2, num_gauss_points=num_gauss_points[1]
            ),
            GVSSegment(
                link=link3, joint=joint, basis=basis3, num_gauss_points=num_gauss_points[2]
            ),
            GVSSegment(
                link=link4, joint=joint, basis=basis4, num_gauss_points=num_gauss_points[3]
            ),
        ]

        self.model = GVS.from_segments(
            segments,
            gravity=jnp.asarray(gravity),
            base_pose=p0,
            actuators=ThreadlikeActuator.tendons(active_tendon_routing),
            scale_rotational_basis_by_length=True,
            backend="auto",
        )

        self.q0 = jnp.zeros((self.model.num_internal_dofs,))
        self.q0dot = jnp.zeros((self.model.num_internal_dofs,))

        def _statics_eq(q, u):
            K = self.model.stiffness_matrix()
            B = self.model.actuation_matrix(q)
            G = self.model.gravitational_force(q)
            return K @ q + G - B @ u

        self._statics_eq_jit = eqx.filter_jit(_statics_eq)
        self._equilibrium_solver = optx.Newton(rtol=1e-6, atol=1e-6)

    def initial_state(self) -> SystemState:
        return SystemState(
            t=jnp.array(0.0),
            y=jnp.concatenate([self.q0, self.q0dot]),
            u=jnp.zeros((self.model.num_actuators,)),
        )

    def equilibrium(self, u: Array) -> Array:
        return optx.root_find(
            self._statics_eq_jit, self._equilibrium_solver, self.q0, u, max_steps=200
        ).value
