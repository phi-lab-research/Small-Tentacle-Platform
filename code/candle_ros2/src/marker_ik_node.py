import sys
from collections import deque
from pathlib import Path

import jax.numpy as jnp
import numpy as np
import rclpy
from mocap_optitrack_interfaces.msg import MarkerArray
from rclpy.node import Node
from sensor_msgs.msg import JointState
from std_msgs.msg import Float64MultiArray

for d in (Path(__file__).resolve().parents[2], Path(__file__).resolve().parent):
    sys.path.insert(0, str(d / "state_estimation"))

from marker_ik import kalman_step, marker_positions, model, qd_finite_difference, solve


class MarkerIKNode(Node):
    def __init__(self):
        super().__init__("marker_ik_node")
        self.arm_ids = self.declare_parameter("arm_ids", [8, 65, 7, 331, 278]).value
        self.rest_s = self.declare_parameter("rest_s", 5.0).value
        self.rate_hz = self.declare_parameter("rate_hz", 40.0).value
        self.pulley_r = self.declare_parameter("pulley_r", 0.0325).value
        self.fd_window = self.declare_parameter("fd_window", 5).value
        self.sigma_a = self.declare_parameter("sigma_a", 100.0).value
        self.sigma_z = self.declare_parameter("sigma_z", 0.011).value

        n = model.num_internal_dofs
        self.q = jnp.zeros(n)
        self.Pf = np.asarray(marker_positions(self.q))
        self.theta = None
        self.t_first = None
        self.t_last = None
        self.rest = []
        self.R = None
        self.t_off = None
        self.buf = deque(maxlen=self.fd_window)
        self.x = None
        self.P = jnp.diag(jnp.array([self.sigma_z**2, 1e2]))

        self.create_subscription(JointState, "/md80/joint_states", self.on_joint_states, 10)
        self.create_subscription(MarkerArray, "/mocap_markers", self.on_markers, 10)
        self.pub_q = self.create_publisher(Float64MultiArray, "/q_hat", 10)
        self.pub_fd = self.create_publisher(Float64MultiArray, "/qd_fd", 10)
        self.pub_kf = self.create_publisher(Float64MultiArray, "/qd_kf", 10)

    def on_joint_states(self, msg):
        self.theta = np.asarray(msg.position)

    def on_markers(self, msg):
        t = msg.header.stamp.sec + 1e-9 * msg.header.stamp.nanosec
        by_id = {int(m.id): [m.position.x, m.position.y, m.position.z] for m in msg.markers}
        if self.theta is None or any(i not in by_id for i in self.arm_ids):
            return
        P = np.array([by_id[i] for i in self.arm_ids])

        if self.R is None:
            self.t_first = t if self.t_first is None else self.t_first
            self.rest.append(P)
            if t - self.t_first >= self.rest_s:
                Pr = np.mean(self.rest, axis=0)
                cf, cr = self.Pf.mean(0), Pr.mean(0)
                U, _, Vt = np.linalg.svd((self.Pf - cf).T @ (Pr - cr))
                self.R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
                self.t_off = cr - self.R @ cf
                self.get_logger().info("rest calibration done")
            return

        if self.t_last is not None and t - self.t_last < 1.0 / self.rate_hz:
            return
        dt = None if self.t_last is None else t - self.t_last
        self.t_last = t

        self.q = solve(self.q, jnp.asarray((P - self.t_off) @ self.R), jnp.asarray(self.pulley_r * self.theta))

        if self.x is None:
            self.x = jnp.stack([self.q, jnp.zeros_like(self.q)], axis=1)
        else:
            self.x, self.P = kalman_step(self.x, self.P, self.q, dt, self.sigma_a, self.sigma_z)

        self.buf.append((t, np.asarray(self.q)))
        qd_fd = jnp.zeros_like(self.q)
        if len(self.buf) == self.fd_window:
            tb, qb = map(np.array, zip(*self.buf))
            qd_fd = qd_finite_difference(jnp.asarray(qb), jnp.asarray(tb), self.fd_window)[-1]

        for pub, arr in ((self.pub_q, self.q), (self.pub_fd, qd_fd), (self.pub_kf, self.x[:, 1])):
            pub.publish(Float64MultiArray(data=np.asarray(arr, dtype=float).tolist()))


def main():
    rclpy.init()
    rclpy.spin(MarkerIKNode())


if __name__ == "__main__":
    main()
