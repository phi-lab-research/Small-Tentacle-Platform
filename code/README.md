# Controller Code

Base version of the controller used to run the soft arm. It is a C++ ROS 2 node
built inside the MAB Robotics motor package
[`candle_ros2`](https://github.com/mabrobotics/candle_ros2).

> **⚠️ Recently updated, not yet validated.**
> The following were recently changed and have **not been tested on the robot yet**.
> Double check them before use:
> - `PIDffd_CollocatedControl_CLIKrtMIXED.cpp`: motor command saturation added; tendon
>   insets, friction coefficients and rest cable lengths (`Lc0`) updated to the robot 3 values.
> - `CMakeLists.txt`: cleaned to build only the files in `src/` (motor node, support
>   library and main controller); entries for removed `.cpp` files were deleted.

## Workspace layout

Both packages go in the `src/` folder of the same ROS 2 workspace:

```
ros2_ws/
└── src/
    ├── candle_ros2/             ← this folder (MAB motor package + controller)
    └── ros2-mocap_optitrack/    ← clone from https://github.com/phi-lab-research/ros2-mocap_optitrack
```

The mocap package is not copied here: clone it directly from its repository
(see [`docs/mocap`](../docs/mocap/README.md)). The controller needs its
`mocap_optitrack_interfaces` messages to build.

## What is in `candle_ros2/`

| File | Role |
|---|---|
| `src/PIDffd_CollocatedControl_CLIKrtMIXED.cpp` | **Main controller node** |
| `src/softrob_gvs.cpp` | Support: GVS soft robot model |
| `src/softrob_CCandCLIK.cpp` | Support: collocated control and CLIK |
| `src/general_fcns.cpp` | Support: general math functions |
| `src/md80_node.cpp` | MAB motor driver node (`candle_ros2_node`) |
| `src/marker_ik_node.py` | State estimation node (Python): `q` and `qdot` from the markers |
| `include/` | Headers for the files above |
| `candle/` | MAB CANdle library |
| `msg/`, `srv/` | MAB motor messages and services |

Next to `candle_ros2/` in `code/`:

| File | Role |
|---|---|
| `state_estimation/marker_ik.py` | Marker IK solver and `qdot` estimators used by the node |
| `model/small_tentacle_platform_soromox.py` | soromox model of the arm |
| `pyproject.toml`, `uv.lock` | Python environment of the state estimation node |

## Build

Always build and source before running:

```bash
cd ~/ros2_ws
colcon build
. install/setup.bash
```

## Run

Start the mocap and the motors first (see [`docs/motors`](../docs/motors/README.md)), then:

```bash
ros2 run candle_ros2 PIDffd_CollocatedControl_CLIKrtMIXED
```

## Before running: check in the code

- **Marker IDs:** the arm's markers are free markers and their IDs change when they are
  lost. Check `base_id`, `marker3_id` and `tip_id` in the node constructor against Motive.
- **Motor order:** `[277, 278, 280, 279]`.
- **Reference trajectory:** the target points `xd` are set in `timerCallback()`.

## Safety: command saturation

Before being sent to the motors, the commands are saturated:

| Command | Limit | Constant |
|---|---|---|
| Position | [−π/2, 0] rad (negative = pulling) | `POS_SAT` |
| Torque | [−0.6, +0.6] N·m | `TAU_SAT` |

If a motor's position command is at its upper bound (0, fully released), its upper torque
limit is set to 0, so it cannot push further in the releasing direction.

> The position range can be extended up to [−π, 0] rad (`POS_SAT = M_PI`) **only if
> strictly necessary**: this is borderline for the arm.

## Physical parameters

The arm parameters used in the code (lengths, radii, material, tendon insets, friction,
rest cable lengths, pulley radius) are documented in
[`data/soft_arm_parameters.yaml`](../data/soft_arm_parameters.yaml).
They are hard-coded in the code: if you change a value in the YAML file, update the code
as well. Most are in `timerCallback()` of the main controller; the number of Gauss
integration points per section (`g_pt`) is in `src/softrob_CCandCLIK.cpp`.
