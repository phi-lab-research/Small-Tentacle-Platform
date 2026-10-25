# Motors (MAB Robotics)

Each tendon is pulled by one **[MA-D-GL60 v2.0](https://www.mabrobotics.pl/product-page/ma-d-gl60)**
actuator from MAB Robotics. The four motors are daisy-chained on CAN-FD and connected to the
PC through a **CANdle** USB adapter. Everything runs in **ROS 2** through the `candle_ros2` node.

## Main reference

**[MAB Ecosystem documentation](https://mabrobotics.github.io/MD80-x-CANdle-Documentation/intro.html)**:
setup, `mdtool` commands, control modes, ROS 2 node, errors and warnings.

This page adds lab-specific notes on top of it.

## Hardware

| Item | Value |
|---|---|
| Motor | MA-D-GL60 v2.0 |
| Motor IDs (current order) | `277`, `278`, `280`, `279` |
| Supply voltage | **24 V** (the motor accepts a wider range, see the website; 24 V is what we use) |
| Interface | CANdle (USB), 1 Mbps |
| Pulley radius | 32.5 mm (see [`data/soft_arm_parameters.yaml`](../../data/soft_arm_parameters.yaml)) |

> The motor order `[277, 278, 280, 279]` matches the current wiring of the motors to the
> tendons. Use this order in every command.

## Control mode

The motors are always used in **impedance mode**. The gains decide which component dominates:

| Configuration | Gains | Behaviour |
|---|---|---|
| Position-dominant | kp, kd > 0 | The motor tracks a target position like a spring-damper |
| Torque-dominant | kp = kd = 0 | Raw torque control |

`max_output` limits the motor output. **0.7** is the value used for this soft robot; lower
values such as 0.3 are too low.

## Scripts

The scripts are in [`scripts/`](scripts/). On the lab PC they are in `~/Desktop/md_motor/steps`.

| Script | What it does |
|---|---|
| [`INIT_4m.sh`](scripts/INIT_4m.sh) | Sets the zero of all four motors and starts the `candle_ros2` node |
| [`Impedance_motion_4motors/impedanceCTRL_motion4.sh`](scripts/Impedance_motion_4motors/impedanceCTRL_motion4.sh) | Adds and enables the motors in impedance mode, position-dominant (kp = 0.5, kd = 0.01) |
| [`Impedance_torque_4motor/impedanceCTRL_torque4.sh`](scripts/Impedance_torque_4motor/impedanceCTRL_torque4.sh) | Adds and enables the motors in impedance mode, torque-dominant (kp = kd = 0) |

---

## Before starting

1. **Power:** motors connected to the **24 V** power supply.
2. **CANdle:** connected to the PC and switched **on**.
3. **Tendons:** all four tendons are inside their pulleys and **in tension**. The tendons
   should be taut, but the **arm must not be deformed**.
4. **Warnings:** check that no motor has active errors or warnings (see below).
   If warnings are not cleared, the motors do not work properly.

> The zero of each motor is set at its **current position** during initialization, so the
> tendons must already be tensioned correctly at that point.

## Running the platform

Use three terminals.

**Terminal 1: mocap**

Launch the mocap package, see [`docs/mocap`](../mocap/README.md).

**Terminal 2: initialize and activate the motors**

```bash
cd ~/Desktop/md_motor/steps && ./INIT_4m.sh
```

This zeroes the motors (`mdtool config zero <ID>`) and starts the CANdle node
(`ros2 run candle_ros2 candle_ros2_node USB 1M`). Leave it running.

**Terminal 3: configure the motors, then run the controller**

First choose the control configuration:

```bash
# Position-dominant impedance
cd ~/Desktop/md_motor/steps/Impedance_motion_4motors && ./impedanceCTRL_motion4.sh

# or: torque-dominant impedance (kp = kd = 0)
cd ~/Desktop/md_motor/steps/Impedance_torque_4motor && ./impedanceCTRL_torque4.sh
```

Then run the controller. It is written in C++, so **always build and source before running it**:

```bash
cd ~/ros2_ws
colcon build
. install/setup.bash
ros2 run candle_ros2 PIDffd_CollocatedControl_CLIKrtMIXED
```

See [`code/`](../../code/README.md) for details.

## Stopping

Press **Ctrl+C in terminal 2** (the CANdle node): this automatically **disables** the motors.

To disable them explicitly:

```bash
ros2 service call /candle_ros2_node/disable_md80s candle_ros2/srv/GenericMd80Msg "{drive_ids:[277,278,280,279]}"
```

To stop the motors **without disabling** them, send zero torque:

```bash
ros2 topic pub /md80/motion_command candle_ros2/MotionCommand "{drive_ids: [277,278,280,279], target_position: [0.0,0.0,0.0,0.0], target_velocity: [0.0,0.0,0.0,0.0], target_torque: [0.0,0.0,0.0,0.0]}" --once
```

## Errors and warnings

Some warnings can appear (see the
[MAB documentation](https://mabrobotics.github.io/MD80-x-CANdle-Documentation/intro.html)
for their meaning). **They must be cleared, otherwise the motors do not work properly.**

```bash
mdtool ping all                  # list connected motors and their IDs
mdtool setup info <ID>           # show errors and warnings of one motor
mdtool clean warnings <ID>       # clear warnings
```
