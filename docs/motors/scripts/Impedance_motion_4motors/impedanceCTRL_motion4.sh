#!/bin/bash
# Impedance mode, position-dominant (kp, kd > 0)
ros2 service call /candle_ros2_node/add_md80s candle_ros2/srv/AddMd80s "{drive_ids:[277,278,280,279]}"
ros2 service call /candle_ros2_node/set_mode_md80s candle_ros2/srv/SetModeMd80s "{drive_ids: [277,278,280,279], mode:["IMPEDANCE","IMPEDANCE","IMPEDANCE","IMPEDANCE"]}"
ros2 topic pub /md80/impedance_command candle_ros2/msg/ImpedanceCommand "{drive_ids: [277,278,280,279], kp: [0.5,0.5,0.5,0.5], kd: [0.01,0.01,0.01,0.01], max_output: [0.7,0.7,0.7,0.7]}" --once
ros2 service call /candle_ros2_node/enable_md80s candle_ros2/srv/GenericMd80Msg "{drive_ids:[277,278,280,279]}"
