#!/bin/bash

mdtool config zero 277
mdtool config zero 278
mdtool config zero 279
mdtool config zero 280
ros2 run candle_ros2 candle_ros2_node USB 1M 
