# ============================================================================
# WARNING: RECENTLY UPDATED, NOT YET VALIDATED
#   Cleaned to build only the files currently in src/ (candle_ros2_node,
#   common_functions, PIDffd_CollocatedControl_CLIKrtMIXED). Entries for
#   .cpp files that no longer exist were removed.
# Double check that the package builds before using it.
# ============================================================================
cmake_minimum_required(VERSION 3.5)
project(candle_ros2)

# Default to C99
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()

# Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Dependencies
find_package(rclcpp REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(ament_cmake REQUIRED)
find_package(ament_index_cpp REQUIRED)
find_package(rosidl_default_generators REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(mocap_optitrack_interfaces REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)

rosidl_generate_interfaces(${PROJECT_NAME}
  "msg/ImpedanceCommand.msg"
  "msg/MotionCommand.msg"
  "msg/Pid.msg"
  "msg/VelocityPidCommand.msg"
  "msg/PositionPidCommand.msg"
  "srv/AddMd80s.srv"
  "srv/GenericMd80Msg.srv"
  "srv/SetModeMd80s.srv"
  "srv/SetLimitsMd80.srv"
  DEPENDENCIES builtin_interfaces std_msgs
)
set(CMAKE_INSTALL_RPATH $ORIGIN)  # sets proper runtime paths for exec to link .so properly

# CANdle library (MAB Robotics)
set(CANDLE_BUILD_STATIC ON)
set(CANDLE_OMIT_EXAMPLES ON)
add_subdirectory(candle)

# Soft robot model and control functions (shared by the controllers)
add_library(common_functions
  src/softrob_gvs.cpp
  src/general_fcns.cpp
  src/softrob_CCandCLIK.cpp
)
target_include_directories(common_functions PUBLIC ${CMAKE_SOURCE_DIR}/include ${EIGEN3_INCLUDE_DIRS})

# ---------------------------------------------------------------------------
# candle_ros2_node: MAB motor driver node
# ---------------------------------------------------------------------------
add_executable(${PROJECT_NAME}_node src/md80_node.cpp)
target_link_libraries(${PROJECT_NAME}_node candle)
target_include_directories(${PROJECT_NAME}_node PUBLIC ./candle/include ${CMAKE_SOURCE_DIR}/include)
rosidl_target_interfaces(${PROJECT_NAME}_node ${PROJECT_NAME} "rosidl_typesupport_cpp")
ament_target_dependencies(${PROJECT_NAME}_node rclcpp sensor_msgs)

# ---------------------------------------------------------------------------
# PIDffd_CollocatedControl_CLIKrtMIXED: soft arm controller
# ---------------------------------------------------------------------------
add_executable(PIDffd_CollocatedControl_CLIKrtMIXED src/PIDffd_CollocatedControl_CLIKrtMIXED.cpp)
target_link_libraries(PIDffd_CollocatedControl_CLIKrtMIXED candle common_functions)
target_include_directories(PIDffd_CollocatedControl_CLIKrtMIXED PUBLIC ./candle/include ${CMAKE_SOURCE_DIR}/include ${EIGEN3_INCLUDE_DIRS})
rosidl_target_interfaces(PIDffd_CollocatedControl_CLIKrtMIXED ${PROJECT_NAME} "rosidl_typesupport_cpp")
ament_target_dependencies(PIDffd_CollocatedControl_CLIKrtMIXED rclcpp sensor_msgs mocap_optitrack_interfaces Eigen3)

# ---------------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------------
install(TARGETS
  ${PROJECT_NAME}_node
  PIDffd_CollocatedControl_CLIKrtMIXED
  DESTINATION lib/${PROJECT_NAME}
)

install(TARGETS common_functions
  DESTINATION lib/${PROJECT_NAME}
)

# marker_ik_node: Python state estimator, installed with the modules it imports
install(PROGRAMS src/marker_ik_node.py DESTINATION lib/${PROJECT_NAME})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../state_estimation ${CMAKE_CURRENT_SOURCE_DIR}/../model
  DESTINATION lib/${PROJECT_NAME}
  PATTERN "__pycache__" EXCLUDE
)

ament_package()
