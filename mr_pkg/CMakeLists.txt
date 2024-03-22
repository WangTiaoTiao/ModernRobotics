cmake_minimum_required(VERSION 3.8)
project(mr_pkg)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(rclcpp_action)
find_package(sensor_msgs REQUIRED)
find_package(mr_interface_pkg REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

# uncomment the following section in order to fill in
# further dependencies manually.
# find_package(<dependency> REQUIRED)

add_executable(exe8_node src/exe8.cpp)
add_executable(exe9_actionS_node src/exe9_actionS.cpp)
add_executable(exe9_actionC_node src/exe9_actionC.cpp)

install(TARGETS 
  exe8_node
  exe9_actionS_node
  exe9_actionC_node
  DESTINATION lib/${PROJECT_NAME}
)

ament_target_dependencies(exe8_node
  "rclcpp"
  "Eigen3"
  "sensor_msgs"
)

ament_target_dependencies(exe9_actionS_node
  "rclcpp"
  "Eigen3"
  "rclcpp_action"
  "mr_interface_pkg"
  "sensor_msgs"
)

ament_target_dependencies(exe9_actionC_node
  "rclcpp"
  "Eigen3"
  "rclcpp_action"
  "mr_interface_pkg"
  "sensor_msgs"
)

install(
  DIRECTORY urdf launch meshes rviz
  DESTINATION share/${PROJECT_NAME}  
)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

ament_package()
