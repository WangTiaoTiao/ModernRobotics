cmake_minimum_required(VERSION 3.8)
project(mr_pkg)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(sensor_msgs REQUIRED)


include_directories(${EIGEN3_INCLUDE_DIR})

# uncomment the following section in order to fill in
# further dependencies manually.
# find_package(<dependency> REQUIRED)

add_executable(exe8_node src/exe8.cpp)

install(TARGETS 
  exe8_node
  DESTINATION lib/${PROJECT_NAME}
)

ament_target_dependencies(
  exe8_node
  "rclcpp"
  "Eigen3"
  "sensor_msgs"
)

target_compile_features(exe8_node PUBLIC c_std_99 cxx_std_17)  # Require C99 and C++17

install(
  DIRECTORY urdf launch meshes rviz
  DESTINATION share/${PROJECT_NAME}  
)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")


ament_package()
