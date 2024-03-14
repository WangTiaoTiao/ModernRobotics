cmake_minimum_required(VERSION 3.8)
project(mr_pkg)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

# uncomment the following section in order to fill in
# further dependencies manually.
# find_package(<dependency> REQUIRED)

add_executable(node5 src/node5.cpp)

target_include_directories(node5 PUBLIC
  Eigen3::Eigen
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_features(node5 PUBLIC c_std_99 cxx_std_17)  # Require C99 and C++17

install(TARGETS 
  node5
  DESTINATION lib/${PROJECT_NAME}
)

install(
  DIRECTORY urdf launch meshes rviz
  DESTINATION share/${PROJECT_NAME}  
)

ament_package()
