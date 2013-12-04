# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Linear_cell_complex_load_off )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

set(CMAKE_MODULE_PATH
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/
)
# add_definition(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)
include("CMakeLCCViewerQt.inc")

include_directories(include)









add_executable(load_off
	src/load_off.cc
)
target_link_libraries(load_off
	${MAP_VIEWER_LIBRARIES}
	CGAL
)

add_executable(load_tetmesh
	src/load_tetmesh.cc
)
target_link_libraries(load_tetmesh
	${MAP_VIEWER_LIBRARIES}
	CGAL
)

install(
  TARGETS
		load_off
    load_tetmesh
  DESTINATION
    bin/
)

