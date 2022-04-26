cmake_minimum_required(VERSION 3.5)
project(NoiseCurve3D)
set (CMAKE_CXX_STANDARD 11)

SET(SOURCE_BASE_HEADER_FILES
  ${PROJECT_SOURCE_DIR}/helpers/GeomHelpers.h
  #${PROJECT_SOURCE_DIR}/helpers/MThelpers.h
  ${PROJECT_SOURCE_DIR}/helpers/maximalblurredsegment3D.h
  ${PROJECT_SOURCE_DIR}/helpers/maximalblurredsegment.h
  ${PROJECT_SOURCE_DIR}/helpers/CLI11.hpp
)

SET(SOURCE_BASE_FILES
  ${PROJECT_SOURCE_DIR}/helpers/GeomHelpers.cpp
  #${PROJECT_SOURCE_DIR}/helpers/MThelpers.cpp
  ${PROJECT_SOURCE_DIR}/helpers/maximalblurredsegment3D.cpp
  ${PROJECT_SOURCE_DIR}/helpers/maximalblurredsegment.cpp
)

#Required in DGtal
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})

add_definitions(-g)
#add_definitions(-DWITH_QT5)

include_directories("${PROJECT_SOURCE_DIR}/helpers/")

LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

add_executable(ATC3D main_ATC3D ${SOURCE_BASE_FILES} ${SOURCE_BASE_HEADER_FILES})
target_link_libraries(ATC3D ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})
