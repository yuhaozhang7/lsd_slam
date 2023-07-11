cmake_minimum_required(VERSION 2.8)

set(lsd_VERSION cpp)

FIND_PACKAGE(GLUT)


find_package(OpenGL )

# --- common includes

include_directories(./lsd_slam_core/src/thirdparty/Sophus/)
include_directories(
  ${EIGEN3_INCLUDE_DIR}
  ${Pangolin_INCLUDE_DIR}
  ${G2O_INCLUDE_DIR}
  ./include
  ${CSPARSE_INCLUDE_DIR} #Has been set by SuiteParse
  ${CHOLMOD_INCLUDE_DIR} #Has been set by SuiteParse
  ${OpenCV_INCLUDE_DIRS}
  ${SLAMBENCH_INCLUDE_DIR}
)



# ------------ COMPILER FLAGS ------------ 
# Also add some useful compiler flag
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -O3 -std=c++11 -g -Wno-error=deprecated-declarations -ftemplate-backtrace-limit=0") 

IF(ANDROID)
ELSE(ANDROID)
# SSE flags
#add_definitions("-DENABLE_SSE")

# Also add some useful compiler flag
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -march=native") 

ENDIF(ANDROID)
# ------------ COMMON  ------------ 
set(G2O_TYPES_SIM3	/home/yuhao/SLAM/slambench/deps/g2o/lib/libg2o_types_sim3.so)
set(G2O_TYPES_SBA	/home/yuhao/SLAM/slambench/deps/g2o/lib/libg2o_types_sba.so)
set(G2O_CORE_LIBRARY	/home/yuhao/SLAM/slambench/deps/g2o/lib/libg2o_core.so)
set(G2O_SOLVER_CSPARSE	/home/yuhao/SLAM/slambench/deps/g2o/lib/libg2o_solver_csparse.so)
set(G2O_SOLVER_CSPARSE_EXTENSION	/home/yuhao/SLAM/slambench/deps/g2o/lib/libg2o_csparse_extension.so)
set(G2O_STUFF_LIBRARY	/home/yuhao/SLAM/slambench/deps/g2o/lib/libg2o_stuff.so)
set(CSPARSE_LIBRARY	/home/yuhao/SLAM/slambench/deps/suitesparse/lib/libcxsparse.a)

set(common_libraries     ${G2O_TYPES_SIM3} ${G2O_TYPES_SBA}   ${G2O_CORE_LIBRARY}  ${G2O_SOLVER_CSPARSE} ${G2O_SOLVER_CSPARSE_EXTENSION} ${G2O_STUFF_LIBRARY}  ${CSPARSE_LIBRARY} )
set(common_compile_flags "")


include_directories(./lsd_slam_core/src/)



# Set source files
set(lsd_SOURCE_FILES
  ./lsdslam_library.cpp
  ./lsd_slam_core/src/front-end-shared.cpp
  ./lsd_slam_core/src/DataStructures/Frame.cpp
  ./lsd_slam_core/src/DataStructures/FramePoseStruct.cpp
  ./lsd_slam_core/src/DataStructures/FrameMemory.cpp
  ./lsd_slam_core/src/SlamSystem.cpp
  ./lsd_slam_core/src/DepthEstimation/DepthMap.cpp
  ./lsd_slam_core/src/DepthEstimation/DepthMapPixelHypothesis.cpp
  ./lsd_slam_core/src/util/SophusUtil.cpp
  ./lsd_slam_core/src/util/settings.cpp
  ./lsd_slam_core/src/util/Undistorter.cpp
  ./lsd_slam_core/src/Tracking/Sim3Tracker.cpp
  ./lsd_slam_core/src/Tracking/Relocalizer.cpp
  ./lsd_slam_core/src/Tracking/SE3Tracker.cpp
  ./lsd_slam_core/src/Tracking/least_squares.cpp
  ./lsd_slam_core/src/Tracking/TrackingReference.cpp
  ./lsd_slam_core/src/GlobalMapping/KeyFrameGraph.cpp
  ./lsd_slam_core/src/GlobalMapping/g2oTypeSim3Sophus.cpp
  ./lsd_slam_core/src/GlobalMapping/TrackableKeyFrameSearch.cpp
)

set(SOURCE_FILES
  ${lsd_SOURCE_FILES}
)

set(CPP_UI_EXTRA_FILES lsd_slam_core/src/IOWrapper/Pangolin/PangolinOutput3DWrapper.cpp)
set(CPP_UI_EXTRA_INCLUDE .)

# build shared library.
add_library(${appname}-${lsd_VERSION} ${SOURCE_FILES} ${CPP_UI_EXTRA_FILES})

target_link_libraries(${appname}-${lsd_VERSION} ${common_libraries})
SET_TARGET_PROPERTIES(${appname}-${lsd_VERSION} PROPERTIES COMPILE_FLAGS "-fPIC")


 install(TARGETS ${appname}-cpp  DESTINATION  /home/yuhao/lib/)

SET(SLAMBENCH_C_WRAPPER                            -Wl,--whole-archive -L/home/yuhao/SLAM/slambench/build/lib slambench-c-wrapper                            -Wl,--no-whole-archive)
SET(NULL_FILE_PATH	"/home/yuhao/SLAM/slambench/build/null.cpp")

function(generate_slam_library library_name )

    add_library(${library_name}-library SHARED ${NULL_FILE_PATH})
    target_link_libraries(${library_name}-library PRIVATE ${SLAMBENCH_C_WRAPPER} ${library_name})
    SET_TARGET_PROPERTIES(${library_name}-library PROPERTIES LINK_FLAGS "${ARGN}")
    install(TARGETS ${library_name}-library  DESTINATION  /home/yuhao/lib/)

endfunction()

generate_slam_library(${appname}-cpp)



