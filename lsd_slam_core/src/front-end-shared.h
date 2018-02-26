//
// Created by andrew on 18/04/15.
//

#ifndef FRONT_END_SHARED_H
#define FRONT_END_SHARED_H

#include <memory>

#include "lsdslamconfig.h"
#include "util/Undistorter.h"
#include "SlamSystem.h"



//Configuation
void copyConfigToVariables(LSDSLAMConfig &lsdslamConfig);
std::unique_ptr<lsd_slam::Undistorter> generateUndistorter(SLAMBenchLibraryHelper &lsdslamConfig);
Sophus::Matrix3f generateCameraMatrix(const lsd_slam::Undistorter * const undistorter);
void logFramePoses(lsd_slam::SlamSystem *system, std::ostream& log_path);
void SlamBenchGrey2Distorded (std::vector<unsigned char>&  src , std::vector<unsigned char>& dest ,const lsd_slam::Undistorter * const undistorter) ;
SE3 getCurrentPose(lsd_slam::SlamSystem *system);

#endif //FRONT_END_SHARED_H
