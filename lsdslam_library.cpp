/**
* This file is part of LSD-SLAM.
*
* Copyright 2013 Jakob Engel <engelj at in dot tum dot de> (Technical University of Munich)
* For more information see <http://vision.in.tum.de/lsdslam>
*
* LSD-SLAM is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* LSD-SLAM is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with LSD-SLAM. If not, see <http://www.gnu.org/licenses/>.
*/


/*
 * This is a modified version of ... to provide a kfusion style interface.
 */

#include "util/settings.h"
#include "util/Parse.h"
#include "util/globalFuncs.h"
#include "front-end-shared.h"

#include "timings.h"

#include <io/SLAMFrame.h>
#include <io/sensor/CameraSensor.h>
#include <io/sensor/CameraSensorFinder.h>

#include <stdlib.h>
#include <GlobalMapping/KeyFrameGraph.h>
#include <DataStructures/FramePoseStruct.h>
#include <DataStructures/Frame.h>
#include <SLAMBenchUI.h>
#include <SLAMBenchAPI.h>
#include <IOWrapper/Pangolin/PangolinOutput3DWrapper.h>

using namespace lsd_slam;

#define LSDSLAM_ID "LSDSLAM"


/***************************
 * Default parameters
***************************** */

const double default_kfusage = 5.0;
const double default_kfdist = 5.0;
const float  default_minUseGrad = 5.0;
const bool   default_kfreactive = false;
const bool   default_subpixelstereo = true;
const double default_maxLoopClosureCandiates = 20;
const bool   default_poseoptim = true;
const bool   default_useFabMap = true;
const int    default_randSeed = 1;
const std::string    default_camera_path = "";
const std::string     default_package_path = "benchmarks/lsdslam/src/cpp/";




/*
 * Accepts commandline parameters and configures the system. Processing happens elsewhere.
 */


static std::unique_ptr<Undistorter>undistorter;
static Sophus::Matrix3f cameraCalibration ;
static SlamSystem * slamsystem;

static std::vector<unsigned char>     grey_raw;
static std::vector<unsigned char>     grey_distorded;

static sb_uint2 inputSize;
PangolinOutput3DWrapper* outputWrapper = NULL;

static slambench::io::CameraSensor *grey_sensor = nullptr;


// algo parameters


double local_kfusage;
double local_kfdist;

bool local_useFabMap;
float local_minUseGrad;
bool local_kfreactive;
bool local_subpixelstereo;
double local_maxLoopClosureCandiates;
bool local_poseoptim;
std::string local_camera_path;
std::string local_package_path;
int local_randSeed;


std::string ui_identifier;

static slambench::TimeStamp last_frame_timestamp;


static slambench::outputs::Output *pose_output;
static slambench::outputs::Output *pointcloud_output;

static slambench::outputs::Output *grey_frame_output;


bool sb_new_slam_configuration(SLAMBenchLibraryHelper * slam_settings) {

    slam_settings->addParameter(TypedParameter<std::string>("cam","camera-path", "This is a camera file (original version only)." , &local_camera_path,&default_camera_path));
    slam_settings->addParameter(TypedParameter<std::string>("pack","package-path", "Path to the package lsdslam/src/cpp/ (original version only)." , &local_package_path,&default_package_path));
    slam_settings->addParameter(TypedParameter<double>("","kfusage", "Determines how often Keyframes are taken, depending on the Overlap to the current Keyframe. Larger -> more Keyframes." , &local_kfusage,&default_kfusage));
    slam_settings->addParameter(TypedParameter<double>("","kfdist",  "Determines how often Keyframes are taken, depending on the Distance to the current Keyframe. Larger -> more Keyframes." , &local_kfdist,&default_kfdist));
    slam_settings->addParameter(TypedParameter<float>("","minUseGrad",  "Minimal Absolute Image Gradient for a Pixel to be used at all. Increase, if your camera has large image noise." ,&local_minUseGrad,&default_minUseGrad));
    slam_settings->addParameter(TypedParameter<bool>("","fabmap",  "Use fabmap (original only)." ,&local_useFabMap,&default_useFabMap));
    slam_settings->addParameter(TypedParameter<bool>("","kfreactive",  "Reactive key-frames from pose-graph." ,&local_kfreactive,&default_kfreactive));
    slam_settings->addParameter(TypedParameter<bool>("","subpixelstereo", "Use subpixel stereo" ,&local_subpixelstereo,&default_subpixelstereo));
    slam_settings->addParameter(TypedParameter<double>("","maxLoopClosureCandiates",   "Maximum number of loop closure candidates" ,&local_maxLoopClosureCandiates,&default_maxLoopClosureCandiates));
    slam_settings->addParameter(TypedParameter<int>("","randSeed",  "Set the random seed" , &local_randSeed,&default_randSeed));
    slam_settings->addParameter(TypedParameter<bool>("", "use-pose-optim", "Enable pose-graph optimisation and loop-closure detection" ,&local_poseoptim,&default_poseoptim));

    return true;
}

bool sb_init_slam_system(SLAMBenchLibraryHelper * slam_settings)  {

    lsd_slam::doSlam = local_poseoptim;

    lsd_slam::KFUsageWeight = local_kfusage;
    lsd_slam::KFDistWeight = local_kfdist;

    lsd_slam::minUseGrad = local_minUseGrad;
    lsd_slam::packagePath = local_package_path;
    lsd_slam::doKFReActivation = local_kfreactive;
    lsd_slam::useSubpixelStereo = local_subpixelstereo;
    lsd_slam::maxLoopClosureCandidates = local_maxLoopClosureCandiates;
    
    undistorter = generateUndistorter(*slam_settings);
    cameraCalibration = generateCameraMatrix(undistorter.get());
    
    
    // Handle to an LSD-SLAM instance
    slamsystem = new SlamSystem(undistorter->getOutputWidth(), undistorter->getOutputHeight(), cameraCalibration, doSlam);

    slamsystem->setVisualization(nullptr);
    
    /**
       * Retrieve RGB and Depth sensors,
       *  - check input_size are the same
       *  - check camera are the same
       *  - get input_file
       */

	int w = undistorter->getOutputWidth();
	int h = undistorter->getOutputHeight();

	grey_raw.resize(w * h);
	grey_distorded.resize(w * h);

	slambench::io::CameraSensorFinder sensor_finder;
	grey_sensor = sensor_finder.FindOne(slam_settings->get_sensors(), {{"camera_type", "grey"}});
	if (grey_sensor == nullptr) {
		std::cerr << "Invalid sensors found, Grey not found." << std::endl;
		
		delete slamsystem;
		
		return false;
	}
	
	// check sensor frame and pixel format
	if(grey_sensor->PixelFormat != slambench::io::pixelformat::G_I_8) {
		std::cerr << "Grey sensor is not in G_I_8 format" << std::endl;
		return false;
	}
	if(grey_sensor->FrameFormat != slambench::io::frameformat::Raster) {
		std::cerr << "Grey sensor is not in raster format" << std::endl;
	}

	inputSize   = make_sb_uint2(grey_sensor->Width,grey_sensor->Height);
 	pose_output = new slambench::outputs::Output("Pose", slambench::values::VT_POSE, true);

  	pointcloud_output = new slambench::outputs::Output("PointCloud", slambench::values::VT_POINTCLOUD, true);
  	pointcloud_output->SetKeepOnlyMostRecent(true);

  	grey_frame_output = new slambench::outputs::Output("Grey Frame", slambench::values::VT_FRAME);
  	grey_frame_output->SetKeepOnlyMostRecent(true);


  	slam_settings->GetOutputManager().RegisterOutput(pose_output);
  	slam_settings->GetOutputManager().RegisterOutput(pointcloud_output);
  	slam_settings->GetOutputManager().RegisterOutput(grey_frame_output);

  	// Initialize Vertices
  	outputWrapper = new PangolinOutput3DWrapper(undistorter->getOutputWidth(), undistorter->getOutputHeight());
  	slamsystem->setVisualization(outputWrapper);


	return true;
}

bool sb_update_frame (SLAMBenchLibraryHelper * , slambench::io::SLAMFrame* s) {
	if(s->FrameSensor == grey_sensor) {
		memcpy(grey_raw.data(), s->GetData(), s->GetSize());
		s->FreeData();
		last_frame_timestamp = s->Timestamp;
		return true;
	}
	
	return false;
}

bool sb_process_once (SLAMBenchLibraryHelper * )  {
    //move this above

    static int runningIDX = 0;
    static float fakeTimeStamp = 0.0;

        SlamBenchGrey2Distorded (grey_raw , grey_distorded ,   undistorter.get());



        if (runningIDX == 0) {
            slamsystem->randomInit(grey_distorded.data(), fakeTimeStamp, runningIDX);
        } else {
            slamsystem->trackFrame(grey_distorded.data(), runningIDX, false, fakeTimeStamp);
        }


          fakeTimeStamp += 0.03;
          runningIDX++;

          return true;
}

bool sb_get_pose     (Eigen::Matrix4f* mat) {
    SE3 pose =  getCurrentPose(slamsystem);
    *mat = pose.matrix().cast<float>();
    return true;
}

bool    sb_get_tracked  (bool* tracked)  {
    *tracked = slamsystem->trackingIsGood;
    return true;
}



bool sb_clean_slam_system() {

    if (!slamsystem->finalized and slamsystem->trackingIsGood) {
            slamsystem->finalize();
    }
    return true;
}


bool sb_update_outputs(SLAMBenchLibraryHelper *lib, const slambench::TimeStamp *ts) {
	(void)lib;
	(void)ts;

	if(pose_output->IsActive()) {
		// Get the current pose as an eigen matrix
		Eigen::Matrix4f matrix;
		sb_get_pose(&matrix);

		std::lock_guard<FastLock> lock (lib->GetOutputManager().GetLock());
		pose_output->AddPoint(last_frame_timestamp, new slambench::values::PoseValue(matrix));
	}





	if(pointcloud_output->IsActive()) {
		auto map = outputWrapper->getMap () ;

		// Take lock only after generating the map
		std::lock_guard<FastLock> lock (lib->GetOutputManager().GetLock());
		if (map)
			pointcloud_output->AddPoint(last_frame_timestamp, map);
	}

	if(grey_frame_output->IsActive()) {
		std::lock_guard<FastLock> lock (lib->GetOutputManager().GetLock());

		grey_frame_output->AddPoint(last_frame_timestamp, new slambench::values::FrameValue(inputSize.x, inputSize.y, slambench::io::pixelformat::G_I_8, (void*) &(grey_raw.at(0))));
	}


	return true;
}




