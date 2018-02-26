
#include "front-end-shared.h"

#include "DataStructures/FramePoseStruct.h"
#include "DataStructures/Frame.h"

#include "timings.h"
#include "io/sensor/CameraSensor.h"
#include "io/sensor/CameraSensorFinder.h"
#include "SLAMBenchConfiguration.h"


using namespace lsd_slam;

void copyConfigToVariables(LSDSLAMConfig &lsdslamConfig)
{


    lsd_slam::doSlam = lsdslamConfig.poseoptim;

    lsd_slam::KFUsageWeight = lsdslamConfig.kfusage;
    lsd_slam::KFDistWeight = lsdslamConfig.kfdist;

    lsd_slam::minUseGrad = lsdslamConfig.minUseGrad;

    lsd_slam::doKFReActivation = lsdslamConfig.kfreactive;
    lsd_slam::useSubpixelStereo = lsdslamConfig.subpixelstereo;
    lsd_slam::maxLoopClosureCandidates = lsdslamConfig.maxLoopClosureCandiates;
}





// get camera calibration in form of an undistorter object.
// if no undistortion is required, the undistorter will just pass images through.
std::unique_ptr<Undistorter> generateUndistorter(SLAMBenchLibraryHelper &slam_settings)
{
	slambench::io::CameraSensorFinder sensor_finder;
    slambench::io::CameraSensor *grey_camera = sensor_finder.FindOne(slam_settings.get_sensors(), {{"camera_type", "grey"}});

	if (grey_camera == nullptr) {
		std::cerr << "Invalid sensors found, Grey not found." << std::endl;
		exit(1);
	}


    sb_float4 camera = {grey_camera->Intrinsics[0],grey_camera->Intrinsics[1],grey_camera->Intrinsics[2],grey_camera->Intrinsics[3]};
    return std::unique_ptr<Undistorter>(Undistorter::getUndistorterForVar(camera,grey_camera->Width,grey_camera->Height));
}

Sophus::Matrix3f generateCameraMatrix(const Undistorter * const undistorter)
{


    float fx = undistorter->getK()(0,0);
    float fy = undistorter->getK()(1,1);
    float cx = undistorter->getK()(2,0);
    float cy = undistorter->getK()(2,1);

    Sophus::Matrix3f cameraCalibration;
    cameraCalibration << fx, 0.0, cx, 0.0, fy, cy, 0.0, 0.0, 1.0;

    return cameraCalibration;
}

enum
{
    yuv_shift = 14,
    R2Y = 4899,
    G2Y = 9617,
    B2Y = 1868,
};

void RGB2GREY (std::vector<sb_uchar3>&  src, std::vector<unsigned char>& dst,  int w , int h) {

    int tab[256*3];
    int n = w * h ;
    const int coeffs[] = { R2Y, G2Y, B2Y };


    int b = 0, g = 0, r = (1 << (yuv_shift-1));
    int db = coeffs[2^2], dg = coeffs[1], dr = coeffs[2];

    for( int i = 0; i < 256; i++, b += db, g += dg, r += dr )
    {
        tab[i] = b;
        tab[i+256] = g;
        tab[i+512] = r;
    }


    const int* _tab = tab;
    for(int i = 0; i < n; i++) {
        dst[i] = (unsigned char)((_tab[src[i].x] + _tab[src[i].y+256] + _tab[src[i].z+512]) >> yuv_shift);



    }


}

void SlamBenchGrey2Distorded (std::vector<unsigned char>& grey_image , std::vector<unsigned char>& dest, const lsd_slam::Undistorter * const undistorter) {

    undistorter->undistort(grey_image.data(), dest.data() , undistorter->getInputWidth(), undistorter->getInputHeight());

}

void SlamBenchDepth2DistordedGrey (std::vector<sb_uchar3>&  src , std::vector<unsigned char>& dest, const lsd_slam::Undistorter * const undistorter) {

    std::vector<unsigned char> grey_image(undistorter->getInputWidth() *  undistorter->getInputHeight());

    RGB2GREY(src,grey_image, undistorter->getInputWidth(), undistorter->getInputHeight() );

    undistorter->undistort(grey_image.data(), dest.data() , undistorter->getInputWidth(), undistorter->getInputHeight());

}


SE3 getCurrentPose(SlamSystem *system) {
    FramePoseStruct * it = system->getAllPoses().back();
    return se3FromSim3(it->getCamToWorld());
}

void writeFramePoseLogLine(std::ostream &out, SE3 &pose, int id)
{
    out << id << " " << pose.translation().transpose()
        << " " << pose.so3().unit_quaternion().x()
        << " " << pose.so3().unit_quaternion().y()
        << " " << pose.so3().unit_quaternion().z()
        << " " << pose.so3().unit_quaternion().w() << std::endl;
}

/*
 * Log the camera poses using the TUM dataset format
 */
void logFramePoses(SlamSystem *system, std::ostream &  log_stream)
{
    //std::vector<FramePoseStruct*> poses = system->getAllPoses();


    //for (int i = 0; i < poses.size(); ++i)

    SE3 prevPose;

    int i = 0;
    for (auto it : system->getAllPoses())
    {
        SE3 pose = se3FromSim3(it->getCamToWorld());

        while (i < it->frameID) {
            std::cout << "Skipped Frame: " << i << std::endl;
            writeFramePoseLogLine(log_stream, prevPose, i);
            ++i;
        }

        writeFramePoseLogLine(log_stream, pose, i);
        ++i;
    }

}
void synchroniseDevices() {};
