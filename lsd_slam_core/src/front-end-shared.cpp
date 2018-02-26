
#include "front-end-shared.h"


#include "DataStructures/FramePoseStruct.h"
#include "DataStructures/Frame.h"

#include "opencv2/opencv.hpp"
#include "timings.h"


using namespace lsd_slam;

void copyConfigToVariables(LSDSLAMConfig &lsdslamConfig)
{



    lsd_slam::doSlam = lsdslamConfig.poseoptim;
    lsd_slam::processEveryFrame = true;
    lsd_slam::packagePath = lsdslamConfig.package_path;

    lsd_slam::KFUsageWeight = lsdslamConfig.kfusage;
    lsd_slam::KFDistWeight = lsdslamConfig.kfdist;

    lsd_slam::useFabMap = lsdslamConfig.useFabMap;
    lsd_slam::minUseGrad = lsdslamConfig.minUseGrad;

    lsd_slam::doKFReActivation = lsdslamConfig.kfreactive;
    lsd_slam::useSubpixelStereo = lsdslamConfig.subpixelstereo;
    lsd_slam::maxLoopClosureCandidates = lsdslamConfig.maxLoopClosureCandiates;
}


// get camera calibration in form of an undistorter object.
// if no undistortion is required, the undistorter will just pass images through.
std::unique_ptr<Undistorter> generateUndistorter(LSDSLAMConfig &lsdslamConfig)
{
    return std::unique_ptr<Undistorter>(Undistorter::getUndistorterForFile(lsdslamConfig.camera_path.c_str()));
}

Sophus::Matrix3f generateCameraMatrix(const Undistorter * const undistorter)
{

    int w = undistorter->getOutputWidth();
    int h = undistorter->getOutputHeight();

    int w_inp = undistorter->getInputWidth();
    int h_inp = undistorter->getInputHeight();

    float fx = undistorter->getK().at<double>(0, 0);
    float fy = undistorter->getK().at<double>(1, 1);
    float cx = undistorter->getK().at<double>(2, 0);
    float cy = undistorter->getK().at<double>(2, 1);

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
void SlamBenchGrey2Distorded (std::vector<unsigned char>&  src , std::vector<unsigned char>& dest, const lsd_slam::Undistorter * const undistorter) {
    cv::Mat image;
    cv::Mat pre_image = cv::Mat(cv::Size2i(undistorter->getInputWidth(), undistorter->getInputHeight()), CV_8UC1, src.data());
    undistorter->undistort(pre_image, image);
    memcpy(dest.data(),image.data,undistorter->getInputWidth() * undistorter->getInputHeight());

}

void SlamBenchDepth2DistordedGrey (std::vector<sb_uchar3>&  src , std::vector<unsigned char>& dest, const lsd_slam::Undistorter * const undistorter) {
    cv::Mat image;
    cv::Mat pre_image = cv::Mat(cv::Size2i(undistorter->getInputWidth(), undistorter->getInputHeight()), CV_8UC3, src.data());
    cv::cvtColor(pre_image, pre_image, CV_RGB2GRAY);
    undistorter->undistort(pre_image, image);
    memcpy(dest.data(),image.data,undistorter->getInputWidth() * undistorter->getInputHeight());

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
