/*
 * PangolinOutput3DWrapper.cpp
 *
 *  Created on: 17 Oct 2014
 *      Author: thomas
 */

#include "PangolinOutput3DWrapper.h"

#include "util/SophusUtil.h"
#include "util/settings.h"
#include "DataStructures/Frame.h"
#include "GlobalMapping/KeyFrameGraph.h"
#include "sophus/sim3.hpp"
#include "GlobalMapping/g2oTypeSim3Sophus.h"

namespace lsd_slam
{

PangolinOutput3DWrapper::PangolinOutput3DWrapper(int width, int height, SLAMBenchUI * gui)
 : width(width),
   height(height),
   publishLvl(0),
   gui(gui)
{

}

PangolinOutput3DWrapper::~PangolinOutput3DWrapper()
{

}

void PangolinOutput3DWrapper::updateImage(unsigned char * data)
{
    //gui->updateRGB(data);
}

void PangolinOutput3DWrapper::publishKeyframe(Frame* f)
{
    Keyframe * fMsg = new Keyframe;

    boost::shared_lock<boost::shared_mutex> lock = f->getActiveLock();

    fMsg->id = f->id();
    fMsg->time = f->timestamp();
    fMsg->isKeyframe = true;

    int w = f->width(publishLvl);
    int h = f->height(publishLvl);

    fMsg->camToWorld = f->getScaledCamToWorld().cast<float>();

    fMsg->fx = f->fx(publishLvl);
    fMsg->fy = f->fy(publishLvl);
    fMsg->cx = f->cx(publishLvl);
    fMsg->cy = f->cy(publishLvl);

    fMsg->width = w;
    fMsg->height = h;

    fMsg->pointData = new unsigned char[w * h * sizeof(InputPointDense)];

    InputPointDense * pc = (InputPointDense*)fMsg->pointData;

    const float* idepth = f->idepth(publishLvl);
    const float* idepthVar = f->idepthVar(publishLvl);
    const float* color = f->image(publishLvl);

    for(int idx = 0;idx < w * h; idx++)
    {
        pc[idx].idepth = idepth[idx];
        pc[idx].idepth_var = idepthVar[idx];
        pc[idx].color[0] = color[idx];
        pc[idx].color[1] = color[idx];
        pc[idx].color[2] = color[idx];
        pc[idx].color[3] = color[idx];
    }

    lock.unlock();

    //Exists
    if(keyframes.find(fMsg->id) != keyframes.end())
    {
        keyframes[fMsg->id]->updatePoints(fMsg);
        gui->addVisibleItem(keyframes[fMsg->id]);
        delete fMsg;
    }
    else
    {
        fMsg->initId = keyframes.size();
        keyframes[fMsg->id] = fMsg;
        if (fMsg->id > 5) gui->addVisibleItem(keyframes[fMsg->id]);
    }

}


void PangolinOutput3DWrapper::publishTrackedFrame(Frame* kf)
{
}

void PangolinOutput3DWrapper::publishKeyframeGraph(KeyFrameGraph* graph)
{
    graph->keyframesAllMutex.lock_shared();

    int num = graph->keyframesAll.size();

    unsigned char * buffer = new unsigned char[num * sizeof(GraphFramePose)];

    GraphFramePose* framePoseData = (GraphFramePose*)buffer;

    for(unsigned int i = 0; i < graph->keyframesAll.size(); i++)
    {
        framePoseData[i].id = graph->keyframesAll[i]->id();
        memcpy(framePoseData[i].camToWorld, graph->keyframesAll[i]->getScaledCamToWorld().cast<float>().data(), sizeof(float) * 7);
    }

    graph->keyframesAllMutex.unlock_shared();

    for(int i = 0; i < num; i++)
     {
         if(keyframes.find(framePoseData[i].id) != keyframes.end())
         {
             memcpy(keyframes[framePoseData[i].id]->camToWorld.data(), &framePoseData[i].camToWorld[0], sizeof(float) * 7);
         }
     }

    delete [] buffer;
}

void PangolinOutput3DWrapper::publishTrajectory(std::vector<Eigen::Matrix<float, 3, 1>> trajectory, std::string identifier)
{
    //TODO
}

void PangolinOutput3DWrapper::publishTrajectoryIncrement(Eigen::Matrix<float, 3, 1> pt, std::string identifier)
{
    //TODO
}

void PangolinOutput3DWrapper::publishDebugInfo(Eigen::Matrix<float, 20, 1> data)
{
    //TODO
}

}
