/*
 * PangolinOutput3DWrapper.cpp
 *
 *  Created on: 17 Oct 2014
 *      Author: thomas
 */

#include "PangolinOutput3DWrapper.h"
#include "util/settings.h"
#include "DataStructures/Frame.h"
#include "GlobalMapping/KeyFrameGraph.h"
#include <values/Value.h>


#include <Eigen/Core>

namespace lsd_slam
{

PangolinOutput3DWrapper::PangolinOutput3DWrapper(int width, int height)
 :
           publishLvl(0),
           width(width),
           height(height)
{

}

PangolinOutput3DWrapper::~PangolinOutput3DWrapper()
{

}



void PangolinOutput3DWrapper::publishKeyframe(Frame* f)
{

    Keyframe * fMsg = new Keyframe;


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



    //Exists
    if(keyframes.find(fMsg->id) != keyframes.end())
    {
        keyframes[fMsg->id]->updatePoints(fMsg);
        delete fMsg;
    }
    else
    {
        fMsg->initId = keyframes.size();
        keyframes[fMsg->id] = fMsg;
    }

}

slambench::values::PointCloudValue *  PangolinOutput3DWrapper::getMap () {

	slambench::values::PointCloudValue *point_cloud = new slambench::values::PointCloudValue();

    for (auto tuple : keyframes) {
        //frame->pointData

        Keyframe* frame = tuple.second;
        InputPointDense * originalInput = (InputPointDense *)frame->pointData;

        float fx = frame->fx;
        float fy = frame->fy;
        float cx = frame->cx;
        float cy = frame->cy;

        float fxi = 1/fx;
        float fyi = 1/fy;
        float cxi = -cx / fx;
        float cyi = -cy / fy;

        float my_scaledTH = 1e-3;
        float my_absTH = 1e-1;
        float my_scale = frame->camToWorld.scale();
        int my_minNearSupport = 9;
        int my_sparsifyFactor = 1;

        for(int y = 1; y < height - 1; y++)
        {
            for(int x = 1; x < width - 1; x++)
            {


                if(originalInput[x + y * width].idepth <= 0)
                            continue;

                        if(my_sparsifyFactor > 1 && rand() % my_sparsifyFactor != 0)
                            continue;

                        float depth = 1 / originalInput[x + y * width].idepth;

                        float depth4 = depth * depth;

                        depth4 *= depth4;

                        if(originalInput[x + y * width].idepth_var * depth4 > my_scaledTH)
                            continue;

                        if(originalInput[x + y * width].idepth_var * depth4 * my_scale * my_scale > my_absTH)
                            continue;

                        if(my_minNearSupport > 1)
                        {
                            int nearSupport = 0;
                            for(int dx = -1; dx < 2; dx++)
                            {
                                for(int dy = -1; dy < 2; dy++)
                                {
                                    int idx = x + dx + (y + dy) * width;
                                    if(originalInput[idx].idepth > 0)
                                    {
                                        float diff = originalInput[idx].idepth - 1.0f / depth;
                                        if(diff * diff < 2 * originalInput[x + y * width].idepth_var)
                                            nearSupport++;
                                    }
                                }
                            }

                            if(nearSupport < my_minNearSupport)
                                continue;
                        }






                Eigen::Vector4f original = { (x * fxi + cxi) * depth , (y * fyi + cyi) * depth, depth, 1};
                Eigen::Matrix4f mat = frame->camToWorld.matrix().cast<float>();
                original =  mat * original;

                point_cloud->AddPoint(slambench::values::Point3DF(original[0] / original[3], original[1] / original[3] , original[2] / original[3]));

            }
        }

    }

    return point_cloud;
}

void PangolinOutput3DWrapper::publishKeyframeGraph(KeyFrameGraph* graph)
{


    int num = graph->keyframesAll.size();

    unsigned char * buffer = new unsigned char[num * sizeof(GraphFramePose)];

    GraphFramePose* framePoseData = (GraphFramePose*)buffer;

    for(unsigned int i = 0; i < graph->keyframesAll.size(); i++)
    {
        framePoseData[i].id = graph->keyframesAll[i]->id();
        memcpy(framePoseData[i].camToWorld, graph->keyframesAll[i]->getScaledCamToWorld().cast<float>().data(), sizeof(float) * 7);
    }



    for(int i = 0; i < num; i++)
     {
         if(keyframes.find(framePoseData[i].id) != keyframes.end())
         {
             memcpy(keyframes[framePoseData[i].id]->camToWorld.data(), &framePoseData[i].camToWorld[0], sizeof(float) * 7);
         }
     }

    delete [] buffer;
}



}
