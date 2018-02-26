/*
 * Keyframe.h
 *
 *  Created on: 17 Oct 2014
 *      Author: thomas
 */

#ifndef KEYFRAME_H_
#define KEYFRAME_H_

#include "util/settings.h"
#include "sophus/sim3.hpp"
#include <iostream>

struct InputPointDense
{
    float idepth;
    float idepth_var;
    unsigned char color[4];
};

struct GraphFramePose
{
    int id;
    float camToWorld[7];
};

class Keyframe
{
    public:
        Keyframe()
         : pointData(0),
           points(0),
           needsUpdate(false)
        {}

        virtual ~Keyframe()
        {
            if(pointData)
                delete [] pointData;
        }

        struct MyVertex
        {
            float point[3];
            unsigned char color[4];
        };

        void updatePoints(Keyframe * newFrame)
        {
            if(pointData == 0)
            {
                pointData = new unsigned char[width * height * sizeof(InputPointDense)];
            }

            memcpy(pointData, newFrame->pointData, width * height * sizeof(InputPointDense));

            needsUpdate = true;
        }

        int id;
        int initId;
        uint64_t time;
        bool isKeyframe;

        Sophus::Sim3f camToWorld;

        float fx;
        float fy;
        float cx;
        float cy;
        int height;
        int width;

        unsigned char * pointData;

        int points;
        bool needsUpdate;
};


#endif /* KEYFRAME_H_ */
