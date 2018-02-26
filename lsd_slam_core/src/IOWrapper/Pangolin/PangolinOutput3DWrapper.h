/*
 * PangolinOutput3DWrapper.h
 *
 *  Created on: 17 Oct 2014
 *      Author: thomas
 */

#ifndef PANGOLINOUTPUT3DWRAPPER_H_
#define PANGOLINOUTPUT3DWRAPPER_H_

#include "IOWrapper/Output3DWrapper.h"
#include "Keyframe.h"
#include <map>
#include <SLAMBenchUI.h>
#include <values/Value.h>

namespace lsd_slam
{

class Frame;
class KeyFrameGraph;

struct GraphConstraint
{
    int from;
    int to;
    float err;
};

class PangolinOutput3DWrapper : public Output3DWrapper
{
    public:
        PangolinOutput3DWrapper(int width, int height);
        virtual ~PangolinOutput3DWrapper();

        virtual void publishKeyframeGraph(KeyFrameGraph* graph);

        // publishes a keyframe. if that frame already existis, it is overwritten, otherwise it is added.
        virtual void publishKeyframe(Frame* f);

        virtual slambench::values::PointCloudValue *  getMap();


        int publishLvl;

    private:
        int width, height;
        std::map<int, Keyframe *> keyframes;

};
}

#endif /* PANGOLINOUTPUT3DWRAPPER_H_ */
