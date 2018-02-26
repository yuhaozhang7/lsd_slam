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

#include "Relocalizer.h"
#include "DataStructures/Frame.h"
#include "Tracking/SE3Tracker.h"

namespace lsd_slam
{


Relocalizer::Relocalizer(int w, int h, Eigen::Matrix3f K)
{


	this->w = w;
	this->h = h;
	this->K = K;

	KFForReloc.clear();
	nextRelocIDX = maxRelocIDX = 0;
	continueRunning = isRunning = false;

	hasResult = false;
	resultKF = 0;
	resultFrameID = 0;
	resultFrameToKeyframe = SE3();
}

Relocalizer::~Relocalizer()
{

}



void Relocalizer::updateCurrentFrame(std::shared_ptr<Frame> currentFrame)
{

	this->CurrentRelocFrame = currentFrame;
	maxRelocIDX = nextRelocIDX + KFForReloc.size();




}


void Relocalizer::run(std::vector<Frame*> &allKeyframesList)
{
    // make KFForReloc List
    KFForReloc.clear();
    for(unsigned int k=0;k < allKeyframesList.size(); k++)
    {
        // insert
        KFForReloc.push_back(allKeyframesList[k]);

        // swap with a random element
        int ridx = rand()%(KFForReloc.size());
        Frame* tmp = KFForReloc.back();
        KFForReloc.back() = KFForReloc[ridx];
        KFForReloc[ridx] = tmp;
    }
    nextRelocIDX=0;
    maxRelocIDX=KFForReloc.size();

    hasResult = false;
    continueRunning = true;
    isRunning = true;

    // start threads
    seqLoop();
}


void Relocalizer::getResult(Frame* &out_keyframe, std::shared_ptr<Frame> &frame, int &out_successfulFrameID, SE3 &out_frameToKeyframe)
{

	if(hasResult)
	{
		out_keyframe = resultKF;
		out_successfulFrameID = resultFrameID;
		out_frameToKeyframe = resultFrameToKeyframe;
		frame = resultRelocFrame;
	}
	else
	{
		out_keyframe = 0;
		out_successfulFrameID = -1;
		out_frameToKeyframe = SE3();
		frame.reset();
	}
}


void Relocalizer::seqLoop()
{


    SE3Tracker* tracker = new SE3Tracker(w,h,K);

    while(continueRunning)
    {
        std::cout << "Reloc iteration ... " << std::endl;
        // if got something: do it (unlock in the meantime)
        if(nextRelocIDX < maxRelocIDX && CurrentRelocFrame)
        {

            std::cout << " something to do ..  " << std::endl;

            Frame* todo = KFForReloc[nextRelocIDX%KFForReloc.size()];
            nextRelocIDX++;
            if(todo->neighbors.size() <= 2) continue;

            std::shared_ptr<Frame> myRelocFrame = CurrentRelocFrame;

            // initial Alignment
            SE3 todoToFrame = tracker->trackFrameOnPermaref(todo, myRelocFrame.get(), SE3());

            // try neighbours
            float todoGoodVal = tracker->pointUsage * tracker->lastGoodCount / (tracker->lastGoodCount+tracker->lastBadCount);
            if(todoGoodVal > relocalizationTH)
            {
                int numGoodNeighbours = 0;
                int numBadNeighbours = 0;

                float bestNeightbourGoodVal = todoGoodVal;

                Frame* bestKF = todo;
                SE3 bestKFToFrame = todoToFrame;
                for(Frame* nkf : todo->neighbors)
                {
                    SE3 nkfToFrame_init = se3FromSim3((nkf->getScaledCamToWorld().inverse() * todo->getScaledCamToWorld() * sim3FromSE3(todoToFrame.inverse(), 1))).inverse();
                    SE3 nkfToFrame = tracker->trackFrameOnPermaref(nkf, myRelocFrame.get(), nkfToFrame_init);

                    float goodVal = tracker->pointUsage * tracker->lastGoodCount / (tracker->lastGoodCount+tracker->lastBadCount);
                    if(goodVal > relocalizationTH*0.8 && (nkfToFrame * nkfToFrame_init.inverse()).log().norm() < 0.1)
                        numGoodNeighbours++;
                    else
                        numBadNeighbours++;

                    if(goodVal > bestNeightbourGoodVal)
                    {
                        bestNeightbourGoodVal = goodVal;
                        bestKF = nkf;
                        bestKFToFrame = nkfToFrame;

                    }
                }

                if(numGoodNeighbours > numBadNeighbours || numGoodNeighbours >= 5)
                {


                    // set everything to stop!
                    continueRunning = false;
                    resultRelocFrame = myRelocFrame;
                    resultFrameID = myRelocFrame->id();
                    resultKF = bestKF;
                    resultFrameToKeyframe = bestKFToFrame.inverse();

                    hasResult = true;

                }
                else
                {

                }
            }

        } else {
            continueRunning = false;
        }
    }

    delete tracker;
}

}
