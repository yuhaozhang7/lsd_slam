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

#include "GlobalMapping/TrackableKeyFrameSearch.h"


#include "GlobalMapping/KeyFrameGraph.h"
#include "DataStructures/Frame.h"
#include "Tracking/SE3Tracker.h"

#include <timings.h>
#include <sys/time.h>
namespace lsd_slam
{


TrackableKeyFrameSearch::TrackableKeyFrameSearch(KeyFrameGraph* graph, int w, int h, Eigen::Matrix3f K)
: graph(graph)
{
	tracker = new SE3Tracker(w,h,K);

	fowX = 2 * atanf((float)((w / K(0,0)) / 2.0f));
	fowY = 2 * atanf((float)((h / K(1,1)) / 2.0f));

	msTrackPermaRef=0;
	nTrackPermaRef=0;
	nAvgTrackPermaRef=0;


}

TrackableKeyFrameSearch::~TrackableKeyFrameSearch()
{
	delete tracker;
}




std::vector<TrackableKFStruct> TrackableKeyFrameSearch::findEuclideanOverlapFrames(Frame* frame, float distanceTH, float angleTH, bool checkBothScales)
{
	// basically the maximal angle-difference in viewing direction is angleTH*(average FoV).
	// e.g. if the FoV is 130°, then it is angleTH*130°.
	float cosAngleTH = cosf(angleTH*0.5f*(fowX + fowY));


	Eigen::Vector3d pos = frame->getScaledCamToWorld().translation();
	Eigen::Vector3d viewingDir = frame->getScaledCamToWorld().rotationMatrix().rightCols<1>();

	std::vector<TrackableKFStruct> potentialReferenceFrames;

	float distFacReciprocal = 1;
	if(checkBothScales)
		distFacReciprocal = frame->meanIdepth / frame->getScaledCamToWorld().scale();

	// for each frame, calculate the rough score, consisting of pose, scale and angle overlap.

	for(unsigned int i=0;i<graph->keyframesAll.size();i++)
	{
		Eigen::Vector3d otherPos = graph->keyframesAll[i]->getScaledCamToWorld().translation();

		// get distance between the frames, scaled to fit the potential reference frame.
		float distFac = graph->keyframesAll[i]->meanIdepth / graph->keyframesAll[i]->getScaledCamToWorld().scale();
		if(checkBothScales && distFacReciprocal < distFac) distFac = distFacReciprocal;
		Eigen::Vector3d dist = (pos - otherPos) * distFac;
		float dNorm2 = dist.dot(dist);
		if(dNorm2 > distanceTH) continue;

		Eigen::Vector3d otherViewingDir = graph->keyframesAll[i]->getScaledCamToWorld().rotationMatrix().rightCols<1>();
		float dirDotProd = otherViewingDir.dot(viewingDir);
		if(dirDotProd < cosAngleTH) continue;

		potentialReferenceFrames.push_back(TrackableKFStruct());
		potentialReferenceFrames.back().ref = graph->keyframesAll[i];
		potentialReferenceFrames.back().refToFrame = se3FromSim3(graph->keyframesAll[i]->getScaledCamToWorld().inverse() * frame->getScaledCamToWorld()).inverse();
		potentialReferenceFrames.back().dist = dNorm2;
		potentialReferenceFrames.back().angle = dirDotProd;
	}


	return potentialReferenceFrames;
}




Frame* TrackableKeyFrameSearch::findRePositionCandidate(Frame* frame, float maxScore)
{

	std::vector<TrackableKFStruct> potentialReferenceFrames = findEuclideanOverlapFrames(frame, maxScore / (KFDistWeight*KFDistWeight), 0.75);



	float bestScore = maxScore;

	Frame* bestFrame = 0;
	SE3 bestRefToFrame = SE3();
	SE3 bestRefToFrame_tracked = SE3();

	int checkedSecondary = 0;
	for(unsigned int i=0;i<potentialReferenceFrames.size();i++)
	{
		if(frame->getTrackingParent() == potentialReferenceFrames[i].ref)
			continue;

		if(potentialReferenceFrames[i].ref->idxInKeyframes < INITIALIZATION_PHASE_COUNT)
			continue;

		struct timeval tv_start, tv_end;
		gettimeofday(&tv_start, NULL);
		tracker->checkPermaRefOverlap(potentialReferenceFrames[i].ref, potentialReferenceFrames[i].refToFrame);
		gettimeofday(&tv_end, NULL);
		msTrackPermaRef = 0.9*msTrackPermaRef + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
		nTrackPermaRef++;

		float score = getRefFrameScore(potentialReferenceFrames[i].dist, tracker->pointUsage);

		if(score < maxScore)
		{
			SE3 RefToFrame_tracked = tracker->trackFrameOnPermaref(potentialReferenceFrames[i].ref, frame, potentialReferenceFrames[i].refToFrame);
			Sophus::Vector3d dist = RefToFrame_tracked.translation() * potentialReferenceFrames[i].ref->meanIdepth;

			float newScore = getRefFrameScore(dist.dot(dist), tracker->pointUsage);
			float poseDiscrepancy = (potentialReferenceFrames[i].refToFrame * RefToFrame_tracked.inverse()).log().norm();
			float goodVal = tracker->pointUsage * tracker->lastGoodCount / (tracker->lastGoodCount+tracker->lastBadCount);
			checkedSecondary++;

			if(tracker->trackingWasGood && goodVal > relocalizationTH && newScore < bestScore && poseDiscrepancy < 0.2)
			{

				bestScore = score;
				bestFrame = potentialReferenceFrames[i].ref;
				bestRefToFrame = potentialReferenceFrames[i].refToFrame;
				bestRefToFrame_tracked = RefToFrame_tracked;

			}
		}
	}


	if(bestFrame != 0)
	{

		return bestFrame;
	}
	else
	{

		return 0;
	}
}

std::unordered_set<Frame*> TrackableKeyFrameSearch::findCandidates(Frame* keyframe, Frame* &, bool , bool closenessTH)
{
	std::unordered_set<Frame*> results;

	// Add all candidates that are similar in an euclidean sense.
	std::vector<TrackableKFStruct> potentialReferenceFrames = findEuclideanOverlapFrames(keyframe, closenessTH * 15 / (KFDistWeight*KFDistWeight), 1.0 - 0.25 * closenessTH, true);
	for(unsigned int i=0;i<potentialReferenceFrames.size();i++)
		results.insert(potentialReferenceFrames[i].ref);

	return results;
}

Frame* TrackableKeyFrameSearch::findAppearanceBasedCandidate(Frame* )
{

	return nullptr;

}


}
