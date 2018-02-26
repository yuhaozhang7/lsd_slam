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

#include "SlamSystem.h"

#include "DataStructures/Frame.h"
#include "Tracking/SE3Tracker.h"
#include "Tracking/Sim3Tracker.h"
#include "DepthEstimation/DepthMap.h"
#include "Tracking/TrackingReference.h"
#include "util/globalFuncs.h"
#include "GlobalMapping/KeyFrameGraph.h"
#include "GlobalMapping/TrackableKeyFrameSearch.h"
#include "GlobalMapping/g2oTypeSim3Sophus.h"
#include "IOWrapper/Output3DWrapper.h"
#include <g2o/core/robust_kernel_impl.h>
#include "DataStructures/FrameMemory.h"
#include "deque"


#include "unistd.h"

// for mkdir
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>


using namespace lsd_slam;


SlamSystem::SlamSystem(int w, int h, Eigen::Matrix3f K, bool enableSLAM)
: SLAMEnabled(enableSLAM), finalized(false), relocalizer(w,h,K)
{
	if(w%16 != 0 || h%16!=0)
	{
		printf("image dimensions must be multiples of 16! Please crop your images / video accordingly.\n");
		assert(false);
	}

	this->width = w;
	this->height = h;
	this->K = K;
	trackingIsGood = true;


	currentKeyFrame =  nullptr;
	trackingReferenceFrameSharedPT = nullptr;
	keyFrameGraph = new KeyFrameGraph();
	createNewKeyFrame = false;

	map =  new DepthMap(w,h,K);
	
	newConstraintAdded = false;
	haveUnmergedOptimizationOffset = false;


	tracker = new SE3Tracker(w,h,K);
	// Do not use more than 4 levels for odometry tracking
	for (int level = 4; level < PYRAMID_LEVELS; ++level)
		tracker->settings.maxItsPerLvl[level] = 0;
	trackingReference = new TrackingReference();
	mappingTrackingReference = new TrackingReference();


	if(SLAMEnabled)
	{
		trackableKeyFrameSearch = new TrackableKeyFrameSearch(keyFrameGraph,w,h,K);
		constraintTracker = new Sim3Tracker(w,h,K);
		constraintSE3Tracker = new SE3Tracker(w,h,K);
		newKFTrackingReference = new TrackingReference();
		candidateTrackingReference = new TrackingReference();
	}
	else
	{
		constraintSE3Tracker = 0;
		trackableKeyFrameSearch = 0;
		constraintTracker = 0;
		newKFTrackingReference = 0;
		candidateTrackingReference = 0;
	}


	outputWrapper = 0;

	doFinalOptimization = false;
	depthMapScreenshotFlag = false;
	lastTrackingClosenessScore = 0;






	msTrackFrame = msOptimizationIteration = msFindConstraintsItaration = msFindReferences = 0;
	nTrackFrame = nOptimizationIteration = nFindConstraintsItaration = nFindReferences = 0;
	nAvgTrackFrame = nAvgOptimizationIteration = nAvgFindConstraintsItaration = nAvgFindReferences = 0;
	gettimeofday(&lastHzUpdate, NULL);

}

SlamSystem::~SlamSystem()
{

	// make sure none is waiting for something.
	printf("... waiting for SlamSystem's threads to exit\n");


	printf("DONE waiting for SlamSystem's threads to exit\n");

	if(trackableKeyFrameSearch != 0) delete trackableKeyFrameSearch;
	if(constraintTracker != 0) delete constraintTracker;
	if(constraintSE3Tracker != 0) delete constraintSE3Tracker;
	if(newKFTrackingReference != 0) delete newKFTrackingReference;
	if(candidateTrackingReference != 0) delete candidateTrackingReference;

	delete mappingTrackingReference;
	delete map;
	delete trackingReference;
	delete tracker;

	// make shure to reset all shared pointers to all frames before deleting the keyframegraph!
	unmappedTrackedFrames.clear();
	latestFrameTriedForReloc.reset();
	latestTrackedFrame.reset();
	currentKeyFrame.reset();
	trackingReferenceFrameSharedPT.reset();

	// delte keyframe graph
	delete keyFrameGraph;

	FrameMemory::getInstance().releaseBuffes();


}

void SlamSystem::setVisualization(Output3DWrapper* outputWrapper)
{
	this->outputWrapper = outputWrapper;
}

void SlamSystem::mergeOptimizationOffset()
{
	// update all vertices that are in the graph!


	bool needPublish = false;
	if(haveUnmergedOptimizationOffset)
	{

		for(unsigned int i=0;i<keyFrameGraph->keyframesAll.size(); i++)
			keyFrameGraph->keyframesAll[i]->pose->applyPoseGraphOptResult();


		haveUnmergedOptimizationOffset = false;
		needPublish = true;
	}


	if(needPublish)
		publishKeyframeGraph();
}





void SlamSystem::finalize()
{
    finalized = true;

	if (!SLAMEnabled)
	{
		std::cout << "Skipping finalisation, since we're not SLAM'ing" << std::endl;
		return;
	}

	printf("Finalizing Graph... finding final constraints!!\n");


	  doFullReConstraintTrack = true;
	  if(doFullReConstraintTrack)
	        {

	            std::cout << "Optimising Full Map " << std::endl;

	            int added = 0;
	            for(unsigned int i=0;i<keyFrameGraph->keyframesAll.size();i++)
	            {
	                if(keyFrameGraph->keyframesAll[i]->pose->isInGraph)
	                    added += findConstraintsForNewKeyFrames(keyFrameGraph->keyframesAll[i], false,  1.0);
	            }

	            printf("Done optizing Full Map! Added %d constraints.\n", added);

	            doFullReConstraintTrack = false;

	            lastNumConstraintsAddedOnFullRetrack = added;



	        }




	            newContraintsOptimised = false;
	            newConstraintAdded = true;
	            while(optimizationIteration(5, 0.02));
	            newContraintsOptimised = true;









	printf("Finalizing Graph... optimizing!!\n");



	 printf("doing final optimization iteration!\n");
	 optimizationIteration(50, 0.001);

	 while(optimizationIteration(5, 0.02));

	 newContraintsOptimised = true;





	while(doFinalOptimization)
	{
		usleep(200000);
	}

	printf("Finalizing Graph... publishing!!\n");

    doMappingIteration();


	while(doFinalOptimization)
	{
		usleep(200000);
	}


	usleep(200000);
	printf("Done Finalizing Graph.!!\n");
}



void SlamSystem::publishKeyframeGraph()
{
	if (outputWrapper != nullptr)
		outputWrapper->publishKeyframeGraph(keyFrameGraph);
}

void SlamSystem::requestDepthMapScreenshot(const std::string& filename)
{
	depthMapScreenshotFilename = filename;
	depthMapScreenshotFlag = true;
}

void SlamSystem::finishCurrentKeyframe()
{

	map->finalizeKeyFrame();

	if(SLAMEnabled)
	{
		mappingTrackingReference->importFrame(currentKeyFrame.get());
		currentKeyFrame->setPermaRef(mappingTrackingReference);
		mappingTrackingReference->invalidate();

		if(currentKeyFrame->idxInKeyframes < 0)
		{

			currentKeyFrame->idxInKeyframes = keyFrameGraph->keyframesAll.size();
			keyFrameGraph->keyframesAll.push_back(currentKeyFrame.get());
			keyFrameGraph->totalPoints += currentKeyFrame->numPoints;
			keyFrameGraph->totalVertices ++;



			newKeyFrames.push_back(currentKeyFrame.get());



			/////////////// INCLUDE CONSTRAINT STUFF

		    bool doneSomething = false;

		    while (newKeyFrames.size() != 0) {

		            Frame* newKF = newKeyFrames.front();
		            newKeyFrames.pop_front();


		            struct timeval tv_start, tv_end;
		            gettimeofday(&tv_start, NULL);


		            findConstraintsForNewKeyFrames(newKF, true, 1.0);

		            gettimeofday(&tv_end, NULL);
		            msFindConstraintsItaration = 0.9*msFindConstraintsItaration + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
		            nFindConstraintsItaration++;

		            FrameMemory::getInstance().pruneActiveFrames();


		             doneSomething = true;
		    }

		        if(doFullReConstraintTrack)
		        {

		            std::cout << "Optimising Full Map " << std::endl;

		            int added = 0;
		            for(unsigned int i=0;i<keyFrameGraph->keyframesAll.size();i++)
		            {
		                if(keyFrameGraph->keyframesAll[i]->pose->isInGraph)
		                    added += findConstraintsForNewKeyFrames(keyFrameGraph->keyframesAll[i], false,  1.0);
		            }

		            printf("Done optizing Full Map! Added %d constraints.\n", added);

		            doFullReConstraintTrack = false;

		            lastNumConstraintsAddedOnFullRetrack = added;

		            doneSomething = true;
		        }


		        if (doneSomething)
		        {

		            newContraintsOptimised = false;
		            newConstraintAdded = true;
		            while(optimizationIteration(5, 0.02));
		            newContraintsOptimised = true;


		        }










		}
	}

	if(outputWrapper!= 0)
		outputWrapper->publishKeyframe(currentKeyFrame.get());
}

void SlamSystem::discardCurrentKeyframe()
{

	if(currentKeyFrame->idxInKeyframes >= 0)
	{
		printf("WARNING: trying to discard a KF that has already been added to the graph... finalizing instead.\n");
		finishCurrentKeyframe();
		return;
	}


	map->invalidate();


	for(FramePoseStruct* p : keyFrameGraph->allFramePoses)
	{
		if(p->trackingParent != 0 && p->trackingParent->frameID == currentKeyFrame->id())
			p->trackingParent = 0;
	}

	keyFrameGraph->idToKeyFrame.erase(currentKeyFrame->id());


}

void SlamSystem::createNewCurrentKeyframe(std::shared_ptr<Frame> newKeyframeCandidate)
{


	if(SLAMEnabled)
	{
		// add NEW keyframe to id-lookup

		keyFrameGraph->idToKeyFrame.insert(std::make_pair(newKeyframeCandidate->id(), newKeyframeCandidate));

	}

	// propagate & make new.
	map->createKeyFrame(newKeyframeCandidate.get());


	//std::cout << __LINE__ << " New Keyframe: "  << newKeyframeCandidate->id() << std::endl;
	currentKeyFrame = newKeyframeCandidate;

}
void SlamSystem::loadNewCurrentKeyframe(Frame* keyframeToLoad)
{

	map->setFromExistingKF(keyframeToLoad);


	//std::cout << __LINE__ << " New Keyframe" << std::endl;
	currentKeyFrame = keyFrameGraph->idToKeyFrame.find(keyframeToLoad->id())->second;
	currentKeyFrame->depthHasBeenUpdatedFlag = false;

}

void SlamSystem::changeKeyframe(bool noCreate, bool force, float maxScore)
{
	Frame* newReferenceKF=0;
	std::shared_ptr<Frame> newKeyframeCandidate = latestTrackedFrame;
	if(doKFReActivation && SLAMEnabled)
	{
		struct timeval tv_start, tv_end;
		gettimeofday(&tv_start, NULL);
		newReferenceKF = trackableKeyFrameSearch->findRePositionCandidate(newKeyframeCandidate.get(), maxScore);
		gettimeofday(&tv_end, NULL);
		msFindReferences = 0.9*msFindReferences + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
		nFindReferences++;
	}

	if(newReferenceKF != 0)
		loadNewCurrentKeyframe(newReferenceKF);
	else
	{
		if(force)
		{
			if(noCreate)
			{
				trackingIsGood = false;
				nextRelocIdx = -1;
				printf("mapping is disabled & moved outside of known map. Starting Relocalizer!\n");
			}
			else
				createNewCurrentKeyframe(newKeyframeCandidate);
		}
	}


	createNewKeyFrame = false;
}

bool SlamSystem::updateKeyframe()
{
	std::shared_ptr<Frame> reference = nullptr;
	std::deque< std::shared_ptr<Frame> > references;


	int count = 0;

	// remove frames that have a different tracking parent.
	while(unmappedTrackedFrames.size() > 0 &&
			(!unmappedTrackedFrames.front()->hasTrackingParent() ||
					unmappedTrackedFrames.front()->getTrackingParent() != currentKeyFrame.get()))
	{
		unmappedTrackedFrames.front()->clear_refPixelWasGood();

		//std::cout << "Dropping: " << unmappedTrackedFrames.front()->id() << std::endl;

		unmappedTrackedFrames.pop_front();
		++count;
	}

	//std::cout << "Removing Frame Count: " << count << std::endl;

	// clone list
	if(unmappedTrackedFrames.size() > 0)
	{
		//std::cout << "updating depth" << std::endl;
		for(unsigned int i=0;i<unmappedTrackedFrames.size(); i++)
			references.push_back(unmappedTrackedFrames[i]);

		std::shared_ptr<Frame> popped = unmappedTrackedFrames.front();
		unmappedTrackedFrames.pop_front();





		map->updateKeyframe(references);

		popped->clear_refPixelWasGood();
		references.clear();
	}
	else
	{
		//std::cout << "NOT updating depth" << std::endl;

		return false;
	}






	if(outputWrapper != 0 && continuousPCOutput && currentKeyFrame != 0)
		outputWrapper->publishKeyframe(currentKeyFrame.get());

	return true;
}


void SlamSystem::addTimingSamples()
{
	map->addTimingSample();
	struct timeval now;
	gettimeofday(&now, NULL);
	float sPassed = ((now.tv_sec-lastHzUpdate.tv_sec) + (now.tv_usec-lastHzUpdate.tv_usec)/1000000.0f);
	if(sPassed > 1.0f)
	{
		nAvgTrackFrame = 0.8*nAvgTrackFrame + 0.2*(nTrackFrame / sPassed); nTrackFrame = 0;
		nAvgOptimizationIteration = 0.8*nAvgOptimizationIteration + 0.2*(nOptimizationIteration / sPassed); nOptimizationIteration = 0;
		nAvgFindReferences = 0.8*nAvgFindReferences + 0.2*(nFindReferences / sPassed); nFindReferences = 0;

		if(trackableKeyFrameSearch != 0)
		{
			trackableKeyFrameSearch->nAvgTrackPermaRef = 0.8*trackableKeyFrameSearch->nAvgTrackPermaRef + 0.2*(trackableKeyFrameSearch->nTrackPermaRef / sPassed); trackableKeyFrameSearch->nTrackPermaRef = 0;
		}
		nAvgFindConstraintsItaration = 0.8*nAvgFindConstraintsItaration + 0.2*(nFindConstraintsItaration / sPassed); nFindConstraintsItaration = 0;
		nAvgOptimizationIteration = 0.8*nAvgOptimizationIteration + 0.2*(nOptimizationIteration / sPassed); nOptimizationIteration = 0;

		lastHzUpdate = now;



	}

}



void SlamSystem::takeRelocalizeResult()
{
	Frame* keyframe;
	int succFrameID;
	SE3 succFrameToKF_init;
	std::shared_ptr<Frame> succFrame;
	relocalizer.getResult(keyframe, succFrame, succFrameID, succFrameToKF_init);
	assert(keyframe != 0);

	loadNewCurrentKeyframe(keyframe);

	trackingReference->importFrame(currentKeyFrame.get());
	trackingReferenceFrameSharedPT = currentKeyFrame;


	tracker->trackFrame(
			trackingReference,
			succFrame.get(),
			succFrameToKF_init);

	if(!tracker->trackingWasGood || tracker->lastGoodCount / (tracker->lastGoodCount + tracker->lastBadCount) < 1-0.75f*(1-MIN_GOODPERGOODBAD_PIXEL))
	{
	    std::cout << "Invalidate !!!!! " << std::endl;
		trackingReference->invalidate();
	}
	else
	{
		keyFrameGraph->addFrame(succFrame.get());

		if(unmappedTrackedFrames.size() < 50)
			unmappedTrackedFrames.push_back(succFrame);


		createNewKeyFrame = false;
		trackingIsGood = true;
	}
}

bool SlamSystem::doMappingIteration()
{

	if(currentKeyFrame == 0)
		return false;

	if(!doMapping && currentKeyFrame->idxInKeyframes < 0)
	{
		if(currentKeyFrame->numMappedOnThisTotal >= MIN_NUM_MAPPED)
			finishCurrentKeyframe();
		else
			discardCurrentKeyframe();

		map->invalidate();
		printf("Finished KF %d as Mapping got disabled!\n",currentKeyFrame->id());

		changeKeyframe(true, true, 1.0f);
	}

	mergeOptimizationOffset();
	addTimingSamples();

	if(dumpMap)
	{
		printf("Error inside doMappingIteration, with dumpMap == true \n");
	    exit(1);
		dumpMap = false;
	}


	// set mappingFrame
	if(trackingIsGood)
	{
		if(!doMapping)
		{
			printf("tryToChange refframe, lastScore %f!\n", lastTrackingClosenessScore);
			if(lastTrackingClosenessScore > 1)
				changeKeyframe(true, false, lastTrackingClosenessScore * 0.75);



			return false;
		}

		bool ret;

		if (createNewKeyFrame)
		{

			finishCurrentKeyframe();
			changeKeyframe(false, true, 1.0f);




			ret = true;
		}
		else
		{

			bool didSomething = updateKeyframe();




			ret = didSomething;
		}

		return ret;
	}
	else
	{
		// invalidate map if it was valid.
		if(map->isValid())
		{
			if(currentKeyFrame->numMappedOnThisTotal >= MIN_NUM_MAPPED)
				finishCurrentKeyframe();
			else
				discardCurrentKeyframe();

			map->invalidate();
		}

        std::cout << "Try the relocalizer .." << std::endl;
		//// start relocalizer if it isnt running already
		relocalizer.run(keyFrameGraph->keyframesAll);
		// did we find a frame to relocalize with?
		if(relocalizer.isGood()) {
	        std::cout << "The relocalizer SAVE US !! " << std::endl;
			takeRelocalizeResult();
		}

 		return true;
	}
}


void SlamSystem::gtDepthInit(unsigned char* image, float* depth, double timeStamp, int id)
{
	printf("Doing GT initialization!\n");


	currentKeyFrame.reset(new Frame(id, width, height, K, timeStamp, image));
	currentKeyFrame->setDepthFromGroundTruth(depth);

	map->initializeFromGTDepth(currentKeyFrame.get());
	keyFrameGraph->addFrame(currentKeyFrame.get());



	if(doSlam)
	{

		keyFrameGraph->idToKeyFrame.insert(std::make_pair(currentKeyFrame->id(), currentKeyFrame));

	}
	if(continuousPCOutput && outputWrapper != 0) outputWrapper->publishKeyframe(currentKeyFrame.get());

	printf("Done GT initialization!\n");
}


void SlamSystem::randomInit(unsigned char* image, double timeStamp, int id)
{
	printf("Doing Random initialization!\n");

	if(!doMapping)
		printf("WARNING: mapping is disabled, but we just initialized... THIS WILL NOT WORK! Set doMapping to true.\n");


	currentKeyFrame.reset(new Frame(id, width, height, K, timeStamp, image));
	map->initializeRandomly(currentKeyFrame.get());
	keyFrameGraph->addFrame(currentKeyFrame.get());


	if(doSlam)
	{

		keyFrameGraph->idToKeyFrame.insert(std::make_pair(currentKeyFrame->id(), currentKeyFrame));

	}
	if(continuousPCOutput && outputWrapper != 0) outputWrapper->publishKeyframe(currentKeyFrame.get());





	printf("Done Random initialization!\n");

}

void SlamSystem::trackFrame(unsigned char* image, unsigned int frameID, bool , double timestamp)
{

	// Create new frame
	std::shared_ptr<Frame> trackingNewFrame(new Frame(frameID, width, height, K, timestamp, image));


	bool my_createNewKeyframe = createNewKeyFrame;	// pre-save here, to make decision afterwards.
	if(trackingReference->keyframe != currentKeyFrame.get() || currentKeyFrame->depthHasBeenUpdatedFlag)
	{
		trackingReference->importFrame(currentKeyFrame.get());
		currentKeyFrame->depthHasBeenUpdatedFlag = false;
		trackingReferenceFrameSharedPT = currentKeyFrame;
	}

	FramePoseStruct* trackingReferencePose = trackingReference->keyframe->pose;



	SE3 frameToReference_initialEstimate = se3FromSim3(
			trackingReferencePose->getCamToWorld().inverse() * keyFrameGraph->allFramePoses.back()->getCamToWorld());



	struct timeval tv_start, tv_end;
	gettimeofday(&tv_start, NULL);

	SE3 newRefToFrame_poseUpdate;
	newRefToFrame_poseUpdate = tracker->trackFrame(
				trackingReference,
				trackingNewFrame.get(),
				frameToReference_initialEstimate);



	gettimeofday(&tv_end, NULL);
	msTrackFrame = 0.9*msTrackFrame + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
	nTrackFrame++;

	tracking_lastResidual = tracker->lastResidual;
	tracking_lastUsage = tracker->pointUsage;
	tracking_lastGoodPerBad = tracker->lastGoodCount / (tracker->lastGoodCount + tracker->lastBadCount);
	tracking_lastGoodPerTotal = tracker->lastGoodCount / (trackingNewFrame->width(SE3TRACKING_MIN_LEVEL)*trackingNewFrame->height(SE3TRACKING_MIN_LEVEL));


	if(manualTrackingLossIndicated || tracker->diverged || (keyFrameGraph->keyframesAll.size() > INITIALIZATION_PHASE_COUNT && !tracker->trackingWasGood))
	{
		printf("TRACKING LOST for frame %d (%1.2f%% good Points, which is %1.2f%% of available points, %s)!\n",
				trackingNewFrame->id(),
				100*tracking_lastGoodPerTotal,
				100*tracking_lastGoodPerBad,
				tracker->diverged ? "DIVERGED" : "NOT DIVERGED");

        trackingIsGood = false;

	    // I guess this is the case when we're lost !
	    if(!trackingIsGood)  {

	        std::cout << "We're gonna crash .." << std::endl;

	        relocalizer.updateCurrentFrame(trackingNewFrame);

            std::cout << "We're gonna crash for real .." << std::endl;

	        doMappingIteration() ;

            std::cout << "We crashed .." << std::endl;

	    }

        //// Stop the tracking if uncommented
		//trackingReference->invalidate();
		//trackingIsGood = true;
		//nextRelocIdx = -1;
		//manualTrackingLossIndicated = false;
        //std::cout << "We crashed for real.." << std::endl;

		//Might as well die, since retracking does not work.
		return;
	}




	keyFrameGraph->addFrame(trackingNewFrame.get());


	//Sim3 lastTrackedCamToWorld = mostCurrentTrackedFrame->getScaledCamToWorld();//  mostCurrentTrackedFrame->TrackingParent->getScaledCamToWorld() * sim3FromSE3(mostCurrentTrackedFrame->thisToParent_SE3TrackingResult, 1.0);
	if (outputWrapper != 0)
	{
		outputWrapper->publishTrackedFrame(trackingNewFrame.get());
	}


	// Keyframe selection
	latestTrackedFrame = trackingNewFrame;
	if (!my_createNewKeyframe && currentKeyFrame->numMappedOnThisTotal > MIN_NUM_MAPPED)
	{
		Sophus::Vector3d dist = newRefToFrame_poseUpdate.translation() * currentKeyFrame->meanIdepth;

		float minVal = fmin(0.2f + keyFrameGraph->keyframesAll.size() * 0.8f / INITIALIZATION_PHASE_COUNT,1.0f);

		if(keyFrameGraph->keyframesAll.size() < INITIALIZATION_PHASE_COUNT)	minVal *= 0.7;

		lastTrackingClosenessScore = trackableKeyFrameSearch->getRefFrameScore(dist.dot(dist), tracker->pointUsage);

		if (lastTrackingClosenessScore > minVal)
		{
			createNewKeyFrame = true;
		}
	}

	if(unmappedTrackedFrames.size() < 50 || (unmappedTrackedFrames.size() < 100 && trackingNewFrame->getTrackingParent()->numMappedOnThisTotal < 10))
		unmappedTrackedFrames.push_back(trackingNewFrame);



    doMappingIteration();



}


float SlamSystem::tryTrackSim3(
		TrackingReference* A, TrackingReference* B,
		int lvlStart, int lvlEnd,
		Sim3 &AtoB, Sim3 &BtoA,
		KFConstraintStruct* e1, KFConstraintStruct* e2 )
{
	BtoA = constraintTracker->trackFrameSim3(
			A,
			B->keyframe,
			BtoA,
			lvlStart,lvlEnd);
	Matrix7x7 BtoAInfo = constraintTracker->lastSim3Hessian;
	float BtoA_meanResidual = constraintTracker->lastResidual;
	float BtoA_meanDResidual = constraintTracker->lastDepthResidual;
	float BtoA_meanPResidual = constraintTracker->lastPhotometricResidual;
	float BtoA_usage = constraintTracker->pointUsage;


	if (constraintTracker->diverged ||
		BtoA.scale() > 1 / Sophus::SophusConstants<sophusType>::epsilon() ||
		BtoA.scale() < Sophus::SophusConstants<sophusType>::epsilon() ||
		BtoAInfo(0,0) == 0 ||
		BtoAInfo(6,6) == 0)
	{
		return 1e20;
	}


	AtoB = constraintTracker->trackFrameSim3(
			B,
			A->keyframe,
			AtoB,
			lvlStart,lvlEnd);
	Matrix7x7 AtoBInfo = constraintTracker->lastSim3Hessian;
	float AtoB_meanResidual = constraintTracker->lastResidual;
	float AtoB_meanDResidual = constraintTracker->lastDepthResidual;
	float AtoB_meanPResidual = constraintTracker->lastPhotometricResidual;
	float AtoB_usage = constraintTracker->pointUsage;


	if (constraintTracker->diverged ||
		AtoB.scale() > 1 / Sophus::SophusConstants<sophusType>::epsilon() ||
		AtoB.scale() < Sophus::SophusConstants<sophusType>::epsilon() ||
		AtoBInfo(0,0) == 0 ||
		AtoBInfo(6,6) == 0)
	{
		return 1e20;
	}

	// Propagate uncertainty (with d(a * b) / d(b) = Adj_a) and calculate Mahalanobis norm
	Matrix7x7 datimesb_db = AtoB.cast<float>().Adj();
	Matrix7x7 diffHesse = (AtoBInfo.inverse() + datimesb_db * BtoAInfo.inverse() * datimesb_db.transpose()).inverse();
	Vector7 diff = (AtoB * BtoA).log().cast<float>();


	float reciprocalConsistency = (diffHesse * diff).dot(diff);


	if(e1 != 0 && e2 != 0)
	{
		e1->firstFrame = A->keyframe;
		e1->secondFrame = B->keyframe;
		e1->secondToFirst = BtoA;
		e1->information = BtoAInfo.cast<double>();
		e1->meanResidual = BtoA_meanResidual;
		e1->meanResidualD = BtoA_meanDResidual;
		e1->meanResidualP = BtoA_meanPResidual;
		e1->usage = BtoA_usage;

		e2->firstFrame = B->keyframe;
		e2->secondFrame = A->keyframe;
		e2->secondToFirst = AtoB;
		e2->information = AtoBInfo.cast<double>();
		e2->meanResidual = AtoB_meanResidual;
		e2->meanResidualD = AtoB_meanDResidual;
		e2->meanResidualP = AtoB_meanPResidual;
		e2->usage = AtoB_usage;

		e1->reciprocalConsistency = e2->reciprocalConsistency = reciprocalConsistency;
	}

	return reciprocalConsistency;
}


void SlamSystem::testConstraint(
		Frame* candidate,
		KFConstraintStruct* &e1_out, KFConstraintStruct* &e2_out,
		Sim3 candidateToFrame_initialEstimate,
		float strictness)
{
	candidateTrackingReference->importFrame(candidate);

	Sim3 FtoC = candidateToFrame_initialEstimate.inverse(), CtoF = candidateToFrame_initialEstimate;
	Matrix7x7 FtoCInfo, CtoFInfo;

	float err_level3 = tryTrackSim3(
			newKFTrackingReference, candidateTrackingReference,	// A = frame; b = candidate
			SIM3TRACKING_MAX_LEVEL-1, 3,
			FtoC, CtoF);

	if(err_level3 > 3000*strictness)
	{

		e1_out = e2_out = 0;

		newKFTrackingReference->keyframe->trackingFailed.insert(std::pair<Frame*,Sim3>(candidate, candidateToFrame_initialEstimate));
		return;
	}

	float err_level2 = tryTrackSim3(
			newKFTrackingReference, candidateTrackingReference,	// A = frame; b = candidate
			2, 2,
			FtoC, CtoF);

	if(err_level2 > 4000*strictness)
	{

		e1_out = e2_out = 0;
		newKFTrackingReference->keyframe->trackingFailed.insert(std::pair<Frame*,Sim3>(candidate, candidateToFrame_initialEstimate));
		return;
	}

	e1_out = new KFConstraintStruct();
	e2_out = new KFConstraintStruct();


	float err_level1 = tryTrackSim3(
			newKFTrackingReference, candidateTrackingReference,	// A = frame; b = candidate
			1, 1,
			FtoC, CtoF, e1_out, e2_out);

	if(err_level1 > 6000*strictness)
	{

		delete e1_out;
		delete e2_out;
		e1_out = e2_out = 0;
		newKFTrackingReference->keyframe->trackingFailed.insert(std::pair<Frame*,Sim3>(candidate, candidateToFrame_initialEstimate));
		return;
	}




	const float kernelDelta = 5 * sqrt(6000*loopclosureStrictness);
	e1_out->robustKernel = new g2o::RobustKernelHuber();
	e1_out->robustKernel->setDelta(kernelDelta);
	e2_out->robustKernel = new g2o::RobustKernelHuber();
	e2_out->robustKernel->setDelta(kernelDelta);
}

int SlamSystem::findConstraintsForNewKeyFrames(Frame* newKeyFrame, bool forceParent, float closeCandidatesTH)
{
	if(!newKeyFrame->hasTrackingParent())
	{

		keyFrameGraph->addKeyFrame(newKeyFrame);

	    if(doFinalOptimization)
	    {
	            printf("doing final optimization iteration!\n");
	            optimizationIteration(50, 0.001);
	            doFinalOptimization = false;
	    }

	    while(optimizationIteration(5, 0.02));
        {
	           newContraintsOptimised = true;
	    }



		return 0;
	}

	if(!forceParent && (newKeyFrame->lastConstraintTrackedCamToWorld * newKeyFrame->getScaledCamToWorld().inverse()).log().norm() < 0.01)
		return 0;


	newKeyFrame->lastConstraintTrackedCamToWorld = newKeyFrame->getScaledCamToWorld();

	// =============== get all potential candidates and their initial relative pose. =================
	std::vector<KFConstraintStruct*> constraints;
	Frame* fabMapResult = 0;


	std::unordered_set<Frame*> candidates = trackableKeyFrameSearch->findCandidates(newKeyFrame, fabMapResult, closeCandidatesTH);



	std::map< Frame*, Sim3 > candidateToFrame_initialEstimateMap;


	// erase the ones that are already neighbours.
	for(std::unordered_set<Frame*>::iterator c = candidates.begin(); c != candidates.end();)
	{
		if(newKeyFrame->neighbors.find(*c) != newKeyFrame->neighbors.end())
		{
			c = candidates.erase(c);
		}
		else {
			++c;
		}
	}

	for (Frame* candidate : candidates)
	{
		Sim3 candidateToFrame_initialEstimate = newKeyFrame->getScaledCamToWorld().inverse() * candidate->getScaledCamToWorld();
		candidateToFrame_initialEstimateMap[candidate] = candidateToFrame_initialEstimate;
	}

	std::unordered_map<Frame*, int> distancesToNewKeyFrame;
	if(newKeyFrame->hasTrackingParent())
		keyFrameGraph->calculateGraphDistancesToFrame(newKeyFrame->getTrackingParent(), &distancesToNewKeyFrame);






	// =============== distinguish between close and "far" candidates in Graph =================
	// Do a first check on trackability of close candidates.
	std::unordered_set<Frame*> closeCandidates;
	std::vector<Frame*> farCandidates;
	Frame* parent = newKeyFrame->hasTrackingParent() ? newKeyFrame->getTrackingParent() : 0;

	int closeFailed = 0;
	int closeInconsistent = 0;

	SO3 disturbance = SO3::exp(Sophus::Vector3d(0.05,0,0));

	for (Frame* candidate : candidates)
	{
		if (candidate->id() == newKeyFrame->id())
			continue;
		if(!candidate->pose->isInGraph)
			continue;
		if(newKeyFrame->hasTrackingParent() && candidate == newKeyFrame->getTrackingParent())
			continue;
		if(candidate->idxInKeyframes < INITIALIZATION_PHASE_COUNT)
			continue;

		SE3 c2f_init = se3FromSim3(candidateToFrame_initialEstimateMap[candidate].inverse()).inverse();
		c2f_init.so3() = c2f_init.so3() * disturbance;
		SE3 c2f = constraintSE3Tracker->trackFrameOnPermaref(candidate, newKeyFrame, c2f_init);
		if(!constraintSE3Tracker->trackingWasGood) {closeFailed++; continue;}


		SE3 f2c_init = se3FromSim3(candidateToFrame_initialEstimateMap[candidate]).inverse();
		f2c_init.so3() = disturbance * f2c_init.so3();
		SE3 f2c = constraintSE3Tracker->trackFrameOnPermaref(newKeyFrame, candidate, f2c_init);
		if(!constraintSE3Tracker->trackingWasGood) {closeFailed++; continue;}

		if((f2c.so3() * c2f.so3()).log().norm() >= 0.09) {closeInconsistent++; continue;}

		closeCandidates.insert(candidate);
	}



	for (Frame* candidate : candidates)
	{
		if (candidate->id() == newKeyFrame->id())
			continue;
		if(!candidate->pose->isInGraph)
			continue;
		if(newKeyFrame->hasTrackingParent() && candidate == newKeyFrame->getTrackingParent())
			continue;
		if(candidate->idxInKeyframes < INITIALIZATION_PHASE_COUNT)
			continue;

		if(candidate == fabMapResult)
		{
			farCandidates.push_back(candidate);
			continue;
		}

		if(distancesToNewKeyFrame.at(candidate) < 4)
			continue;

		farCandidates.push_back(candidate);
	}





	// erase the ones that we tried already before (close)
	for(std::unordered_set<Frame*>::iterator c = closeCandidates.begin(); c != closeCandidates.end();)
	{
		if(newKeyFrame->trackingFailed.find(*c) == newKeyFrame->trackingFailed.end())
		{
			++c;
			continue;
		}
		auto range = newKeyFrame->trackingFailed.equal_range(*c);

		bool skip = false;
		Sim3 f2c = candidateToFrame_initialEstimateMap[*c].inverse();
		for (auto it = range.first; it != range.second; ++it)
		{
			if((f2c * it->second).log().norm() < 0.1)
			{
				skip=true;
				break;
			}
		}

		if(skip)
		{
			c = closeCandidates.erase(c);
		}
		else {
			++c;
		}
	}

	// erase the ones that are already neighbours (far)
	for(unsigned int i=0;i<farCandidates.size();i++)
	{
		if(newKeyFrame->trackingFailed.find(farCandidates[i]) == newKeyFrame->trackingFailed.end())
			continue;

		auto range = newKeyFrame->trackingFailed.equal_range(farCandidates[i]);

		bool skip = false;
		for (auto it = range.first; it != range.second; ++it)
		{
			if((it->second).log().norm() < 0.2)
			{
				skip=true;
				break;
			}
		}

		if(skip)
		{
			farCandidates[i] = farCandidates.back();
			farCandidates.pop_back();
			i--;
		}
	}

	// =============== limit number of close candidates ===============
	// while too many, remove the one with the highest connectivity.
	while((int)closeCandidates.size() > maxLoopClosureCandidates)
	{
		Frame* worst = 0;
		int worstNeighbours = 0;
		for(Frame* f : closeCandidates)
		{
			int neightboursInCandidates = 0;
			for(Frame* n : f->neighbors)
				if(closeCandidates.find(n) != closeCandidates.end())
					neightboursInCandidates++;

			if(neightboursInCandidates > worstNeighbours || worst == 0)
			{
				worst = f;
				worstNeighbours = neightboursInCandidates;
			}
		}

		closeCandidates.erase(worst);
	}


	// =============== limit number of far candidates ===============
	// delete randomly
	int maxNumFarCandidates = (maxLoopClosureCandidates +1) / 2;
	if(maxNumFarCandidates < 5) maxNumFarCandidates = 5;
	while((int)farCandidates.size() > maxNumFarCandidates)
	{
		int toDelete = rand() % farCandidates.size();
		if(farCandidates[toDelete] != fabMapResult)
		{
			farCandidates[toDelete] = farCandidates.back();
			farCandidates.pop_back();
		}
	}



	// =============== TRACK! ===============

	// make tracking reference for newKeyFrame.
	newKFTrackingReference->importFrame(newKeyFrame);


	for (Frame* candidate : closeCandidates)
	{
		KFConstraintStruct* e1=0;
		KFConstraintStruct* e2=0;

		testConstraint(
				candidate, e1, e2,
				candidateToFrame_initialEstimateMap[candidate],
				loopclosureStrictness);


		if(e1 != 0)
		{
			constraints.push_back(e1);
			constraints.push_back(e2);

			// delete from far candidates if it's in there.
			for(unsigned int k=0;k<farCandidates.size();k++)
			{
				if(farCandidates[k] == candidate)
				{

					farCandidates[k] = farCandidates.back();
					farCandidates.pop_back();
				}
			}
		}
	}


	for (Frame* candidate : farCandidates)
	{
		KFConstraintStruct* e1=0;
		KFConstraintStruct* e2=0;

		testConstraint(
				candidate, e1, e2,
				Sim3(),
				loopclosureStrictness);


		if(e1 != 0)
		{
			constraints.push_back(e1);
			constraints.push_back(e2);
		}
	}



	if(parent != 0 && forceParent)
	{
		KFConstraintStruct* e1=0;
		KFConstraintStruct* e2=0;
		testConstraint(
				parent, e1, e2,
				candidateToFrame_initialEstimateMap[parent],
				100);

		if(e1 != 0)
		{
			constraints.push_back(e1);
			constraints.push_back(e2);
		}
		else
		{
			float downweightFac = 5;
			const float kernelDelta = 5 * sqrt(6000*loopclosureStrictness) / downweightFac;
			printf("warning: reciprocal tracking on new frame failed badly, added odometry edge (Hacky).\n");

			constraints.push_back(new KFConstraintStruct());
			constraints.back()->firstFrame = newKeyFrame;
			constraints.back()->secondFrame = newKeyFrame->getTrackingParent();
			constraints.back()->secondToFirst = constraints.back()->firstFrame->getScaledCamToWorld().inverse() * constraints.back()->secondFrame->getScaledCamToWorld();
			constraints.back()->information  <<
					0.8098,-0.1507,-0.0557, 0.1211, 0.7657, 0.0120, 0,
					-0.1507, 2.1724,-0.1103,-1.9279,-0.1182, 0.1943, 0,
					-0.0557,-0.1103, 0.2643,-0.0021,-0.0657,-0.0028, 0.0304,
					 0.1211,-1.9279,-0.0021, 2.3110, 0.1039,-0.0934, 0.0005,
					 0.7657,-0.1182,-0.0657, 0.1039, 1.0545, 0.0743,-0.0028,
					 0.0120, 0.1943,-0.0028,-0.0934, 0.0743, 0.4511, 0,
					0,0, 0.0304, 0.0005,-0.0028, 0, 0.0228;
			constraints.back()->information *= (1e9/(downweightFac*downweightFac));

			constraints.back()->robustKernel = new g2o::RobustKernelHuber();
			constraints.back()->robustKernel->setDelta(kernelDelta);

			constraints.back()->meanResidual = 10;
			constraints.back()->meanResidualD = 10;
			constraints.back()->meanResidualP = 10;
			constraints.back()->usage = 0;

		}
	}



	keyFrameGraph->addKeyFrame(newKeyFrame);
	for(unsigned int i=0;i<constraints.size();i++)
		keyFrameGraph->insertConstraint(constraints[i]);



	newKFTrackingReference->invalidate();
	candidateTrackingReference->invalidate();

	return constraints.size();
}




bool SlamSystem::optimizationIteration(int itsPerTry, float minChange)
{
	struct timeval tv_start, tv_end;
	gettimeofday(&tv_start, NULL);


	// lock new elements buffer & take them over.

	keyFrameGraph->addElementsFromBuffer();



	// Do the optimization. This can take quite some time!
	int its = keyFrameGraph->optimize(itsPerTry);

	// save the optimization result.

	float maxChange = 0;
	float sumChange = 0;
	float sum = 0;
	for(size_t i=0;i<keyFrameGraph->keyframesAll.size(); i++)
	{
		// set edge error sum to zero
		keyFrameGraph->keyframesAll[i]->edgeErrorSum = 0;
		keyFrameGraph->keyframesAll[i]->edgesNum = 0;

		if(!keyFrameGraph->keyframesAll[i]->pose->isInGraph) continue;



		// get change from last optimization
		Sim3 a = keyFrameGraph->keyframesAll[i]->pose->graphVertex->estimate();
		Sim3 b = keyFrameGraph->keyframesAll[i]->getScaledCamToWorld();
		Sophus::Vector7f diff = (a*b.inverse()).log().cast<float>();


		for(int j=0;j<7;j++)
		{
			float d = fabsf((float)(diff[j]));
			if(d > maxChange) maxChange = d;
			sumChange += d;
		}
		sum +=7;

		// set change
		keyFrameGraph->keyframesAll[i]->pose->setPoseGraphOptResult(
				keyFrameGraph->keyframesAll[i]->pose->graphVertex->estimate());

		// add error
		for(auto edge : keyFrameGraph->keyframesAll[i]->pose->graphVertex->edges())
		{
			keyFrameGraph->keyframesAll[i]->edgeErrorSum += ((EdgeSim3*)(edge))->chi2();
			keyFrameGraph->keyframesAll[i]->edgesNum++;
		}
	}

	haveUnmergedOptimizationOffset = true;



	gettimeofday(&tv_end, NULL);
	msOptimizationIteration = 0.9*msOptimizationIteration + 0.1*((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
	nOptimizationIteration++;


	return maxChange > minChange && its == itsPerTry;
}

void SlamSystem::optimizeGraph()
{
	keyFrameGraph->optimize(1000);

	mergeOptimizationOffset();
}


SE3 SlamSystem::getCurrentPoseEstimate()
{
	SE3 camToWorld = SE3();
	if(keyFrameGraph->allFramePoses.size() > 0)
		camToWorld = se3FromSim3(keyFrameGraph->allFramePoses.back()->getCamToWorld());

	return camToWorld;
}

Sophus::Sim3f SlamSystem::getCurrentPoseEstimateScale()
{
    Sophus::Sim3f camToWorld = Sophus::Sim3f();

    if(keyFrameGraph->allFramePoses.size() > 0)
        camToWorld = keyFrameGraph->allFramePoses.back()->getCamToWorld().cast<float>();

    return camToWorld;
}

std::vector<FramePoseStruct*> SlamSystem::getAllPoses()
{
	return keyFrameGraph->allFramePoses;
}

bool SlamSystem::relocaliserHasRun() {
	return !map->isValid();
}
