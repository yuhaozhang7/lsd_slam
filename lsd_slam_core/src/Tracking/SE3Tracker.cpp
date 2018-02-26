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

#include "SE3Tracker.h"


#include "DataStructures/Frame.h"
#include "Tracking/TrackingReference.h"
#include "util/globalFuncs.h"
#include "Tracking/least_squares.h"

namespace lsd_slam
{


#define callOptimized(function, arguments) function arguments



SE3Tracker::SE3Tracker(int w, int h, Eigen::Matrix3f K)
{
	width = w;
	height = h;

	this->K = K;
	fx = K(0,0);
	fy = K(1,1);
	cx = K(0,2);
	cy = K(1,2);

	settings = DenseDepthTrackerSettings();
	//settings.maxItsPerLvl[0] = 2;

	KInv = K.inverse();
	fxi = KInv(0,0);
	fyi = KInv(1,1);
	cxi = KInv(0,2);
	cyi = KInv(1,2);


	buf_warped_residual = new float[w*h];
	buf_warped_dx = new float[w*h];
	buf_warped_dy = new float[w*h];
	buf_warped_x = new float[w*h];
	buf_warped_y = new float[w*h];
	buf_warped_z = new float[w*h];

	buf_d = new float[w*h];
	buf_idepthVar = new float[w*h];
	buf_weight_p = new float[w*h];

	buf_warped_size = 0;



	lastResidual = 0;
	iterationNumber = 0;
	pointUsage = 0;
	lastGoodCount = lastBadCount = 0;

	diverged = false;
}

SE3Tracker::~SE3Tracker()
{


	delete[] buf_warped_residual;
	delete[] buf_warped_dx;
	delete[] buf_warped_dy;
	delete[] buf_warped_x;
	delete[] buf_warped_y;
	delete[] buf_warped_z;

	delete[] buf_d;
	delete[] buf_idepthVar;
	delete[] buf_weight_p;
}



// tracks a frame.
// first_frame has depth, second_frame DOES NOT have depth.
float SE3Tracker::checkPermaRefOverlap(
		Frame* reference,
		SE3 referenceToFrameOrg)
{
	Sophus::SE3f referenceToFrame = referenceToFrameOrg.cast<float>();

	int w2 = reference->width(QUICK_KF_CHECK_LVL)-1;
	int h2 = reference->height(QUICK_KF_CHECK_LVL)-1;
	Eigen::Matrix3f KLvl = reference->K(QUICK_KF_CHECK_LVL);
	float fx_l = KLvl(0,0);
	float fy_l = KLvl(1,1);
	float cx_l = KLvl(0,2);
	float cy_l = KLvl(1,2);

	Eigen::Matrix3f rotMat = referenceToFrame.rotationMatrix();
	Eigen::Vector3f transVec = referenceToFrame.translation();

	const Eigen::Vector3f* refPoint_max = reference->permaRef_posData + reference->permaRefNumPts;
	const Eigen::Vector3f* refPoint = reference->permaRef_posData;

	float usageCount = 0;
	for(;refPoint<refPoint_max; refPoint++)
	{
		Eigen::Vector3f Wxp = rotMat * (*refPoint) + transVec;
		float u_new = (Wxp[0]/Wxp[2])*fx_l + cx_l;
		float v_new = (Wxp[1]/Wxp[2])*fy_l + cy_l;
		if((u_new > 0 && v_new > 0 && u_new < w2 && v_new < h2))
		{
			float depthChange = (*refPoint)[2] / Wxp[2];
			usageCount += depthChange < 1 ? depthChange : 1;
		}
	}

	pointUsage = usageCount / (float)reference->permaRefNumPts;
	return pointUsage;
}


// tracks a frame.
// first_frame has depth, second_frame DOES NOT have depth.
SE3 SE3Tracker::trackFrameOnPermaref(
		Frame* reference,
		Frame* frame,
		SE3 referenceToFrameOrg)
{

	Sophus::SE3f referenceToFrame = referenceToFrameOrg.cast<float>();


	affineEstimation_a = 1; affineEstimation_b = 0;

	NormalEquationsLeastSquares ls;
	diverged = false;
	trackingWasGood = true;

	{

			callOptimized(calcResidualAndBuffers, (reference->permaRef_posData, reference->permaRef_colorAndVarData, 0, reference->permaRefNumPts, frame, referenceToFrame, QUICK_KF_CHECK_LVL, false));

	}

	if(buf_warped_size < MIN_GOODPERALL_PIXEL_ABSMIN * (width>>QUICK_KF_CHECK_LVL)*(height>>QUICK_KF_CHECK_LVL))
	{
		diverged = true;
		trackingWasGood = false;
		return SE3();
	}
	if(useAffineLightningEstimation)
	{
		affineEstimation_a = affineEstimation_a_lastIt;
		affineEstimation_b = affineEstimation_b_lastIt;
	}

	float lastErr;
	{

			lastErr =callOptimized(calcWeightsAndResidual, (referenceToFrame));

	}

	float LM_lambda = settings.lambdaInitialTestTrack;

	for(int iteration=0; iteration < settings.maxItsTestTrack; iteration++)
	{
		{

			callOptimized(calculateWarpUpdate, (ls));

		}


		int incTry=0;
		while(true)
		{
			// solve LS system with current lambda
			Vector6 b = -ls.b;
			Matrix6x6 A = ls.A;
			for(int i=0;i<6;i++) A(i,i) *= 1+LM_lambda;
			Vector6 inc = A.ldlt().solve(b);
			incTry++;

			// apply increment. pretty sure this way round is correct, but hard to test.
			Sophus::SE3f new_referenceToFrame = Sophus::SE3f::exp((inc)) * referenceToFrame;

			// re-evaluate residual
			{

					callOptimized(calcResidualAndBuffers, (reference->permaRef_posData, reference->permaRef_colorAndVarData, 0, reference->permaRefNumPts, frame, new_referenceToFrame, QUICK_KF_CHECK_LVL, false));
			}


			if(buf_warped_size < MIN_GOODPERALL_PIXEL_ABSMIN * (width>>QUICK_KF_CHECK_LVL)*(height>>QUICK_KF_CHECK_LVL))
			{
				diverged = true;
				trackingWasGood = false;
				return SE3();
			}

			float error;
			{

				error = callOptimized(calcWeightsAndResidual,(new_referenceToFrame));
			}

			// accept inc?
			if(error < lastErr)
			{
				// accept inc
				referenceToFrame = new_referenceToFrame;
				if(useAffineLightningEstimation)
				{
					affineEstimation_a = affineEstimation_a_lastIt;
					affineEstimation_b = affineEstimation_b_lastIt;
				}
				// converged?
				if(error / lastErr > settings.convergenceEpsTestTrack)
					iteration = settings.maxItsTestTrack;


				lastErr = error;


				if(LM_lambda <= 0.2)
					LM_lambda = 0;
				else
					LM_lambda *= settings.lambdaSuccessFac;

				break;
			}
			else
			{
				if(!(inc.dot(inc) > settings.stepSizeMinTestTrack))
				{
					iteration = settings.maxItsTestTrack;
					break;
				}

				if(LM_lambda == 0)
					LM_lambda = 0.2;
				else
					LM_lambda *= std::pow(settings.lambdaFailFac, incTry);
			}
		}
	}

	lastResidual = lastErr;

	trackingWasGood = !diverged
			&& lastGoodCount / (frame->width(QUICK_KF_CHECK_LVL)*frame->height(QUICK_KF_CHECK_LVL)) > MIN_GOODPERALL_PIXEL
			&& lastGoodCount / (lastGoodCount + lastBadCount) > MIN_GOODPERGOODBAD_PIXEL;

	return toSophus(referenceToFrame);
}





// tracks a frame.
// first_frame has depth, second_frame DOES NOT have depth.
SE3 SE3Tracker::trackFrame(
		TrackingReference* reference,
		Frame* frame,
		const SE3& frameToReference_initialEstimate)
{


	diverged = false;
	trackingWasGood = true;
	affineEstimation_a = 1; affineEstimation_b = 0;

	if(saveAllTrackingStages)
	{
		saveAllTrackingStages = false;
		saveAllTrackingStagesInternal = true;
	}
	

	// ============ track frame ============
	Sophus::SE3f referenceToFrame = frameToReference_initialEstimate.inverse().cast<float>();
	NormalEquationsLeastSquares ls;


	int numCalcResidualCalls[PYRAMID_LEVELS];
	int numCalcWarpUpdateCalls[PYRAMID_LEVELS];

	float last_residual = 0;


	for(int lvl=SE3TRACKING_MAX_LEVEL-1;lvl >= SE3TRACKING_MIN_LEVEL;lvl--)
	{
		numCalcResidualCalls[lvl] = 0;
		numCalcWarpUpdateCalls[lvl] = 0;

		reference->makePointCloud(lvl);

		{

			callOptimized(calcResidualAndBuffers,
						  (reference->posData[lvl], reference->colorAndVarData[lvl], SE3TRACKING_MIN_LEVEL == lvl
																					 ? reference->pointPosInXYGrid[lvl]
																					 : 0, reference->numData[lvl], frame, referenceToFrame, lvl, false ));
		}

		if(buf_warped_size < MIN_GOODPERALL_PIXEL_ABSMIN * (width>>lvl)*(height>>lvl))
		{
			diverged = true;
			trackingWasGood = false;
			return SE3();
		}

		if(useAffineLightningEstimation)
		{
			affineEstimation_a = affineEstimation_a_lastIt;
			affineEstimation_b = affineEstimation_b_lastIt;
		}

		float lastErr;
		{

			lastErr = callOptimized(calcWeightsAndResidual, (referenceToFrame));

		}

		numCalcResidualCalls[lvl]++;


		float LM_lambda = settings.lambdaInitial[lvl];

		for(int iteration=0; iteration < settings.maxItsPerLvl[lvl]; iteration++)
		{

			{
				callOptimized(calculateWarpUpdate,(ls));

			}

			numCalcWarpUpdateCalls[lvl]++;

			iterationNumber = iteration;

			int incTry=0;
			while(true)
			{
				// solve LS system with current lambda
				Vector6 b = -ls.b;
				Matrix6x6 A = ls.A;
				for(int i=0;i<6;i++) A(i,i) *= 1+LM_lambda;
				Vector6 inc = A.ldlt().solve(b);
				incTry++;


				// apply increment. pretty sure this way round is correct, but hard to test.
				Sophus::SE3f new_referenceToFrame = Sophus::SE3f::exp((inc)) * referenceToFrame;
				//Sophus::SE3f new_referenceToFrame = referenceToFrame * Sophus::SE3f::exp((inc));


				{
					// re-evaluate residual
					callOptimized(calcResidualAndBuffers,
								  (reference->posData[lvl], reference->colorAndVarData[lvl], SE3TRACKING_MIN_LEVEL ==
																							 lvl
																							 ? reference->pointPosInXYGrid[lvl]
																							 : 0, reference->numData[lvl], frame, new_referenceToFrame, lvl, (false)));

				}

				if(buf_warped_size < MIN_GOODPERALL_PIXEL_ABSMIN* (width>>lvl)*(height>>lvl))
				{
					diverged = true;
					trackingWasGood = false;
					return SE3();
				}

				float error;
				{

					error = callOptimized(calcWeightsAndResidual,(new_referenceToFrame));

					}

				numCalcResidualCalls[lvl]++;


				// accept inc?
				if(error < lastErr)
				{
					// accept inc
					referenceToFrame = new_referenceToFrame;
					if(useAffineLightningEstimation)
					{
						affineEstimation_a = affineEstimation_a_lastIt;
						affineEstimation_b = affineEstimation_b_lastIt;
					}




					// converged?
					if(error / lastErr > settings.convergenceEps[lvl])
					{

						iteration = settings.maxItsPerLvl[lvl];
					}

					last_residual = lastErr = error;


					if(LM_lambda <= 0.2)
						LM_lambda = 0;
					else
						LM_lambda *= settings.lambdaSuccessFac;

					break;
				}
				else
				{


					if(!(inc.dot(inc) > settings.stepSizeMin[lvl]))
					{

						iteration = settings.maxItsPerLvl[lvl];
						break;
					}

					if(LM_lambda == 0)
						LM_lambda = 0.2;
					else
						LM_lambda *= std::pow(settings.lambdaFailFac, incTry);
				}
			}
		}
	}






	saveAllTrackingStagesInternal = false;

	lastResidual = last_residual;

	trackingWasGood = !diverged
			&& lastGoodCount / (frame->width(SE3TRACKING_MIN_LEVEL)*frame->height(SE3TRACKING_MIN_LEVEL)) > MIN_GOODPERALL_PIXEL
			&& lastGoodCount / (lastGoodCount + lastBadCount) > MIN_GOODPERGOODBAD_PIXEL;

	if(trackingWasGood)
		reference->keyframe->numFramesTrackedOnThis++;

	frame->initialTrackedResidual = lastResidual / pointUsage;
	frame->pose->thisToParent_raw = sim3FromSE3(toSophus(referenceToFrame.inverse()),1);
	frame->pose->trackingParent = reference->keyframe->pose;
	return toSophus(referenceToFrame.inverse());
}




float SE3Tracker::calcWeightsAndResidual(
		const Sophus::SE3f& referenceToFrame)
{
	float tx = referenceToFrame.translation()[0];
	float ty = referenceToFrame.translation()[1];
	float tz = referenceToFrame.translation()[2];

	float sumRes = 0;

	for(int i=0;i<buf_warped_size;i++)
	{
		float px = *(buf_warped_x+i);	// x'
		float py = *(buf_warped_y+i);	// y'
		float pz = *(buf_warped_z+i);	// z'
		float d = *(buf_d+i);	// d
		float rp = *(buf_warped_residual+i); // r_p
		float gx = *(buf_warped_dx+i);	// \delta_x I
		float gy = *(buf_warped_dy+i);  // \delta_y I
		float s = settings.var_weight * *(buf_idepthVar+i);	// \sigma_d^2


		// calc dw/dd (first 2 components):
		float g0 = (tx * pz - tz * px) / (pz*pz*d);
		float g1 = (ty * pz - tz * py) / (pz*pz*d);


		// calc w_p
		float drpdd = gx * g0 + gy * g1;	// ommitting the minus
		float w_p = 1.0f / ((cameraPixelNoise2) + s * drpdd * drpdd);

		float weighted_rp = fabs(rp*sqrtf(w_p));

		float wh = fabs(weighted_rp < (settings.huber_d/2) ? 1 : (settings.huber_d/2) / weighted_rp);

		sumRes += wh * w_p * rp*rp;


		*(buf_weight_p+i) = wh * w_p;
	}

	return sumRes / buf_warped_size;
}




float SE3Tracker::calcResidualAndBuffers(
		const Eigen::Vector3f* refPoint,
		const Eigen::Vector2f* refColVar,
		int* idxBuf,
		int refNum,
		Frame* frame,
		const Sophus::SE3f& referenceToFrame,
		int level,
		bool )
{



	int w = frame->width(level);
	int h = frame->height(level);
	Eigen::Matrix3f KLvl = frame->K(level);
	float fx_l = KLvl(0,0);
	float fy_l = KLvl(1,1);
	float cx_l = KLvl(0,2);
	float cy_l = KLvl(1,2);

	Eigen::Matrix3f rotMat = referenceToFrame.rotationMatrix();
	Eigen::Vector3f transVec = referenceToFrame.translation();
	
	const Eigen::Vector3f* refPoint_max = refPoint + refNum;


	const Eigen::Vector4f* frame_gradients = frame->gradients(level);

	int idx=0;

	float sumResUnweighted = 0;

	bool* isGoodOutBuffer = idxBuf != 0 ? frame->refPixelWasGood() : 0;

	int goodCount = 0;
	int badCount = 0;

	float sumSignedRes = 0;



	float sxx=0,syy=0,sx=0,sy=0,sw=0;

	float usageCount = 0;

	for(;refPoint<refPoint_max; refPoint++, refColVar++, idxBuf++)
	{

		Eigen::Vector3f Wxp = rotMat * (*refPoint) + transVec;
		float u_new = (Wxp[0]/Wxp[2])*fx_l + cx_l;
		float v_new = (Wxp[1]/Wxp[2])*fy_l + cy_l;

		// step 1a: coordinates have to be in image:
		// (inverse test to exclude NANs)
		if(!(u_new > 1 && v_new > 1 && u_new < w-2 && v_new < h-2))
		{
			if(isGoodOutBuffer != 0)
				isGoodOutBuffer[*idxBuf] = false;
			continue;
		}

		Eigen::Vector3f resInterp = getInterpolatedElement43(frame_gradients, u_new, v_new, w);

		float c1 = affineEstimation_a * (*refColVar)[0] + affineEstimation_b;
		float c2 = resInterp[2];
		float residual = c1 - c2;

		float weight = fabsf(residual) < 5.0f ? 1 : 5.0f / fabsf(residual);
		sxx += c1*c1*weight;
		syy += c2*c2*weight;
		sx += c1*weight;
		sy += c2*weight;
		sw += weight;

		bool isGood = residual*residual / (MAX_DIFF_CONSTANT + MAX_DIFF_GRAD_MULT*(resInterp[0]*resInterp[0] + resInterp[1]*resInterp[1])) < 1;

		if(isGoodOutBuffer != 0)
			isGoodOutBuffer[*idxBuf] = isGood;

		*(buf_warped_x+idx) = Wxp(0);
		*(buf_warped_y+idx) = Wxp(1);
		*(buf_warped_z+idx) = Wxp(2);

		*(buf_warped_dx+idx) = fx_l * resInterp[0];
		*(buf_warped_dy+idx) = fy_l * resInterp[1];
		*(buf_warped_residual+idx) = residual;

		*(buf_d+idx) = 1.0f / (*refPoint)[2];
		*(buf_idepthVar+idx) = (*refColVar)[1];
		idx++;


		if(isGood)
		{
			sumResUnweighted += residual*residual;
			sumSignedRes += residual;
			goodCount++;
		}
		else
			badCount++;

		float depthChange = (*refPoint)[2] / Wxp[2];	// if depth becomes larger: pixel becomes "smaller", hence count it less.
		usageCount += depthChange < 1 ? depthChange : 1;



	}

	buf_warped_size = idx;

	pointUsage = usageCount / (float)refNum;
	lastGoodCount = goodCount;
	lastBadCount = badCount;
	lastMeanRes = sumSignedRes / goodCount;

	affineEstimation_a_lastIt = sqrtf((syy - sy*sy/sw) / (sxx - sx*sx/sw));
	affineEstimation_b_lastIt = (sy - affineEstimation_a_lastIt*sx)/sw;


	return sumResUnweighted / goodCount;
}



Vector6 SE3Tracker::calculateWarpUpdate(
		NormalEquationsLeastSquares &ls) {
//	weightEstimator.reset();
//	weightEstimator.estimateDistribution(buf_warped_residual, buf_warped_size);
//	weightEstimator.calcWeights(buf_warped_residual, buf_warped_weights, buf_warped_size);
//

//	{
//	float px = *(buf_warped_x);
//	float py = *(buf_warped_y);
//	float pz = *(buf_warped_z);
//	float r = *(buf_warped_residual);
//	float gx = *(buf_warped_dx);
//	float gy = *(buf_warped_dy);
//
//		std::cout << "px: " << px << std::endl;
//		std::cout << "py: " << py << std::endl;
//		std::cout << "pz: " << pz << std::endl;
//		std::cout << "r: " << r << std::endl;
//		std::cout << "gx: " << gx << std::endl;
//		std::cout << "gy: " << gy << std::endl;
//
//	}


	ls.initialize(width*height);
	for(int i=0;i<buf_warped_size;i++)
	{
		float px = *(buf_warped_x+i);
		float py = *(buf_warped_y+i);
		float pz = *(buf_warped_z+i);
		float r =  *(buf_warped_residual+i);
		float gx = *(buf_warped_dx+i);
		float gy = *(buf_warped_dy+i);
		// step 3 + step 5 comp 6d error vector

		float z = 1.0f / pz;
		float z_sqr = z * z ; //1.0f / (pz * pz); //z * z; //1.0f / (pz*pz);
		Vector6 v;
//
//		if (fabs(z_sqr - (z*z)) > 0.00000001) {
//			std::cout.precision(17);
//			std::cout << std::fixed << fabs(z_sqr - (z * z)) << "pz: " << std::fixed << pz << ", z: " << std::fixed << z << ", z_sqr: " << std::fixed << z_sqr << ", z^2: " <<
//					std::fixed <<  z * z << std::endl;
//		}

		v[0] = z*gx + 0;
		v[1] = 0 +         z*gy;
		v[2] = (-px * z_sqr) * gx +
			  (-py * z_sqr) * gy;
		v[3] = (-px * py * z_sqr) * gx +
			  (-(1.0 + py * py * z_sqr)) * gy;
		v[4] = (1.0 + px * px * z_sqr) * gx +
			  (px * py * z_sqr) * gy;
		v[5] = (-py * z) * gx +
			  (px * z) * gy;

//		v[0] = z*gx;
//		v[1] = z*gy;
//		v[5] =  (px * gy) * z - (py * gx) * z;
//
//		v[2] = -((px * gx) * z_sqr +
//			     (py * gy) * z_sqr);
//
//		v[3] = -(gy + (((px * gx) * z_sqr) * py) +
//				 ((py * gy) * z_sqr) * py);
//
//		//	   (-(py * py * z_sqr)) * gy;
//
//
//		v[4] = (gx + ((px * gx) * z_sqr) * px) +
//				(((py * gy) * z_sqr) * px);

//
//		v[0] = z*gx;
//		v[1] = z*gy;
//		v[2] = (-px * z_sqr) * gx + (-py * z_sqr) * gy;
//
//		v[3] = (-px * py * z_sqr) * gx + (-(py * py * z_sqr)) * gy - gy;
//
//		v[4] = gx + (px * px * z_sqr) * gx + (px * py * z_sqr) * gy;
//
//		v[5] = (px * gy) * z - (py * gx) * z; //(-py * z) * gx + (px * z) * gy;



		// step 6: integrate into A and b:
		ls.update(v, r, *(buf_weight_p+i));
	}
	Vector6 result;

	// solve ls
	ls.finish();
	ls.solve(result);

	return result;
}



}

