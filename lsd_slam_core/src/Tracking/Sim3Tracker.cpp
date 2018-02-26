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

#include "Sim3Tracker.h"

#include <timings.h>
#include "DataStructures/Frame.h"
#include "Tracking/TrackingReference.h"
#include "util/globalFuncs.h"
#include "Tracking/least_squares.h"

namespace lsd_slam
{



#define callOptimized(function, arguments) function arguments



Sim3Tracker::Sim3Tracker(int w, int h, Eigen::Matrix3f K)
{
	width = w;
	height = h;

	this->K = K;
	fx = K(0,0);
	fy = K(1,1);
	cx = K(0,2);
	cy = K(1,2);

	settings = DenseDepthTrackerSettings();


	KInv = K.inverse();
	fxi = KInv(0,0);
	fyi = KInv(1,1);
	cxi = KInv(0,2);
	cyi = KInv(1,2);


	buf_warped_residual = new float[w*h];
	buf_warped_weights = new float[w*h];
	buf_warped_dx = new float[w*h];
	buf_warped_dy = new float[w*h];
	buf_warped_x = new float[w*h];
	buf_warped_y = new float[w*h];
	buf_warped_z = new float[w*h];

	buf_d = new float[w*h];
	buf_residual_d = new float[w*h];
	buf_idepthVar = new float[w*h];
	buf_warped_idepthVar = new float[w*h];
	buf_weight_p = new float[w*h];
	buf_weight_d = new float[w*h];

	buf_weight_Huber = new float[w*h];
	buf_weight_VarP = new float[w*h];
	buf_weight_VarD = new float[w*h];

	buf_warped_size = 0;


	
	lastResidual = 0;
	iterationNumber = 0;
	lastDepthResidual = lastPhotometricResidual = lastDepthResidualUnweighted = lastPhotometricResidualUnweighted = lastResidualUnweighted = 0;
	pointUsage = 0;

}

Sim3Tracker::~Sim3Tracker()
{


	delete[] buf_warped_residual;
	delete[] buf_warped_weights;
	delete[] buf_warped_dx;
	delete[] buf_warped_dy;
	delete[] buf_warped_x;
	delete[] buf_warped_y;
	delete[] buf_warped_z;

	delete[] buf_d;
	delete[] buf_residual_d;
	delete[] buf_idepthVar;
	delete[] buf_warped_idepthVar;
	delete[] buf_weight_p;
	delete[] buf_weight_d;

	delete[] buf_weight_Huber;
	delete[] buf_weight_VarP;
	delete[] buf_weight_VarD;
}


Sim3 Sim3Tracker::trackFrameSim3(
		TrackingReference* reference,
		Frame* frame,
		const Sim3& frameToReference_initialEstimate,
		int startLevel, int finalLevel)
{



	diverged = false;

	affineEstimation_a = 1; affineEstimation_b = 0;


	// ============ track frame ============
    Sim3 referenceToFrame = frameToReference_initialEstimate.inverse();
	NormalEquationsLeastSquares7 ls7;


	int numCalcResidualCalls[PYRAMID_LEVELS];
	int numCalcWarpUpdateCalls[PYRAMID_LEVELS];

	Sim3ResidualStruct finalResidual;

	bool warp_update_up_to_date = false;

	for(int lvl=startLevel;lvl >= finalLevel;lvl--)
	{
		numCalcResidualCalls[lvl] = 0;
		numCalcWarpUpdateCalls[lvl] = 0;

		if(settings.maxItsPerLvl[lvl] == 0)
			continue;

		reference->makePointCloud(lvl);

		// evaluate baseline-residual.

		callOptimized(calcSim3Buffers, (reference, frame, referenceToFrame, lvl));




		if(buf_warped_size < 0.5 * MIN_GOODPERALL_PIXEL_ABSMIN * (width>>lvl)*(height>>lvl) || buf_warped_size < 10)
		{
			diverged = true;
			return Sim3();
		}


		Sim3ResidualStruct lastErr = callOptimized(calcSim3WeightsAndResidual,(referenceToFrame));






		numCalcResidualCalls[lvl]++;

		if(useAffineLightningEstimation)
		{
			affineEstimation_a = affineEstimation_a_lastIt;
			affineEstimation_b = affineEstimation_b_lastIt;
		}

		float LM_lambda = settings.lambdaInitial[lvl];

		warp_update_up_to_date = false;
		for(int iteration=0; iteration < settings.maxItsPerLvl[lvl]; iteration++)
		{

			// calculate LS System, result is saved in ls.

			callOptimized(calcSim3LGS,(ls7));

			warp_update_up_to_date = true;
			numCalcWarpUpdateCalls[lvl]++;

			iterationNumber = iteration;


			int incTry=0;
			while(true)
			{
				// solve LS system with current lambda
				Vector7 b = - ls7.b / ls7.num_constraints;
				Matrix7x7 A = ls7.A / ls7.num_constraints;
				for(int i=0;i<7;i++) A(i,i) *= 1+LM_lambda;
				Vector7 inc = A.ldlt().solve(b);
				incTry++;

				float absInc = inc.dot(inc);
				if(!(absInc >= 0 && absInc < 1))
				{
					// ERROR tracking diverged.
					lastSim3Hessian.setZero();
					return Sim3();
				}

				// apply increment. pretty sure this way round is correct, but hard to test.
				Sim3 new_referenceToFrame =Sim3::exp(inc.cast<sophusType>()) * referenceToFrame;
				//Sim3 new_referenceToFrame = referenceToFrame * Sim3::exp((inc));


				// re-evaluate residual

				callOptimized(calcSim3Buffers,(reference, frame, new_referenceToFrame, lvl));




				if(buf_warped_size < 0.5 * MIN_GOODPERALL_PIXEL_ABSMIN * (width>>lvl)*(height>>lvl) || buf_warped_size < 10)
				{
					diverged = true;
					return Sim3();
				}


				Sim3ResidualStruct error = callOptimized(calcSim3WeightsAndResidual,(new_referenceToFrame));






				numCalcResidualCalls[lvl]++;






				// accept inc?
				if(error.mean < lastErr.mean)
				{
					// accept inc
					referenceToFrame = new_referenceToFrame;
					warp_update_up_to_date = false;

					if(useAffineLightningEstimation)
					{
						affineEstimation_a = affineEstimation_a_lastIt;
						affineEstimation_b = affineEstimation_b_lastIt;
					}



					// converged?
					if(error.mean / lastErr.mean > settings.convergenceEps[lvl])
					{

						iteration = settings.maxItsPerLvl[lvl];
					}

					finalResidual = lastErr = error;

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






	// Make sure that there is a warp update at the final position to get the correct information matrix
	if (!warp_update_up_to_date)
	{
		reference->makePointCloud(finalLevel);

		callOptimized(calcSim3Buffers,(reference, frame, referenceToFrame, finalLevel));





	    finalResidual = callOptimized(calcSim3WeightsAndResidual,(referenceToFrame));





	    callOptimized(calcSim3LGS,(ls7));



	}

	lastSim3Hessian = ls7.A;


	if(referenceToFrame.scale() <= 0 )
	{
		diverged = true;
		return Sim3();
	}

	lastResidual = finalResidual.mean;
	lastDepthResidual = finalResidual.meanD;
	lastPhotometricResidual = finalResidual.meanP;


	return referenceToFrame.inverse();
}






void Sim3Tracker::calcSim3Buffers(
		const TrackingReference* reference,
		Frame* frame,
		const Sim3& referenceToFrame,
		int level, bool )
{


	// get static values
	int w = frame->width(level);
	int h = frame->height(level);
	Eigen::Matrix3f KLvl = frame->K(level);
	float fx_l = KLvl(0,0);
	float fy_l = KLvl(1,1);
	float cx_l = KLvl(0,2);
	float cy_l = KLvl(1,2);

	Eigen::Matrix3f rotMat = referenceToFrame.rxso3().matrix().cast<float>();
	Eigen::Matrix3f rotMatUnscaled = referenceToFrame.rotationMatrix().cast<float>();
	Eigen::Vector3f transVec = referenceToFrame.translation().cast<float>();

	// Calculate rotation around optical axis for rotating source frame gradients
	Eigen::Vector3f forwardVector(0, 0, -1);
	Eigen::Vector3f rotatedForwardVector = rotMatUnscaled * forwardVector;
	Eigen::Quaternionf shortestBackRotation;
	shortestBackRotation.setFromTwoVectors(rotatedForwardVector, forwardVector);
	Eigen::Matrix3f rollMat = shortestBackRotation.toRotationMatrix() * rotMatUnscaled;
	float xRoll0 = rollMat(0, 0);
	float xRoll1 = rollMat(0, 1);
	float yRoll0 = rollMat(1, 0);
	float yRoll1 = rollMat(1, 1);


	const Eigen::Vector3f* refPoint_max = reference->posData[level] + reference->numData[level];
	const Eigen::Vector3f* refPoint = reference->posData[level];
	const Eigen::Vector2f* refColVar = reference->colorAndVarData[level];
	const Eigen::Vector2f* refGrad = reference->gradData[level];

	const float* 			frame_idepth = frame->idepth(level);
	const float* 			frame_idepthVar = frame->idepthVar(level);
	const Eigen::Vector4f* 	frame_intensityAndGradients = frame->gradients(level);


	float sxx=0,syy=0,sx=0,sy=0,sw=0;

	float usageCount = 0;

	int idx=0;
	for(;refPoint<refPoint_max; refPoint++, refGrad++, refColVar++)
	{
		Eigen::Vector3f Wxp = rotMat * (*refPoint) + transVec;
		float u_new = (Wxp[0]/Wxp[2])*fx_l + cx_l;
		float v_new = (Wxp[1]/Wxp[2])*fy_l + cy_l;

		// step 1a: coordinates have to be in image:
		// (inverse test to exclude NANs)
		if(!(u_new > 1 && v_new > 1 && u_new < w-2 && v_new < h-2))
			continue;

		*(buf_warped_x+idx) = Wxp(0);
		*(buf_warped_y+idx) = Wxp(1);
		*(buf_warped_z+idx) = Wxp(2);

		Eigen::Vector3f resInterp = getInterpolatedElement43(frame_intensityAndGradients, u_new, v_new, w);


		// save values
#if USE_ESM_TRACKING == 1
		// get rotated gradient of point
		float rotatedGradX = xRoll0 * (*refGrad)[0] + xRoll1 * (*refGrad)[1];
		float rotatedGradY = yRoll0 * (*refGrad)[0] + yRoll1 * (*refGrad)[1];

		*(buf_warped_dx+idx) = fx_l * 0.5f * (resInterp[0] + rotatedGradX);
		*(buf_warped_dy+idx) = fy_l * 0.5f * (resInterp[1] + rotatedGradY);
#else
		*(buf_warped_dx+idx) = fx_l * resInterp[0];
		*(buf_warped_dy+idx) = fy_l * resInterp[1];
#endif


		float c1 = affineEstimation_a * (*refColVar)[0] + affineEstimation_b;
		float c2 = resInterp[2];
		float residual_p = c1 - c2;

		float weight = fabsf(residual_p) < 2.0f ? 1 : 2.0f / fabsf(residual_p);
		sxx += c1*c1*weight;
		syy += c2*c2*weight;
		sx += c1*weight;
		sy += c2*weight;
		sw += weight;


		*(buf_warped_residual+idx) = residual_p;
		*(buf_idepthVar+idx) = (*refColVar)[1];


		// new (only for Sim3):
		int idx_rounded = (int)(u_new+0.5f) + w*(int)(v_new+0.5f);
		float var_frameDepth = frame_idepthVar[idx_rounded];
		float ref_idepth = 1.0f / Wxp[2];
		*(buf_d+idx) = 1.0f / (*refPoint)[2];
		if(var_frameDepth > 0)
		{
			float residual_d = ref_idepth - frame_idepth[idx_rounded];
			*(buf_residual_d+idx) = residual_d;
			*(buf_warped_idepthVar+idx) = var_frameDepth;
		}
		else
		{
			*(buf_residual_d+idx) = -1;
			*(buf_warped_idepthVar+idx) = -1;
		}



		idx++;

		float depthChange = (*refPoint)[2] / Wxp[2];
		usageCount += depthChange < 1 ? depthChange : 1;
	}
	buf_warped_size = idx;


	pointUsage = usageCount / (float)reference->numData[level];

	affineEstimation_a_lastIt = sqrtf((syy - sy*sy/sw) / (sxx - sx*sx/sw));
	affineEstimation_b_lastIt = (sy - affineEstimation_a_lastIt*sx)/sw;





}



Sim3ResidualStruct Sim3Tracker::calcSim3WeightsAndResidual(
		const Sim3& referenceToFrame)
{
	float tx = referenceToFrame.translation()[0];
	float ty = referenceToFrame.translation()[1];
	float tz = referenceToFrame.translation()[2];

	Sim3ResidualStruct sumRes;
	memset(&sumRes, 0, sizeof(Sim3ResidualStruct));



	for(int i=0;i<buf_warped_size;i++)
	{
		float px = *(buf_warped_x+i);	// x'
		float py = *(buf_warped_y+i);	// y'
		float pz = *(buf_warped_z+i);	// z'

		float d = *(buf_d+i);	// d

		float rp = *(buf_warped_residual+i); // r_p
		float rd = *(buf_residual_d+i);	 // r_d

		float gx = *(buf_warped_dx+i);	// \delta_x I
		float gy = *(buf_warped_dy+i);  // \delta_y I

		float s = settings.var_weight * *(buf_idepthVar+i);	// \sigma_d^2
		float sv = settings.var_weight * *(buf_warped_idepthVar+i);	// \sigma_d^2'


		// calc dw/dd (first 2 components):
		float g0 = (tx * pz - tz * px) / (pz*pz*d);
		float g1 = (ty * pz - tz * py) / (pz*pz*d);
		float g2 = (pz - tz) / (pz*pz*d);

		// calc w_p
		float drpdd = gx * g0 + gy * g1;	// ommitting the minus
		float w_p = 1.0f / (cameraPixelNoise2 + s * drpdd * drpdd);

		float w_d = 1.0f / (sv + g2*g2*s);

		float weighted_rd = fabs(rd*sqrtf(w_d));
		float weighted_rp = fabs(rp*sqrtf(w_p));


		float weighted_abs_res = sv > 0 ? weighted_rd+weighted_rp : weighted_rp;
		float wh = fabs(weighted_abs_res < settings.huber_d ? 1 : settings.huber_d / weighted_abs_res);

		if(sv > 0)
		{
			sumRes.sumResD += wh * w_d * rd*rd;
			sumRes.numTermsD++;
		}

		sumRes.sumResP += wh * w_p * rp*rp;
		sumRes.numTermsP++;



		*(buf_weight_p+i) = wh * w_p;

		if(sv > 0) {
			*(buf_weight_d+i) = wh * w_d;
		}
		else {
			*(buf_weight_d+i) = 0;
		}

	}

	sumRes.mean = (sumRes.sumResD + sumRes.sumResP) / (sumRes.numTermsD + sumRes.numTermsP);
	sumRes.meanD = (sumRes.sumResD) / (sumRes.numTermsD);
	sumRes.meanP = (sumRes.sumResP) / (sumRes.numTermsP);


	return sumRes;
}



void Sim3Tracker::calcSim3LGS(NormalEquationsLeastSquares7 &ls7)
{
	NormalEquationsLeastSquares4 ls4;
	NormalEquationsLeastSquares ls6;
	ls6.initialize(width*height);
	ls4.initialize(width*height);

	for(int i=0;i<buf_warped_size;i++)
	{
		float px = *(buf_warped_x+i);	// x'
		float py = *(buf_warped_y+i);	// y'
		float pz = *(buf_warped_z+i);	// z'

		float wp = *(buf_weight_p+i);	// wr/wp
		float wd = *(buf_weight_d+i);	// wr/wd

		float rp = *(buf_warped_residual+i); // r_p
		float rd = *(buf_residual_d+i);	 // r_d

		float gx = *(buf_warped_dx+i);	// \delta_x I
		float gy = *(buf_warped_dy+i);  // \delta_y I


		float z = 1.0f / pz;
		float z_sqr = 1.0f / (pz*pz);
		Vector6 v;
		Vector4 v4;
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

		// new:
		v4[0] = z_sqr;
		v4[1] = z_sqr * py;
		v4[2] = -z_sqr * px;
		v4[3] = z;

		ls6.update(v, rp, wp);		// Jac = - v
		ls4.update(v4, rd, wd);	// Jac = v4

	}

	ls4.finishNoDivide();
	ls6.finishNoDivide();


	ls7.initializeFrom(ls6, ls4);

}


}
