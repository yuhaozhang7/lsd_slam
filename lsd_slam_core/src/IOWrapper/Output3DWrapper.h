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

#pragma once
#include <string>
#include <vector>
#include <Eigen/Eigen>



namespace lsd_slam
{

class Frame;

class KeyFrameGraph;
class Frame;



/**
 * Virtual 3D display object.
 */
class Output3DWrapper
{
public:
	virtual ~Output3DWrapper() {};


	virtual void publishKeyframeGraph(KeyFrameGraph* ) {};

	// publishes a keyframe. if that frame already exists, it is overwritten, otherwise it is added.
	virtual void publishKeyframe(Frame* ) {};


	// published a tracked frame that did not become a keyframe (yet; i.e. has no depth data)
	virtual void publishTrackedFrame(Frame* ) {};

	// publishes graph and all constraints, as well as updated KF poses.
	virtual void publishTrajectory(std::vector<Eigen::Matrix<float, 3, 1>> , std::string ) {};
	virtual void publishTrajectoryIncrement(Eigen::Matrix<float, 3, 1> , std::string ) {};

    virtual void publishDebugInfo(Eigen::Matrix<float, 20, 1> ) {};

};
}
