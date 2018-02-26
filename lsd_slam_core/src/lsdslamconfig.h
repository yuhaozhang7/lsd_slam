#ifndef __LSDSLAMConfig_H__
#define __LSDSLAMConfig_H__

#include <map>
#include <string>

#include <iostream>
#include <vector>

#include <SLAMBenchConfiguration.h>



class LSDSLAMConfig : public SLAMBenchConfiguration {



public :

    double kfusage;
    double kfdist;

    bool useFabMap;
    float minUseGrad;
    bool kfreactive;
    bool subpixelstereo;
    double maxLoopClosureCandiates;
    bool poseoptim;
    std::string camera_path;
    std::string package_path;
    int randSeed;

    LSDSLAMConfig()  : SLAMBenchConfiguration()   {

   }




};





#endif
