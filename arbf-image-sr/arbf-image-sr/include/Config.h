//
//  Config.h
//  arbf-image-super_resolution
//
//  Created by Ke Liu on 2/26/15.
//  Copyright (c) 2015 Ke Liu. All rights reserved.
//

#ifndef __arbf_image_super_resolution__Config__
#define __arbf_image_super_resolution__Config__

#include <string>

//////////////////////////////////////////////////////
// choose one of the two
//#define NEIGHBOR_1RING // use 1-ring neighborhood
#define NEIGHBOR_2RING // use 2-ring neighborhood

//////////////////////////////////////////////////////
//#define ISO_RBF // use isotropic RBF interpolation

//////////////////////////////////////////////////////
//#define WEIGHTED // take weighted average

class Config
{
public:
    static const double EPSILON;
    static const double C0; // C0 parameter in paper
    static bool isUsing2Ring; // if use 2-ring neighborhood
    static bool isUsingWeightedInterpolation; // if use weighted interpolation
    static bool isUsingISORBF; // if use ISO-RBF for eigenvalues
    
//    static std::string loggerName;
//    static const std::string loggerINFO;
//    static const std::string loggerWARNING;
//    static const std::string loggerERROR;
//    static const std::string loggerFATAL;
};

#endif /* defined(__arbf_image_super_resolution__Config__) */
