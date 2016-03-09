//
//  Config.h
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#ifndef __arbf_test_Config_h__
#define __arbf_test_Config_h__

#include <string>

//////////////////////////////////////////////////////
// choose one of the two
//#define NEIGHBOR_1RING // use 1-ring neighborhood
//#define NEIGHBOR_2RING // use 2-ring neighborhood

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
    static bool isUsingISORBF; // if use ISO-RBF for eigenvalues
    static bool isPrintingDebugInfo; // if printing out debug info
    static int dim; // model dimension
};

#endif /* defined(__arbf_test_Config_h__) */