//
//  Common.cpp
//  arbf-image-super_resolution
//
//  Created by Ke Liu on 2/26/15.
//  Copyright (c) 2015 Ke Liu. All rights reserved.
//

#include "../include/Config.h"

using namespace std;

int Config::neighborhoodSize = 0;
//bool Config::isUsingISORBF = false;
bool Config::isDebugEnabled = false;
int Config::problemDim = 2;
unsigned Config::numEvalPoints = 0;
double Config::epsilon = 0.0;
double Config::c0 = 0.0;
