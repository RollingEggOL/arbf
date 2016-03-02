//
//  Common.cpp
//  arbf-image-super_resolution
//
//  Created by Ke Liu on 2/26/15.
//  Copyright (c) 2015 Ke Liu. All rights reserved.
//

#include "../include/Config.h"

using namespace std;

const double Config::EPSILON = 1e-5;
const double Config::C0 = 0.2;
bool Config::isUsing2Ring = false;
bool Config::isUsingWeightedInterpolation = false;
bool Config::isUsingISORBF = false;

//string Config::loggerName = string();
//const string Config::loggerINFO = "_info_";
//const string Config::loggerWARNING = "_warning_";
//const string Config::loggerERROR = "_error_";
//const string Config::loggerFATAL = "_fatal_";