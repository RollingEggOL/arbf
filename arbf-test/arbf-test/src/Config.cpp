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
bool Config::isUsingISORBF = false;
bool Config::isPrintingDebugInfo = false;
int Config::dim = 2; // 2D model by default
