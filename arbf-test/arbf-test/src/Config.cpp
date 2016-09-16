//
//  Common.cpp
//  arbf-image-super_resolution
//
//  Created by Ke Liu on 2/26/15.
//  Copyright (c) 2015 Ke Liu. All rights reserved.
//

#include <cstdlib>
#include <map>
#include <string>
#include <utility>
#include "../include/Config.h"

int Config::neighborhoodSize = 0;
//bool Config::isUsingISORBF = false;
Config::InterpolationSchemes Config::scheme = Config::InterpolationSchemes::global;
bool Config::isDebugEnabled = false;
int Config::problemDim = 2;
unsigned Config::numEvalPoints = 0;
double Config::epsilon = 0.0;
double Config::c0 = 0.0;

std::map<std::string, Config::InterpolationSchemes> Config::schemeMap = {
        std::make_pair("local", InterpolationSchemes::local),
        std::make_pair("global", InterpolationSchemes::global)
};

void Config::setInterpolationScheme(const std::string &scheme) {
    if (scheme != "local" and scheme != "global") {
        exit(EXIT_FAILURE);
    }

    Config::scheme = Config::schemeMap[scheme];
}
