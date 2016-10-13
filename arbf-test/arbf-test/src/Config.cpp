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
Config::BasisFunctionType Config::basisType = Config::BasisFunctionType::IMQ;
Config::InterpolationScheme Config::interpolationScheme = Config::InterpolationScheme::global;
Config::MeshType Config::meshType = Config::MeshType::Tetrahedron;
float Config::disturbanceFactor = 0.0f;
bool Config::isDisturbanceEnabled = false;
bool Config::isDebugEnabled = false;
int Config::problemDim = 2;
unsigned Config::numEvalPoints = 0;
double Config::epsilon = 0.0;
double Config::c0 = 0.0;

std::map<std::string, Config::BasisFunctionType> Config::basisFunctionMap = {
        std::make_pair("MQ", BasisFunctionType::MQ),
        std::make_pair("IMQ", BasisFunctionType::IMQ),
        std::make_pair("Gaussian", BasisFunctionType::Gaussian),
        std::make_pair("TPS", BasisFunctionType::TPS)
};

std::map<std::string, Config::InterpolationScheme> Config::schemeMap = {
        std::make_pair("local", InterpolationScheme::local),
        std::make_pair("global", InterpolationScheme::global)
};

std::map<std::string, Config::MeshType> Config::meshMap = {
        std::make_pair("Tetrahedron", MeshType::Tetrahedron),
        std::make_pair("Hexahedron", MeshType::Hexahedron)
};

void Config::setBasisFunctionType(const std::string &basis) {
    if (basis != "MQ" && basis != "IMQ" && basis != "Gaussian" && basis != "TPS") {
        exit(EXIT_FAILURE);
    }

    basisType = basisFunctionMap[basis];
}

void Config::setInterpolationScheme(const std::string &scheme) {
    if (scheme != "local" && scheme != "global") {
        exit(EXIT_FAILURE);
    }

    interpolationScheme = schemeMap[scheme];
}

void Config::setMeshType(const std::string &type) {
    if (type != "Tetrahedron" && type != "Hexahedron") {
        exit(EXIT_FAILURE);
    }

    meshType = meshMap[type];
}
