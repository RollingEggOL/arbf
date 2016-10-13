//
//  Config.h
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#ifndef __arbf_test_Config_h__
#define __arbf_test_Config_h__

#include <map>
#include <string>

class Config {
public:
    enum BasisFunctionType { MQ, IMQ, Gaussian, TPS };
    enum InterpolationScheme { global, local };
    enum MeshType { Tetrahedron, Hexahedron };

public:
    static void setBasisFunctionType(const std::string &basis);
    static void setInterpolationScheme(const std::string &scheme);
    static void setMeshType(const std::string &type);

    static int neighborhoodSize;
//    static bool isUsingISORBF; // if use ISO-RBF for eigenvalues
    static BasisFunctionType basisType;
    static InterpolationScheme interpolationScheme;
    static MeshType meshType;
    static bool isDisturbanceEnabled; // if disturbance of vertices, face centers, cell centers enabled
    static float disturbanceFactor; // disturbance factor, defaults to 0f
    static bool isDebugEnabled; // if printing out debug info
    static int problemDim; // model dimension
    static unsigned numEvalPoints; // # of evaluation points
    static double epsilon;
    static double c0; // C0 parameter used in paper
private:
    static std::map<std::string, BasisFunctionType> basisFunctionMap;
    static std::map<std::string, InterpolationScheme> schemeMap;
    static std::map<std::string, MeshType> meshMap;
};

#endif /* defined(__arbf_test_Config_h__) */
