//
//  Common.h
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#ifndef __arbf_test_Common_h__
#define __arbf_test_Common_h__

#include <string>
#include <array>

#define SQ(x) ((x) * (x))
//#define EPSILON (1e-5)
//#define C0 (0.2) // C0 parameter in paper

typedef struct _vertex {
    double x;
    double y;
    double z;
    double intensity;
    double eig1;
    double eig2;
    double eig3;
    std::array<double, 3> eigVec1;
    std::array<double, 3> eigVec2;
    std::array<double, 3> eigVec3;
} Vertex;

typedef struct _face {
    int a;
    int b;
    int c;
    double intensity;
} Face;

typedef struct _edge {
    int a;
    int b;
    double intensity;
} Edge;

struct edgeHasher {
    size_t operator() (const Edge& edge) const {
        std::string temp = std::to_string(edge.a) + std::to_string(edge.b) + std::to_string(edge.a + edge.b) + std::to_string(edge.a * edge.b);
        return temp.length();
    }
};

struct edgeComparator {
    bool operator() (const Edge& lhs, const Edge& rhs) const {
        if ((lhs.a == rhs.a && lhs.b == rhs.b) ||
            (lhs.a == rhs.b && lhs.b == rhs.a)) {
            return true;
        } else {
            return false;
        }
    }
};

enum BasisType { MQ, IMQ, Gaussian, TPS }; // type of basis functions

#endif /* defined(__arbf_test_Common_h__) */
