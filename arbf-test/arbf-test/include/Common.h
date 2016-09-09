//
//  Common.h
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#ifndef ARBF_TEST_COMMON_H
#define ARBF_TEST_COMMON_H

#include <string>
#include <cmath>
#include <array>
#include <functional>
#include "Config.h"

#define SQ(x) ((x) * (x))
//#define EPSILON (1e-5)
//#define C0 (0.2) // C0 parameter in paper

struct Vertex {
    Vertex(): x(0.0), y(0.0), z(0.0), intensity(1.0),
              eigVec1(std::array<double, 3>({1, 0, 0})), eigVec2(std::array<double, 3>({0, 1, 0})),
              eigVec3(std::array<double, 3>({0, 0, 1})) {
        // nothing here
    }

    Vertex(double x, double y, double z, double intensity=1.0): x(x), y(y), z(z), intensity(intensity),
                                                                eig1(0.0), eig2(0.0), eig3(0.0),
                                                                eigVec1(std::array<double, 3>({1, 0, 0})),
                                                                eigVec2(std::array<double, 3>({0, 1, 0})),
                                                                eigVec3(std::array<double, 3>({0, 0, 1})) {
        // nothing here
    }

    bool operator==(const Vertex &other) const {
        return (std::abs(x - other.x) < Config::epsilon) &&
                (std::abs(y - other.y) < Config::epsilon) &&
                (std::abs(z - other.z) < Config::epsilon);
    }

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
};

struct Face {
    Face(int a, int b, int c, double intensity=-1.0): a(a), b(b), c(c), intensity(intensity) {
        // nothing here
    }

    bool operator==(const Face &other) const {
        bool cond1 = (a == other.a) || (a == other.b) || (a == other.c);
        bool cond2 = (b == other.a) || (b == other.b) || (b == other.c);
        bool cond3 = (c == other.a) || (c == other.b) || (c == other.c);
        return cond1 && cond2 && cond3;
    }

    void computeCenter(const Vertex &v0, const Vertex &v1, const Vertex &v2) {
        double x = (v0.x + v1.x + v2.x) / 3.0;
        double y = (v0.y + v1.y + v2.y) / 3.0;
        double z = (v0.z + v1.z + v2.z) / 3.0;
        center.x = x;
        center.y = y;
        center.z = z;
        center.intensity = -1.0;
    }

    int a;
    int b;
    int c;
    double intensity;
    Vertex center; // face center
};

struct Edge {
    Edge(int a, int b, double intensity=-1.0): a(a), b(b), intensity(intensity) {
        // nothing here
    }

    bool operator==(const Edge &other) const {
        return (a == other.a && b == other.b) || (a == other.b && b == other.a);
    }

    void computeCenter(const Vertex &v0, const Vertex &v1) {
        double x = (v0.x + v1.x) / 2.0;
        double y = (v0.y + v1.y) / 2.0;
        double z = (v0.z + v1.z) / 2.0;
        center.x = x;
        center.y = y;
        center.z = z;
        center.intensity = -1.0;
    }

    int a;
    int b;
    double intensity;
    Vertex center; // edge center
};

struct Tetrahedron {
    Tetrahedron(int a, int b, int c, int d, double intensity=-1.0):
            a(a), b(b), c(c), d(d), intensity(intensity), f1(b, c, d), f2(a, c, d), f3(a, b, d), f4(a, b, c) {
        // nothing here
    }

    bool operator==(const Tetrahedron &other) const {
        bool cond1 = (a == other.a) || (a == other.b) || (a == other.c) || (a == other.d);
        bool cond2 = (b == other.a) || (b == other.b) || (b == other.c) || (b == other.d);
        bool cond3 = (c == other.a) || (c == other.b) || (c == other.c) || (c == other.d);
        bool cond4 = (d == other.a) || (d == other.b) || (d == other.c) || (d == other.d);
        return cond1 && cond2 && cond3 && cond4;
    }

    void computeCenter(const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3) {
        double x = (v0.x + v1.x + v2.x + v3.x) / 4.0;
        double y = (v0.y + v1.y + v2.y + v3.y) / 4.0;
        double z = (v0.z + v1.z + v2.z + v3.z) / 4.0;
        center.x = x;
        center.y = y;
        center.z = z;
        center.intensity = -1.0;
    }

    int a;
    int b;
    int c;
    int d;
    Face f1;
    Face f2;
    Face f3;
    Face f4;
    double intensity;
    Vertex center;
};

namespace std {
    template<>
    struct hash<Face> {
        size_t operator()(const Face &f) const {
            size_t hash1 = std::hash<int>()(f.a);
            size_t hash2 = std::hash<int>()(f.b) >> 1;
            size_t hash3 = std::hash<int>()(f.c) << 1;
            return hash1 ^ hash2 ^ hash3;
        }
    };

    template<>
    struct std::hash<Edge> {
        size_t operator()(const Edge &e) const {
            size_t hash1 = std::hash<int>()(e.a);
            size_t hash2 = std::hash<int>()(e.b) >> 1;
            return hash1 ^ hash2;
        }
    };
}

#endif //ARBF_TEST_COMMON_H
