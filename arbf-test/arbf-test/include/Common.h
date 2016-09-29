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
#include <algorithm>
#include <functional>
#include <stdexcept>
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

struct TriangleFace {
    TriangleFace(int a, int b, int c, double intensity=-1.0): a(a), b(b), c(c), intensity(intensity) {
        // nothing here
    }

    bool operator==(const TriangleFace &other) const {
        std::array<int, 3> arr1 = { a, b, c };
        std::array<int, 3> arr2 = { other.a, other.b, other.c };
        std::sort(arr1.begin(), arr1.end());
        std::sort(arr2.begin(), arr2.end());

        for (int i = 0; i < 3; i++) {
            if (arr1[i] != arr2[i]) {
                return false;
            }
        }

        return true;
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

struct QuadrangleFace {
    QuadrangleFace(int a, int b, int c, int d, double intensity=-1.0): a(a), b(b), c(c), d(d), intensity(intensity) {
        // nothing here
    }

    bool operator==(const QuadrangleFace &other) const {
        std::array<int, 4> arr1 = { a, b, c, d };
        std::array<int, 4> arr2 = { other.a, other.b, other.c, other.d };
        std::sort(arr1.begin(), arr1.end());
        std::sort(arr2.begin(), arr2.end());

        for (int i = 0; i < 4; i++) {
            if (arr1[i] != arr2[i]) {
                return false;
            }
        }

        return true;
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
    double intensity;
    Vertex center; // quadrangle center
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
        std::array<int, 4> arr1 = { a, b, c, d };
        std::array<int, 4> arr2 = { other.a, other.b, other.c, other.d };
        std::sort(arr1.begin(), arr1.end());
        std::sort(arr2.begin(), arr2.end());

        for (int i = 0; i < 4; i++) {
            if (arr1[i] != arr2[i]) {
                return false;
            }
        }

        return true;
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
    TriangleFace f1;
    TriangleFace f2;
    TriangleFace f3;
    TriangleFace f4;
    double intensity;
    Vertex center;
};

struct Hexahedron {
    Hexahedron(int a, int b, int c, int d, int e, int f, int g, int h, double intensity=-1.0):
        a(a), b(b), c(c), d(d), e(e), f(f), g(g), h(h), intensity(intensity),
        f1(a, d, c, b), f2(e, f, g, h), f3(a, b, f, e), f4(b, c, g, f), f5(c, d, h, g), f6(a, e, h, d) {
        // nothing here
    }

    bool operator==(const Hexahedron &other) const {
        std::array<int, 8> arr1 = { a, b, c, d, e, f, g, h };
        std::array<int, 8> arr2 = { other.a, other.b, other.c, other.d, other.e, other.f, other.g, other.h };
        std::sort(arr1.begin(), arr1.end());
        std::sort(arr2.begin(), arr2.end());

        for (int i = 0; i < 8; i++) {
            if (arr1[i] != arr2[i]) {
                return false;
            }
        }

        return true;
    }

    void computeCenter(const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3,
                       const Vertex &v4, const Vertex &v5, const Vertex &v6, const Vertex &v7) {
        double x = (v0.x + v1.x + v2.x + v3.x + v4.x + v5.x + v6.x + v7.x) / 8.0;
        double y = (v0.y + v1.y + v2.y + v3.y + v4.y + v5.y + v6.y + v7.y) / 8.0;
        double z = (v0.z + v1.z + v2.z + v3.z + v4.z + v5.z + v6.z + v7.z) / 8.0;
        center.x = x;
        center.y = y;
        center.z = z;
        center.intensity = -1.0;
    }

    int a;
    int b;
    int c;
    int d;
    int e;
    int f;
    int g;
    int h;
    QuadrangleFace f1;
    QuadrangleFace f2;
    QuadrangleFace f3;
    QuadrangleFace f4;
    QuadrangleFace f5;
    QuadrangleFace f6;
    double intensity;
    Vertex center;
};

namespace std {
    template<>
    struct hash<TriangleFace> {
        size_t operator()(const TriangleFace &f) const {
            size_t hash1 = std::hash<int>()(f.a);
            size_t hash2 = std::hash<int>()(f.b);
            size_t hash3 = std::hash<int>()(f.c);
            return ((hash1 ^ hash2) >> 1) ^ hash3;
        }
    };

    template<>
    struct hash<QuadrangleFace> {
        size_t operator()(const QuadrangleFace &f) const {
            size_t hash1 = std::hash<int>()(f.a);
            size_t hash2 = std::hash<int>()(f.b);
            size_t hash3 = std::hash<int>()(f.c);
            size_t hash4 = std::hash<int>()(f.d);
            return hash1 ^ hash2 ^ hash3 ^ hash4;
        }
    };

    template<>
    struct hash<Edge> {
        size_t operator()(const Edge &e) const {
            size_t hash1 = std::hash<int>()(e.a);
            size_t hash2 = std::hash<int>()(e.b);
            return hash1 ^ hash2;
        }
    };
}

class NotImplementedException: public std::logic_error {
public:
    NotImplementedException(): std::logic_error("Error: not implemented.") {
        // nothing here
    }
};

#endif //ARBF_TEST_COMMON_H
