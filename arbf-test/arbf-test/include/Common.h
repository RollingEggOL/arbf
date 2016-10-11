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
#include <vector>
#include <algorithm>
#include <functional>
#include <stdexcept>

#define SQ(x) ((x) * (x))

struct DirectionSamplings {
    DirectionSamplings();
    DirectionSamplings(const char *filename);
    DirectionSamplings(const DirectionSamplings &other);
    DirectionSamplings &operator=(const DirectionSamplings &other);
    void importFromFile(const char *filename);
    std::array<float, 3> operator[](unsigned index) const;
    std::vector<std::array<float, 3>> getDirections() const;
    unsigned getNumDirections() const;

private:
    std::vector<std::array<float, 3>> m_data; // spherical direction samplings
};

struct Vertex {
    Vertex();
    Vertex(double x, double y, double z, double intensity=1.0);
    Vertex &operator=(const Vertex &other);
    bool operator==(const Vertex &other) const;

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
    TriangleFace(int a, int b, int c, double intensity=-1.0);
    bool operator==(const TriangleFace &other) const;
    void computeCenter(const Vertex &v0, const Vertex &v1, const Vertex &v2);

    int a;
    int b;
    int c;
    double intensity;
    Vertex center; // face center
};

struct QuadrangleFace {
    QuadrangleFace(int a, int b, int c, int d, double intensity=-1.0);
    bool operator==(const QuadrangleFace &other) const;
    void computeCenter(const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3);

    int a;
    int b;
    int c;
    int d;
    double intensity;
    Vertex center; // quadrangle center
};

struct Edge {
    Edge(int a, int b, double intensity=-1.0);
    bool operator==(const Edge &other) const;
    void computeCenter(const Vertex &v0, const Vertex &v1);

    int a;
    int b;
    double intensity;
    Vertex center; // edge center
};

struct Tetrahedron {
    Tetrahedron(int a, int b, int c, int d, double intensity=-1.0);
    bool operator==(const Tetrahedron &other) const;
    void computeCenter(const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3);

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
    Hexahedron(int a, int b, int c, int d, int e, int f, int g, int h, double intensity=-1.0);
    bool operator==(const Hexahedron &other) const;
    void computeCenter(const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3,
                       const Vertex &v4, const Vertex &v5, const Vertex &v6, const Vertex &v7);

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
