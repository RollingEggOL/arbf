//
// Created by Ke Liu on 10/4/16.
//

#include <cstdlib>
#include "../include/Config.h"
#include "../include/Common.h"

DirectionSamplings::DirectionSamplings() {}

DirectionSamplings::DirectionSamplings(const char *filename) {
    this->importFromFile(filename);
}

DirectionSamplings::DirectionSamplings(const DirectionSamplings &other) {
    this->m_data.clear();

    for (auto dir: other.m_data) {
        this->m_data.push_back(dir);
    }
}

DirectionSamplings &DirectionSamplings::operator=(const DirectionSamplings &other) {
    if (this == &other) {
        return *this;
    }

    this->m_data.clear();

    for (auto dir: other.m_data) {
        this->m_data.push_back(dir);
    }

    return *this;
}

std::array<float, 3> DirectionSamplings::operator[](unsigned index) const {
    return this->m_data[index];
};

void DirectionSamplings::importFromFile(const char *filename) {
    printf("Reading direction sampling file %s\n", filename);
    int num;
    float x, y, z;
    char line[256];
    FILE *fin = nullptr;

    if ((fin = fopen(filename, "r")) == nullptr) {
        fprintf(stderr, "Error: Reading direction sampling file %s \n", filename);
        exit(EXIT_FAILURE);
    }

    fscanf(fin, "%d\n", &num);

    for (int i = 0; i < num; i++) {
        fscanf(fin, "%f %f %f\n", &x, &y, &z);
        m_data.push_back({ x, y, z });
    }

    fclose(fin);
    printf("%d direction samplings read.\n", num);
}

std::vector<std::array<float, 3>> DirectionSamplings::getDirections() const {
    return m_data;
}

unsigned DirectionSamplings::getNumDirections() const {
    return m_data.size();
}

Vertex::Vertex():
        x(0.0), y(0.0), z(0.0), intensity(1.0),
        eigVec1({ 1, 0, 0 }),
        eigVec2({ 0, 1, 0 }),
        eigVec3({ 0, 0, 1 }) {
    // nothing here
}

Vertex::Vertex(double x, double y, double z, double intensity):
        x(x), y(y), z(z), intensity(intensity),
        eig1(0.0), eig2(0.0), eig3(0.0),
        eigVec1({ 1, 0, 0 }),
        eigVec2({ 0, 1, 0 }),
        eigVec3({ 0, 0, 1 }) {
    // nothing here
}

Vertex &Vertex::operator=(const Vertex &other) {
    if (this == &other) {
        return *this;
    }

    x = other.x;
    y = other.y;
    z = other.z;
    intensity = other.intensity;
    eig1 = other.eig1;
    eig2 = other.eig2;
    eig3 = other.eig3;
    eigVec1 = other.eigVec1;
    eigVec2 = other.eigVec2;
    eigVec3 = other.eigVec3;
    return *this;
}

bool Vertex::operator==(const Vertex &other) const {
    return (std::abs(x - other.x) < Config::epsilon) &&
           (std::abs(y - other.y) < Config::epsilon) &&
           (std::abs(z - other.z) < Config::epsilon);
}

TriangleFace::TriangleFace(int a, int b, int c, double intensity):
        a(a), b(b), c(c), intensity(intensity) {
    // nothing here
}

bool TriangleFace::operator==(const TriangleFace &other) const {
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

void TriangleFace::computeCenter(const Vertex &v0, const Vertex &v1, const Vertex &v2) {
    double x = (v0.x + v1.x + v2.x) / 3.0;
    double y = (v0.y + v1.y + v2.y) / 3.0;
    double z = (v0.z + v1.z + v2.z) / 3.0;
    center.x = x;
    center.y = y;
    center.z = z;
    center.intensity = -1.0;
}

QuadrangleFace::QuadrangleFace(int a, int b, int c, int d, double intensity):
        a(a), b(b), c(c), d(d), intensity(intensity) {
    // nothing here
}

bool QuadrangleFace::operator==(const QuadrangleFace &other) const {
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

void QuadrangleFace::computeCenter(const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3) {
    double x = (v0.x + v1.x + v2.x + v3.x) / 4.0;
    double y = (v0.y + v1.y + v2.y + v3.y) / 4.0;
    double z = (v0.z + v1.z + v2.z + v3.z) / 4.0;
    center.x = x;
    center.y = y;
    center.z = z;
    center.intensity = -1.0;
}

Edge::Edge(int a, int b, double intensity): a(a), b(b), intensity(intensity) {
    // nothing here
}

bool Edge::operator==(const Edge &other) const {
    return (a == other.a && b == other.b) || (a == other.b && b == other.a);
}

void Edge::computeCenter(const Vertex &v0, const Vertex &v1) {
    double x = (v0.x + v1.x) / 2.0;
    double y = (v0.y + v1.y) / 2.0;
    double z = (v0.z + v1.z) / 2.0;
    center.x = x;
    center.y = y;
    center.z = z;
    center.intensity = -1.0;
}

Tetrahedron::Tetrahedron(int a, int b, int c, int d, double intensity):
        a(a), b(b), c(c), d(d), intensity(intensity),
        f1(b, c, d), f2(a, c, d), f3(a, b, d), f4(a, b, c) {
    // nothing here
}

bool Tetrahedron::operator==(const Tetrahedron &other) const {
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

void Tetrahedron::computeCenter(const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3) {
    double x = (v0.x + v1.x + v2.x + v3.x) / 4.0;
    double y = (v0.y + v1.y + v2.y + v3.y) / 4.0;
    double z = (v0.z + v1.z + v2.z + v3.z) / 4.0;
    center.x = x;
    center.y = y;
    center.z = z;
    center.intensity = -1.0;
}

Hexahedron::Hexahedron(int a, int b, int c, int d, int e, int f, int g, int h, double intensity):
        a(a), b(b), c(c), d(d), e(e), f(f), g(g), h(h), intensity(intensity),
        f1(a, d, c, b), f2(e, f, g, h), f3(a, b, f, e), f4(b, c, g, f), f5(c, d, h, g), f6(a, e, h, d) {
    // nothing here
}

bool Hexahedron::operator==(const Hexahedron &other) const {
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

void Hexahedron::computeCenter(const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3,
                               const Vertex &v4, const Vertex &v5, const Vertex &v6, const Vertex &v7) {
    double x = (v0.x + v1.x + v2.x + v3.x + v4.x + v5.x + v6.x + v7.x) / 8.0;
    double y = (v0.y + v1.y + v2.y + v3.y + v4.y + v5.y + v6.y + v7.y) / 8.0;
    double z = (v0.z + v1.z + v2.z + v3.z + v4.z + v5.z + v6.z + v7.z) / 8.0;
    center.x = x;
    center.y = y;
    center.z = z;
    center.intensity = -1.0;
}
