//
//  Mesh.cpp
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#include <cstdio>
#include <cassert>
#include <cmath>
#include <limits>
#include <string>
#include <cstdlib>
#include <ctime>
#include <array>
#include <vector>
#include <unordered_set>
#include "../include/Mesh.h"

using namespace std;

TriMesh::TriMesh(): m_nv(0), m_nf(0), m_ne(0),
m_minX(numeric_limits<double>::max()), m_maxX(numeric_limits<double>::lowest()),
m_minY(numeric_limits<double>::max()), m_maxY(numeric_limits<double>::lowest()),
m_minZ(numeric_limits<double>::max()), m_maxZ(numeric_limits<double>::lowest()) {}

TriMesh::TriMesh(const TriMesh &other) {
    m_vertex = other.m_vertex;
    m_face = other.m_face;
    m_edge = other.m_edge;
    m_nv = other.m_nv;
    m_nf = other.m_nf;
    m_ne = other.m_ne;
    m_minX = other.m_minX;
    m_maxX = other.m_maxX;
    m_minY = other.m_minY;
    m_maxY = other.m_maxY;
    m_minZ = other.m_minZ;
    m_maxZ = other.m_maxZ;
}

TriMesh::~TriMesh() {}

TriMesh &TriMesh::operator=(const TriMesh &other) {
    if (this == &other) {
        return *this;
    }
    
    this->~TriMesh();
    m_vertex = other.m_vertex;
    m_face = other.m_face;
    m_edge = other.m_edge;
    m_nv = other.m_nv;
    m_nf = other.m_nf;
    m_ne = other.m_ne;
    m_minX = other.m_minX;
    m_maxX = other.m_maxX;
    m_minY = other.m_minY;
    m_maxY = other.m_maxY;
    m_minZ = other.m_minZ;
    m_maxZ = other.m_maxZ;
    return *this;
}

unsigned TriMesh::getNumVertices() const {
    return m_nv;
}

unsigned TriMesh::getNumFaces() const {
    return m_nf;
}

unsigned TriMesh::getNumEdges() const {
    return m_ne;
}

const vector<Vertex> &TriMesh::getVertices() const {
    return m_vertex;
}

const vector<Face> &TriMesh::getFaces() const {
    return m_face;
}

const std::unordered_set<Edge, edgeHasher, edgeComparator> &TriMesh::getEdges() const {
    return m_edge;
}

double TriMesh::getMinX() const {
    return m_minX;
}

double TriMesh::getMaxX() const {
    return m_maxX;
}

double TriMesh::getMinY() const {
    return m_minY;
}

double TriMesh::getMaxY() const {
    return m_maxY;
}

double TriMesh::getMinZ() const {
    return m_minZ;
}

double TriMesh::getMaxZ() const {
    return m_maxZ;
}

array<double, 2> TriMesh::getMinMaxX() const {
    return { m_minX, m_maxX };
}

array<double, 2> TriMesh::getMinMaxY() const {
    return { m_minY, m_maxY };
}

array<double, 2> TriMesh::getMinMaxZ() const {
    return { m_minZ, m_maxZ };
}

void TriMesh::setNumVertices(unsigned value) {
    m_nv = value;
}

void TriMesh::setNumFaces(unsigned value) {
    m_nf = value;
}

void TriMesh::setNumEdges(unsigned value) {
    m_ne = value;
}

void TriMesh::setMinX(double value) {
    m_minX = value;
}

void TriMesh::setMaxX(double value) {
    m_maxX = value;
}

void TriMesh::setMinY(double value) {
    m_minY = value;
}

void TriMesh::setMaxY(double value) {
    m_maxY = value;
}

void TriMesh::setMinZ(double value) {
    m_minZ = value;
}

void TriMesh::setMaxZ(double value) {
    m_maxZ = value;
}

void TriMesh::setMinMaxX(const std::array<double, 2> &value) {
    m_minX = value[0];
    m_maxX = value[1];
}

void TriMesh::setMinMaxY(const std::array<double, 2> &value) {
    m_minY = value[0];
    m_maxY = value[1];
}

void TriMesh::setMinMaxZ(const std::array<double, 2> &value) {
    m_minZ = value[0];
    m_maxZ = value[1];
}

void TriMesh::addVertex(const Vertex &vertex) {
    m_vertex.push_back(vertex);
}

void TriMesh::addFace(const Face &face) {
    m_face.push_back(face);
}

void TriMesh::addEdge(const Edge &edge) {
    m_edge.insert(edge);
}

void TriMesh::write(const char *path) {
    printf("Writing mesh file %s\n", path);
    FILE *fout = nullptr;
    
    if ((fout = fopen(path, "w")) == nullptr) {
        fprintf(stderr, "Error: Writing mesh file %s \n", path);
        exit(EXIT_FAILURE);
    }
    
    fputs("OFF\n", fout);
    fprintf(fout, "%d %d %d\n", m_nv, m_nf, m_ne);
    
    for (int v = 0; v < m_nv; v++) {
        fprintf(fout, "%lf %lf %lf\n", m_vertex[v].x, m_vertex[v].y, m_vertex[v].z);
    }
    
    for (int f = 0; f < m_nf; f++) {
        fprintf(fout, "3 %d %d %d\n", m_face[f].a, m_face[f].b, m_face[f].c);
    }
    
    fclose(fout);
}

void TriMesh::deform() {
    srand(static_cast<unsigned int>(time(NULL)));
//    double t1 = 0.5, t2 = 0.5, t3 = 0.5;
    double t1 = (double) rand() / RAND_MAX;
    double t2 = (double) rand() / RAND_MAX;
    double t3 = (double) rand() / RAND_MAX;
    
    // fixed mesh nodes: 0, 3, 5
    // moving mesh nodes: 1, 4, 2
    // x = x0 + t * (x1 - x0)
    // y = y0 + t * (y1 - y0)
    m_vertex[1].x = m_vertex[0].x + t1 * (m_vertex[3].x - m_vertex[0].x);
    m_vertex[1].y = m_vertex[0].y + t1 * (m_vertex[3].y - m_vertex[0].y);
    m_vertex[4].x = m_vertex[3].x + t2 * (m_vertex[5].x - m_vertex[3].x);
    m_vertex[4].y = m_vertex[3].y + t2 * (m_vertex[5].y - m_vertex[3].y);
    m_vertex[2].x = m_vertex[5].x + t3 * (m_vertex[0].x - m_vertex[5].x);
    m_vertex[2].y = m_vertex[5].y + t3 * (m_vertex[0].y - m_vertex[5].y);
}
