//
//  Mesh.cpp
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <array>
#include "../include/Common.h"
#include "../include/Mesh.h"

Mesh::Mesh(): m_nv(0), m_nft(0), m_nfq(0), m_ne(0),
              m_minX(std::numeric_limits<double>::max()), m_maxX(std::numeric_limits<double>::lowest()),
              m_minY(std::numeric_limits<double>::max()), m_maxY(std::numeric_limits<double>::lowest()),
              m_minZ(std::numeric_limits<double>::max()), m_maxZ(std::numeric_limits<double>::lowest()) {
    // nothing here
}

Mesh::Mesh(const Mesh &other) {
    m_vertices = other.m_vertices;
    m_triangleFaces = other.m_triangleFaces;
    m_edges = other.m_edges;
    m_nv = other.m_nv;
    m_nft = other.m_nft;
    m_nfq = other.m_nfq;
    m_ne = other.m_ne;
    m_minX = other.m_minX;
    m_maxX = other.m_maxX;
    m_minY = other.m_minY;
    m_maxY = other.m_maxY;
    m_minZ = other.m_minZ;
    m_maxZ = other.m_maxZ;
}

Mesh::~Mesh() {
    // nothing here
}

Mesh &Mesh::operator=(const Mesh &other) {
    if (this == &other) {
        return *this;
    }

    m_vertices = other.m_vertices;
    m_triangleFaces = other.m_triangleFaces;
    m_edges = other.m_edges;
    m_nv = other.m_nv;
    m_nft = other.m_nft;
    m_nfq = other.m_nfq;
    m_ne = other.m_ne;
    m_minX = other.m_minX;
    m_maxX = other.m_maxX;
    m_minY = other.m_minY;
    m_maxY = other.m_maxY;
    m_minZ = other.m_minZ;
    m_maxZ = other.m_maxZ;
    return *this;
}

unsigned Mesh::getNumVertices() const {
    return m_nv;
}

unsigned Mesh::getNumTriangleFaces() const {
    return m_nft;
}

unsigned Mesh::getNumQuadrangleFaces() const {
    return m_nfq;
}

unsigned Mesh::getNumEdges() const {
    return m_ne;
}

const std::vector<Vertex> &Mesh::getVertices() const {
    return m_vertices;
}

const std::unordered_set<TriangleFace> &Mesh::getTriangleFaces() const {
    return m_triangleFaces;
}

const std::vector<TriangleFace> Mesh::getTriangleFacesAsList() const {
    return std::vector<TriangleFace>(m_triangleFaces.begin(), m_triangleFaces.end());
}

const std::unordered_set<QuadrangleFace> &Mesh::getQuadrangleFaces() const {
    return m_quadrangleFaces;
}

const std::vector<QuadrangleFace> Mesh::getQuadrangleFacesAsList() const {
    return std::vector<QuadrangleFace>(m_quadrangleFaces.begin(), m_quadrangleFaces.end());
}

const std::unordered_set<Edge> &Mesh::getEdges() const {
    return m_edges;
}

const std::vector<Edge> Mesh::getEdgesAsList() const {
    return std::vector<Edge>(m_edges.begin(), m_edges.end());
}

double Mesh::getMinX() const {
    return m_minX;
}

double Mesh::getMaxX() const {
    return m_maxX;
}

double Mesh::getMinY() const {
    return m_minY;
}

double Mesh::getMaxY() const {
    return m_maxY;
}

double Mesh::getMinZ() const {
    return m_minZ;
}

double Mesh::getMaxZ() const {
    return m_maxZ;
}

std::array<double, 2> Mesh::getMinMaxX() const {
    return { m_minX, m_maxX };
}

std::array<double, 2> Mesh::getMinMaxY() const {
    return { m_minY, m_maxY };
}

std::array<double, 2> Mesh::getMinMaxZ() const {
    return { m_minZ, m_maxZ };
}

void Mesh::setNumVertices(unsigned value) {
    m_nv = value;
}

void Mesh::setNumTriangleFaces(unsigned value) {
    m_nft = value;
}

void Mesh::setNumQuadrangleFaces(unsigned value) {
    m_nfq = value;
}

void Mesh::setNumEdges(unsigned value) {
    m_ne = value;
}

void Mesh::setMinX(double value) {
    m_minX = value;
}

void Mesh::setMaxX(double value) {
    m_maxX = value;
}

void Mesh::setMinY(double value) {
    m_minY = value;
}

void Mesh::setMaxY(double value) {
    m_maxY = value;
}

void Mesh::setMinZ(double value) {
    m_minZ = value;
}

void Mesh::setMaxZ(double value) {
    m_maxZ = value;
}

void Mesh::setMinMaxX(const std::array<double, 2> &value) {
    m_minX = value[0];
    m_maxX = value[1];
}

void Mesh::setMinMaxY(const std::array<double, 2> &value) {
    m_minY = value[0];
    m_maxY = value[1];
}

void Mesh::setMinMaxZ(const std::array<double, 2> &value) {
    m_minZ = value[0];
    m_maxZ = value[1];
}

void Mesh::addVertex(const Vertex &vertex) {
    m_vertices.push_back(vertex);
}

void Mesh::addTriangleFace(const TriangleFace &face) {
    m_triangleFaces.insert(face);
}

void Mesh::addQuadrangleFace(const QuadrangleFace &face) {
    m_quadrangleFaces.insert(face);
}

void Mesh::addEdge(const Edge &edge) {
    m_edges.insert(edge);
}

void Mesh::exportToFile(const char *filename) {
    throw new NotImplementedException();
}
