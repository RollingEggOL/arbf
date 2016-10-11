//
//  TriMesh.h
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#ifndef ARBF_TEST_MESH_H
#define ARBF_TEST_MESH_H

#include <vector>
#include <array>
#include <unordered_set>
#include "Common.h"

class Mesh {
public:
    Mesh();
    Mesh(const Mesh &other);
    virtual ~Mesh();
    virtual void exportToFile(const char *filename);
    virtual void update();

    Mesh &operator=(const Mesh &other);
    unsigned getNumVertices() const;
    unsigned getNumEdges() const;
    unsigned getNumTriangleFaces() const;
    unsigned getNumQuadrangleFaces() const;
    const std::vector<Vertex> &getVertices() const;
    const std::unordered_set<Edge> &getEdges() const;
    const std::vector<Edge> getEdgesAsList() const;
    const std::unordered_set<TriangleFace> &getTriangleFaces() const;
    const std::vector<TriangleFace> getTriangleFacesAsList() const;
    const std::unordered_set<QuadrangleFace> &getQuadrangleFaces() const;
    const std::vector<QuadrangleFace> getQuadrangleFacesAsList() const;
    double getMinX() const;
    double getMaxX() const;
    double getMinY() const;
    double getMaxY() const;
    double getMinZ() const;
    double getMaxZ() const;
    std::array<double, 2> getMinMaxX() const;
    std::array<double, 2> getMinMaxY() const;
    std::array<double, 2> getMinMaxZ() const;
    void setNumVertices(unsigned value);
    void setNumTriangleFaces(unsigned value);
    void setNumQuadrangleFaces(unsigned value);
    void setNumEdges(unsigned value);
    void setMinX(double value);
    void setMaxX(double value);
    void setMinY(double value);
    void setMaxY(double value);
    void setMinZ(double value);
    void setMaxZ(double value);
    void setMinMaxX(const std::array<double, 2> &value);
    void setMinMaxY(const std::array<double, 2> &value);
    void setMinMaxZ(const std::array<double, 2> &value);
    void addVertex(const Vertex &vertex);
    void addTriangleFace(const TriangleFace &face);
    void addQuadrangleFace(const QuadrangleFace &face);
    void addEdge(const Edge &edge);
    void updateVertexAt(unsigned i, const Vertex &vertex);

protected:
    std::vector<Vertex> m_vertices; // vertices
    std::unordered_set<TriangleFace> m_triangleFaces; // triangle faces
    std::unordered_set<QuadrangleFace> m_quadrangleFaces; // quadrangle faces
    std::unordered_set<Edge> m_edges; // edges
    unsigned m_nv; // # of vertices
    unsigned m_nft; // # of triangle faces
    unsigned m_nfq; // # of quadrangle faces
    unsigned m_ne; // # of edges
    double m_minX; // minimum x coord
    double m_maxX; // maximum x coord
    double m_minY; // minimum y coord
    double m_maxY; // maximum y coord
    double m_minZ; // minimum z coord
    double m_maxZ; // maximum z coord
};

#endif //ARBF_TEST_MESH_H
