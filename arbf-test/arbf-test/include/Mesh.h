//
//  TriMesh.h
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#ifndef __arbf_test_TriMesh_h__
#define __arbf_test_TriMesh_h__

#include <vector>
#include <array>
#include <unordered_set>
#include "Common.h"

class Mesh {
public:
    virtual ~Mesh() {}

    virtual unsigned getNumVertices() const = 0;
    virtual unsigned getNumFaces() const = 0;
    virtual unsigned getNumEdges() const = 0;
    virtual const std::vector<Vertex> &getVertices() const = 0;
    virtual const std::vector<Face> &getFaces() const = 0;
    virtual const std::unordered_set<Edge, edgeHasher, edgeComparator> &getEdges() const = 0;
    virtual double getMinX() const = 0;
    virtual double getMaxX() const = 0;
    virtual double getMinY() const = 0;
    virtual double getMaxY() const = 0;
    virtual double getMinZ() const = 0;
    virtual double getMaxZ() const = 0;
    virtual std::array<double, 2> getMinMaxX() const = 0;
    virtual std::array<double, 2> getMinMaxY() const = 0;
    virtual std::array<double, 2> getMinMaxZ() const = 0;

    virtual void setNumVertices(unsigned value) = 0;
    virtual void setNumFaces(unsigned value) = 0;
    virtual void setNumEdges(unsigned value) = 0;
    virtual void setMinX(double value) = 0;
    virtual void setMaxX(double value) = 0;
    virtual void setMinY(double value) = 0;
    virtual void setMaxY(double value) = 0;
    virtual void setMinZ(double value) = 0;
    virtual void setMaxZ(double value) = 0;
    virtual void setMinMaxX(const std::array<double, 2> &value) = 0;
    virtual void setMinMaxY(const std::array<double, 2> &value) = 0;
    virtual void setMinMaxZ(const std::array<double, 2> &value) = 0;

    virtual void addVertex(const Vertex &vertex) = 0;
    virtual void addFace(const Face &face) = 0;
    virtual void addEdge(const Edge &edge) = 0;
};

class TriMesh: public Mesh {
public:
    TriMesh();
    TriMesh(const TriMesh &other);
    virtual ~TriMesh();
    TriMesh &operator=(const TriMesh &other);

    virtual unsigned getNumVertices() const;
    virtual unsigned getNumFaces() const;
    virtual unsigned getNumEdges() const;
    virtual const std::vector<Vertex> &getVertices() const;
    virtual const std::vector<Face> &getFaces() const;
    virtual const std::unordered_set<Edge, edgeHasher, edgeComparator> &getEdges() const;
    virtual double getMinX() const;
    virtual double getMaxX() const;
    virtual double getMinY() const;
    virtual double getMaxY() const;
    virtual double getMinZ() const;
    virtual double getMaxZ() const;
    virtual std::array<double, 2> getMinMaxX() const;
    virtual std::array<double, 2> getMinMaxY() const;
    virtual std::array<double, 2> getMinMaxZ() const;

    virtual void setNumVertices(unsigned value);
    virtual void setNumFaces(unsigned value);
    virtual void setNumEdges(unsigned value);
    virtual void setMinX(double value);
    virtual void setMaxX(double value);
    virtual void setMinY(double value);
    virtual void setMaxY(double value);
    virtual void setMinZ(double value);
    virtual void setMaxZ(double value);
    virtual void setMinMaxX(const std::array<double, 2> &value);
    virtual void setMinMaxY(const std::array<double, 2> &value);
    virtual void setMinMaxZ(const std::array<double, 2> &value);

    virtual void addVertex(const Vertex &vertex);
    virtual void addFace(const Face &face);
    virtual void addEdge(const Edge &edge);
    
private:
    void write(const char *path);
    void deform(); // randomize input to make deformed mesh
    
private:
    std::vector<Vertex> m_vertex; // vertices
    std::vector<Face> m_face; // faces
    std::unordered_set<Edge, edgeHasher, edgeComparator> m_edge; // edges
    unsigned m_nv; // # of vertices
    unsigned m_nf; // # of faces
    unsigned m_ne; // # of edges
    double m_minX; // minimum x coord
    double m_maxX; // maximum x coord
    double m_minY; // minimum y coord
    double m_maxY; // maximum y coord
    double m_minZ; // minimum z coord
    double m_maxZ; // maximum z coord
};

#endif /* defined(__arbf_test_TriMesh_h__) */
