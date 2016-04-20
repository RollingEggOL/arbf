//
//  TriMesh.h
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#ifndef __arbf_test_TriMesh_h__
#define __arbf_test_TriMesh_h__

#include <string>
#include <unordered_set>
#include "Common.h"

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

class TriMesh
{
public:
    /*
     * Default constructor
     */
    TriMesh();
    
    /*
     * Create mesh by reading mesh file.
     * \param path: mesh file path
     */
    TriMesh(const char *path);
    
    /*
     * Copy constructor
     */
    TriMesh(const TriMesh &other);
    
    /*
     * Destructor
     */
    virtual ~TriMesh();
    
    /*
     * Assignment operator
     */
    TriMesh& operator=(const TriMesh &other);
    
    int getNumVertices() const {
        return m_nv;
    }
    
    void setNumVertices(int nv) {
        m_nv = nv;
    }
    
    int getNumFaces() const {
        return m_nf;
    }
    
    int getNumEdges() const {
        return m_ne;
    }
    
    void setNumFaces(int nf) {
        m_nf = nf;
    }
    
    const Vertex* getVertices() const {
        return m_vertex;
    }
    
    const Face* getFaces() const {
        return m_face;
    }
    
    const std::unordered_set<Edge, edgeHasher, edgeComparator>& getEdges() const {
        return m_edge;
    }
    
    double getMinX() const {
        return m_minX;
    }
    
    double getMaxX() const {
        return m_maxX;
    }
    
    double getMinY() const {
        return m_minY;
    }
    
    double getMaxY() const {
        return m_maxY;
    }
    
    double getMinZ() const {
        return m_minZ;
    }
    
    double getMaxZ() const {
        return m_maxZ;
    }
    
private:
    void read(const char *path);
    void write(const char *path);
    void deform(); // randomize input to make deformed mesh
    
private:
    Vertex *m_vertex; // vertices
    Face *m_face; // faces
    std::unordered_set<Edge, edgeHasher, edgeComparator> m_edge; // edges
    int m_nv; // # of vertices
    int m_nf; // # of faces
    int m_ne; // # of edges
    double m_minX; // minimum x coord
    double m_maxX; // maximum x coord
    double m_minY; // minimum y coord
    double m_maxY; // maximum y coord
    double m_minZ; // minimum z coord
    double m_maxZ; // maximum z coord
};

#endif /* defined(__arbf_test_TriMesh_h__) */
