//
//  TriMesh.cpp
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#include <cstdio>
#include <cstring>
#include <cassert>
#include <limits>
#include "../include/TriMesh.h"

using namespace std;
//using namespace boost;

TriMesh::TriMesh() : m_vertex(nullptr), m_face(nullptr), m_nv(0), m_nf(0), m_ne(0),
m_minX(numeric_limits<double>::max()), m_maxX(numeric_limits<double>::lowest()),
m_minY(numeric_limits<double>::max()), m_maxY(numeric_limits<double>::lowest()),
m_minZ(numeric_limits<double>::max()), m_maxZ(numeric_limits<double>::lowest()) {}

TriMesh::TriMesh(const char *path): m_nv(0), m_nf(0), m_ne(0),
m_minX(numeric_limits<double>::max()), m_maxX(numeric_limits<double>::lowest()),
m_minY(numeric_limits<double>::max()), m_maxY(numeric_limits<double>::lowest()),
m_minZ(numeric_limits<double>::max()), m_maxZ(numeric_limits<double>::lowest())
{
    assert(path != nullptr);
    read(path);
}

TriMesh::TriMesh(const TriMesh &other) {
    m_vertex = new Vertex[other.m_nv];
    memcpy(m_vertex, other.m_vertex, sizeof(Vertex) * other.m_nv);
    m_face = new Face[other.m_nf];
    memcpy(m_face, other.m_face, sizeof(Face) * other.m_nf);
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

TriMesh::~TriMesh() {
    delete [] m_vertex;
    m_vertex = nullptr;
    delete [] m_face;
    m_face = nullptr;
}

TriMesh& TriMesh::operator=(const TriMesh &other) {
    if (this == &other)
        return *this;
    
    this->~TriMesh();
    m_vertex = new Vertex[other.m_nv];
    memcpy(m_vertex, other.m_vertex, sizeof(Vertex) * other.m_nv);
    m_face = new Face[other.m_nf];
    memcpy(m_face, other.m_face, sizeof(Face) * other.m_nf);
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

void TriMesh::read(const char *path) {
    printf("Reading mesh file %s\n", path);
    int dim;
    int a, b, c, d;
    double x, y, z;
    char line[256];
    FILE *fin = nullptr;
    
    if ((fin = fopen(path, "r")) == nullptr) {
        fprintf(stderr, "Error: Reading mesh file %s \n", path);
        exit(EXIT_FAILURE);
    };
    
    // TRI format
    while (fgets(line, 256, fin) != nullptr) {
        if (line[0] == 'T' && line[1] == 'R' && line[2] == 'I') {
            break;
        }
        
        if (line[0] == '#') {
            break;
        }
    }
    
    fscanf(fin, "%d\n", &dim);

    // support 2D and 3D TRI file
    if (dim != 2 && dim != 3) {
        fprintf(stderr, "ERROR: corrupted TRI format ...\n");
        exit(EXIT_FAILURE);
    }
    
    fscanf(fin, "%d %d %d\n", &m_nv, &m_nf, &m_ne); // # of vertices, # of faces, # of edges
    m_vertex = new Vertex[m_nv];
    m_face = new Face[m_nf];
    
    // read vertices
    for (int v = 0; v < m_nv; v++) {
        fscanf(fin, "%lf %lf %lf\n", &x, &y, &z);
        m_vertex[v].x = x;
        m_vertex[v].y = y;
        m_vertex[v].z = z;
        m_vertex[v].intensity = 1;
        m_vertex[v].eig1 = 0;
        m_vertex[v].eig2 = 0;
        m_vertex[v].eig3 = 0;
        m_vertex[v].eigVec1 = {1, 0, 0};
        m_vertex[v].eigVec2 = {0, 1, 0};
        m_vertex[v].eigVec3 = {0, 0, 1};
        
        if (x < m_minX) {
            m_minX = x;
        }
        
        if (x > m_maxX) {
            m_maxX = x;
        }
        
        if (y < m_minY) {
            m_minY = y;
        }
        
        if (y > m_maxY) {
            m_maxY = y;
        }
        
        if (z < m_minZ) {
            m_minZ = z;
        }
        
        if (z > m_maxZ) {
            m_maxZ = z;
        }
    }
    
    // read faces
    for (int f = 0; f < m_nf; f++) {
        fscanf(fin,"%d %d %d %d\n",&a, &b, &c, &d);
        m_face[f].a = b;
        m_face[f].b = c;
        m_face[f].c = d;
        m_face[f].intensity = -1;
        
        Edge e1, e2, e3;
        e1.a = m_face[f].a;
        e1.b = m_face[f].b;
        e1.intensity = -1;
        e2.a = m_face[f].b;
        e2.b = m_face[f].c;
        e2.intensity = -1;
        e3.a = m_face[f].c;
        e3.b = m_face[f].a;
        e3.intensity = -1;
        m_edge.insert(e1);
        m_edge.insert(e2);
        m_edge.insert(e3);
    }
    fclose(fin);
    assert(m_edge.size() == m_ne);
}
