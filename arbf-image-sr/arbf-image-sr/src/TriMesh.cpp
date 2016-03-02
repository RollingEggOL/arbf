//
//  TriMesh.cpp
//  arbf-image-super_resolution
//
//  Created by Ke Liu on 2/26/15.
//  Copyright (c) 2015 Ke Liu. All rights reserved.
//

#include <cstdio>
#include <cstring>
#include <cassert>
//#include <glog/logging.h>
//#include <boost/format.hpp>
#include "../include/TriMesh.h"

using namespace std;
//using namespace boost;

const int TriMesh::BUF_SIZE = 256;

TriMesh::TriMesh() : m_vertex(nullptr), m_face(nullptr), m_nv(0), m_nf(0) {}

TriMesh::TriMesh(const char *path)
{
    assert(path != nullptr);
    read(path);
}

TriMesh::TriMesh(const TriMesh &other)
{
    m_vertex = new DBL7VECT[other.m_nv];
    memcpy(m_vertex, other.m_vertex, sizeof(DBL7VECT) * other.m_nv);
    m_face = new INT3VECT[other.m_nf];
    memcpy(m_face, other.m_face, sizeof(INT3VECT) * other.m_nf);
    m_nv = other.m_nv;
    m_nf = other.m_nf;
}

TriMesh::~TriMesh()
{
    delete [] m_vertex;
    m_vertex = nullptr;
    delete [] m_face;
    m_face = nullptr;
}

TriMesh& TriMesh::operator=(const TriMesh &other)
{
    if (this == &other)
        return *this;
    
    delete [] m_vertex;
    delete [] m_face;
    m_vertex = new DBL7VECT[other.m_nv];
    memcpy(m_vertex, other.m_vertex, sizeof(DBL7VECT) * other.m_nv);
    m_face = new INT3VECT[other.m_nf];
    memcpy(m_face, other.m_face, sizeof(INT3VECT) * other.m_nf);
    m_nv = other.m_nv;
    m_nf = other.m_nf;
    
    return *this;
}

int TriMesh::numVertices() const
{
    return m_nv;
}

void TriMesh::setNumVertices(int nv)
{
    m_nv = nv;
}

int TriMesh::numFaces() const
{
    return m_nf;
}

void TriMesh::setNumFaces(int nf)
{
    m_nf = nf;
}

DBL7VECT* TriMesh::getVertices()
{
    return m_vertex;
}

const DBL7VECT* TriMesh::getVertices() const
{
    return m_vertex;
}

INT3VECT* TriMesh::getFaces()
{
    return m_face;
}

const INT3VECT* TriMesh::getFaces() const
{
    return m_face;
}

int TriMesh::read(const char *path)
{
    printf("Reading mesh file %s\n", path);
    int n, m, dim;
    int a, b, c, d;
    double x, y, z, lambda1, lambda2, vx1, vy1, gradx, grady;
    char line[BUF_SIZE];
    FILE *fin = nullptr;
    
#ifdef _WIN32
    if (fopen_s(&fin, path, "r")) {
#else
    if ((fin = fopen(path, "r")) == nullptr) {
#endif
        fprintf(stderr, "Error: Reading mesh file %s \n", path);
        return 0;
    };
    
    /* TRI format */
    while (fgets(line, BUF_SIZE, fin) != nullptr) {
        if (line[0] == 'T' && line[1] == 'R' && line[2] == 'I') {
            break;
        }
        
        if (line[0] == '#') {
            break;
        }
    }
    
#ifdef _WIN32
    fscanf_s(fin, "%d\n", &dim);
#else
    fscanf(fin, "%d\n", &dim);
#endif
    if(dim != 2) { // support 2D TRI
        fprintf(stderr, "ERROR: corrupted TRI format ...\n");
        return 0;
    }
    
#ifdef _WIN32
    fscanf_s(fin, "%d %d\n", &m, &n); // # of vertices, # of faces
#else
    fscanf(fin, "%d %d\n", &m, &n); // # of vertices, # of faces
#endif
    m_nv = m;
    m_nf = n;
    m_vertex = new DBL7VECT[m_nv];
    m_face = new INT3VECT[m_nf];
    
    /* read vertices */
    for (n = 0; n < m_nv; n++) {
#ifdef _WIN32
    fscanf_s(fin, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
             &x, &y, &lambda1, &lambda2, &vx1, &vy1, &gradx, &grady);
#else
    fscanf(fin, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
           &x, &y, &lambda1, &lambda2, &vx1, &vy1, &gradx, &grady);
#endif
        m_vertex[n].x = x;
        m_vertex[n].y = y;
        m_vertex[n].intensity = 0;
        m_vertex[n].lambda1 = lambda1;
        m_vertex[n].lambda2 = lambda2;
        m_vertex[n].v1.x = vx1;
        m_vertex[n].v1.y = vy1;
        m_vertex[n].v2.x = -vy1;
        m_vertex[n].v2.y = vx1;
        m_vertex[n].gradient.x = gradx;
        m_vertex[n].gradient.y = grady;
    }
    
    /* read faces */
    for (n = 0; n < m_nf; n++) {
#ifdef _WIN32
        fscanf_s(fin, "%d %d %d %d %lf\n", &a, &b, &c, &d, &z);
#else
        fscanf(fin,"%d %d %d %d %lf\n",&a, &b, &c, &d, &z);
#endif
        m_face[n].a = b;
        m_face[n].b = c;
        m_face[n].c = d;
        m_face[n].intensity = z;
    }
    fclose(fin);
        
    return 1;
}
