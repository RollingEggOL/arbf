//
// Created by Ke Liu on 9/4/16.
//

#include <cstdlib>
#include <ctime>
#include <vector>
#include <unordered_set>
#include "../include/Common.h"
#include "../include/TriMesh.h"

TriMesh::TriMesh(): Mesh() {
    // nothing here
}

TriMesh::TriMesh(const TriMesh &other): Mesh(other) {
    // nothing here
}

TriMesh &TriMesh::operator=(const TriMesh &other) {
    if (this == &other) {
        return *this;
    }

    Mesh::operator=(other);
    return *this;
}

TriMesh::~TriMesh() {
    // nothing here
}

void TriMesh::exportToFile(const char *filename) {
    printf("Writing mesh file %s\n", filename);
    FILE *fout = nullptr;

    if ((fout = fopen(filename, "w")) == nullptr) {
        fprintf(stderr, "Error: Writing mesh file %s \n", filename);
        exit(EXIT_FAILURE);
    }

    fputs("OFF\n", fout);
    fprintf(fout, "%d %d %d\n", m_nv, m_nft, m_ne);

    for (int v = 0; v < m_nv; v++) {
        fprintf(fout, "%lf %lf %lf\n", m_vertices[v].x, m_vertices[v].y, m_vertices[v].z);
    }

    auto faces = getTriangleFacesAsList();
    for (int f = 0; f < m_nft; f++) {
        fprintf(fout, "3 %d %d %d\n", faces[f].a, faces[f].b, faces[f].c);
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
    m_vertices[1].x = m_vertices[0].x + t1 * (m_vertices[3].x - m_vertices[0].x);
    m_vertices[1].y = m_vertices[0].y + t1 * (m_vertices[3].y - m_vertices[0].y);
    m_vertices[4].x = m_vertices[3].x + t2 * (m_vertices[5].x - m_vertices[3].x);
    m_vertices[4].y = m_vertices[3].y + t2 * (m_vertices[5].y - m_vertices[3].y);
    m_vertices[2].x = m_vertices[5].x + t3 * (m_vertices[0].x - m_vertices[5].x);
    m_vertices[2].y = m_vertices[5].y + t3 * (m_vertices[0].y - m_vertices[5].y);
}
