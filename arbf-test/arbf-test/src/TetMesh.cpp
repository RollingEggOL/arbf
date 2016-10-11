//
// Created by Ke Liu on 9/4/16.
//

#include <cstdlib>
#include <array>
#include <vector>
#include <unordered_set>
#include "../include/Common.h"
#include "../include/TetMesh.h"

using namespace std;

TetMesh::TetMesh(): Mesh(), m_nt(0) {
    // nothing here
}

TetMesh::TetMesh(const TetMesh &other): Mesh(other) {
    m_tetrahedrons = other.m_tetrahedrons;
    m_nt = other.m_nt;
}

TetMesh::~TetMesh() {
    // nothing here
}

TetMesh & TetMesh::operator=(const TetMesh &other) {
    if (this == &other) {
        return *this;
    }

    Mesh::operator=(other);
    m_tetrahedrons = other.m_tetrahedrons;
    m_nt = other.m_nt;
    return *this;
}

unsigned TetMesh::getNumTetrahedrons() const {
    return m_nt;
}

void TetMesh::setNumTetrahedrons(unsigned value) {
    m_nt = value;
}

const std::vector<Tetrahedron> &TetMesh::getTetrahedrons() const {
    return m_tetrahedrons;
}

void TetMesh::addTetrahedron(const Tetrahedron &tetrahedron) {
    m_tetrahedrons.push_back(tetrahedron);
}

void TetMesh::exportToFile(const char *filename) {
    printf("Writing mesh file %s\n", filename);
    FILE *fout = nullptr;

    if ((fout = fopen(filename, "w")) == nullptr) {
        fprintf(stderr, "Error: Writing mesh file %s \n", filename);
        exit(EXIT_FAILURE);
    }

    fputs("TET\n", fout);
    fprintf(fout, "%d %d %d\n", m_nv, m_nft, m_ne);

    for (int v = 0; v < m_nv; v++) {
        fprintf(fout, "%lf %lf %lf\n", m_vertices[v].x, m_vertices[v].y, m_vertices[v].z);
    }

    for (int t = 0; t < m_nt; t++) {
        fprintf(fout, "4 %d %d %d %d\n", m_tetrahedrons[t].a, m_tetrahedrons[t].b, m_tetrahedrons[t].c, m_tetrahedrons[t].d);
    }

    fclose(fout);
}

void TetMesh::update() {
    Mesh::update();

    for (auto &tet: m_tetrahedrons) {
        tet.computeCenter(m_vertices[tet.a], m_vertices[tet.b], m_vertices[tet.c], m_vertices[tet.d]);
    }
}
