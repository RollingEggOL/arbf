//
// Created by Ke Liu on 9/22/16.
//

#include <cstdlib>
#include <array>
#include <vector>
#include <unordered_set>
#include "../include/Common.h"
#include "../include/HexMesh.h"

using namespace std;

HexMesh::HexMesh(): Mesh(), m_nh(0) {
    // nothing here
}

HexMesh::HexMesh(const HexMesh &other): Mesh(other) {
    m_hexahedrons = other.m_hexahedrons;
    m_nh = other.m_nh;
}

HexMesh::~HexMesh() {
    // nothing here
}

HexMesh & HexMesh::operator=(const HexMesh &other) {
    if (this == &other) {
        return *this;
    }

    Mesh::operator=(other);
    m_hexahedrons = other.m_hexahedrons;
    m_nh = other.m_nh;
    return *this;
}

unsigned HexMesh::getNumHexahedrons() const {
    return m_nh;
}

void HexMesh::setNumHexahedrons(unsigned value) {
    m_nh = value;
}

const std::vector<Hexahedron> &HexMesh::getHexahedrons() const {
    return m_hexahedrons;
}

void HexMesh::addHexahedron(const Hexahedron &hexahedron) {
    m_hexahedrons.push_back(hexahedron);
}

void HexMesh::exportToFile(const char *filename) {
    printf("Writing mesh file %s\n", filename);
    FILE *fout = nullptr;

    if ((fout = fopen(filename, "w")) == nullptr) {
        fprintf(stderr, "Error: Writing mesh file %s \n", filename);
        exit(EXIT_FAILURE);
    }

    fputs("HEX\n", fout);
    fprintf(fout, "%d %d %d\n", m_nv, m_nfq, m_ne);

    for (int v = 0; v < m_nv; v++) {
        fprintf(fout, "%lf %lf %lf\n", m_vertices[v].x, m_vertices[v].y, m_vertices[v].z);
    }

    for (int h = 0; h < m_nh; h++) {
        fprintf(fout, "8 %d %d %d %d %d %d %d %d\n", m_hexahedrons[h].a, m_hexahedrons[h].b, m_hexahedrons[h].c,
                m_hexahedrons[h].d, m_hexahedrons[h].e, m_hexahedrons[h].f, m_hexahedrons[h].g, m_hexahedrons[h].f);
    }

    fclose(fout);
}

unsigned HexMesh::getNumTriangleFaces() const {
    throw new NotImplementedException();
}

const std::unordered_set<TriangleFace> &HexMesh::getTriangleFaces() const {
    throw new NotImplementedException();
}

const std::vector<TriangleFace> HexMesh::getTriangleFacesAsList() const {
    throw new NotImplementedException();
}
