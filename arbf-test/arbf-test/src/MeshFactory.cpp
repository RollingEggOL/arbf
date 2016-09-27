//
// Created by Ke Liu on 8/3/16.
//

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <array>
#include <unordered_set>
#include "../include/Common.h"
#include "../include/Mesh.h"
#include "../include/TriMesh.h"
#include "../include/TetMesh.h"
#include "../include/HexMesh.h"
#include "../include/MeshFactory.h"

// MeshFactory implementation
MeshFactory::~MeshFactory() {}


// TriMeshFactory implementation
TriMeshFactory::~TriMeshFactory() {}

std::unique_ptr<Mesh> TriMeshFactory::createMeshFromFile(const char* filename) {
    printf("Reading mesh file %s\n", filename);
    unsigned a, b, c, d;
    unsigned dim;
    double x, y, z;
    char line[256];
    FILE *fin = nullptr;

    if ((fin = fopen(filename, "r")) == nullptr) {
        fprintf(stderr, "Error: Reading mesh file %s \n", filename);
        exit(EXIT_FAILURE);
    }

    // TRI format
    while (fgets(line, 256, fin) != nullptr) {
        if (line[0] == '#') {
            continue;
        }

        if (line[0] == 'T' && line[1] == 'R' && line[2] == 'I') {
            break;
        }
    }

    fscanf(fin, "%d\n", &dim);

    // support 2D and 3D TRI file
    if (dim != 2 && dim != 3) {
        fprintf(stderr, "ERROR: corrupted TRI format ...\n");
        exit(EXIT_FAILURE);
    }

    Mesh* mesh = new TriMesh();

    fscanf(fin, "%d %d %d\n", &a, &b, &c); // # of vertices, # of faces, # of edges
    mesh->setNumVertices(a);
    mesh->setNumTriangleFaces(b);
    mesh->setNumEdges(c);

    // read vertices
    for (int v = 0; v < mesh->getNumVertices(); v++) {
        fscanf(fin, "%lf %lf %lf\n", &x, &y, &z);
        mesh->addVertex(Vertex(x, y, z, 1.0));

        std::array<double, 2> minMaxX = mesh->getMinMaxX();
        std::array<double, 2> minMaxY = mesh->getMinMaxY();
        std::array<double, 2> minMaxZ = mesh->getMinMaxZ();
        if (x < minMaxX[0]) {
            minMaxX[0] = x;
        }

        if (x > minMaxX[1]) {
            minMaxX[1] = x;
        }

        if (y < minMaxY[0]) {
            minMaxY[0] = y;
        }

        if (y > minMaxY[1]) {
            minMaxY[1] = y;
        }

        if (z < minMaxZ[0]) {
            minMaxZ[0] = z;
        }

        if (z > minMaxZ[1]) {
            minMaxZ[1] = z;
        }
        mesh->setMinMaxX(minMaxX);
        mesh->setMinMaxY(minMaxY);
        mesh->setMinMaxZ(minMaxZ);
    }

    // read faces
    for (int f = 0; f < mesh->getNumTriangleFaces(); f++) {
        fscanf(fin,"%d %d %d %d\n",&a, &b, &c, &d);
        assert(a == 3);
        TriangleFace face(b, c, d, -1.0);
        const Vertex *v0 = &mesh->getVertices()[face.a];
        const Vertex *v1 = &mesh->getVertices()[face.b];
        const Vertex *v2 = &mesh->getVertices()[face.c];
        face.computeCenter(*v0, *v1, *v2);
        mesh->addTriangleFace(face);

        Edge e1(face.a, face.b, -1.0);
        Edge e2(face.b, face.c, -1.0);
        Edge e3(face.c, face.a, -1.0);
        v0 = &mesh->getVertices()[e1.a];
        v1 = &mesh->getVertices()[e1.b];
        e1.computeCenter(*v0, *v1);
        mesh->addEdge(e1);
        mesh->addEdge(e2);
        mesh->addEdge(e3);
    }

    fclose(fin);
    assert(mesh->getNumEdges() == mesh->getEdges().size());
    return std::unique_ptr<Mesh>(mesh);
}


// TetMeshFactory implementation
TetMeshFactory::~TetMeshFactory() {}

std::unique_ptr<Mesh> TetMeshFactory::createMeshFromFile(const char *filename) {
    printf("Reading mesh file %s\n", filename);
    unsigned a, b, c, d;
    unsigned dim;
    double x, y, z;
    char line[256];
    FILE *fin = nullptr;

    if ((fin = fopen(filename, "r")) == nullptr) {
        fprintf(stderr, "Error: Reading mesh file %s \n", filename);
        exit(EXIT_FAILURE);
    }

    // TET format
    while (fgets(line, 256, fin) != nullptr) {
        if (line[0] == '#') {
            continue;
        }

        if (line[0] == 'T' && line[1] == 'E' && line[2] == 'T') {
            break;
        }
    }

    fscanf(fin, "%d\n", &dim);

    // support 3D TET file
    if (dim != 3) {
        fprintf(stderr, "ERROR: corrupted TET format ...\n");
        exit(EXIT_FAILURE);
    }

    TetMesh* mesh = new TetMesh();

    fscanf(fin, "%d %d\n", &a, &b); // # of vertices, # of tetrahedrons
    mesh->setNumVertices(a);
    mesh->setNumTetrahedrons(b);

    // read vertices
    for (int v = 0; v < mesh->getNumVertices(); v++) {
        fscanf(fin, "%lf %lf %lf\n", &x, &y, &z);
        mesh->addVertex(Vertex(x, y, z, 1.0));

        std::array<double, 2> minMaxX = mesh->getMinMaxX();
        std::array<double, 2> minMaxY = mesh->getMinMaxY();
        std::array<double, 2> minMaxZ = mesh->getMinMaxZ();
        if (x < minMaxX[0]) {
            minMaxX[0] = x;
        }

        if (x > minMaxX[1]) {
            minMaxX[1] = x;
        }

        if (y < minMaxY[0]) {
            minMaxY[0] = y;
        }

        if (y > minMaxY[1]) {
            minMaxY[1] = y;
        }

        if (z < minMaxZ[0]) {
            minMaxZ[0] = z;
        }

        if (z > minMaxZ[1]) {
            minMaxZ[1] = z;
        }
        mesh->setMinMaxX(minMaxX);
        mesh->setMinMaxY(minMaxY);
        mesh->setMinMaxZ(minMaxZ);
    }

    // read tetrahedrons
    std::unordered_set<TriangleFace> set;
    for (int t = 0; t < mesh->getNumTetrahedrons(); t++) {
        fscanf(fin,"4 %d %d %d %d\n",&a, &b, &c, &d);
        Tetrahedron tet(a, b, c, d, -1.0);
        const Vertex *v0 = &mesh->getVertices()[tet.a];
        const Vertex *v1 = &mesh->getVertices()[tet.b];
        const Vertex *v2 = &mesh->getVertices()[tet.c];
        const Vertex *v3 = &mesh->getVertices()[tet.d];
        tet.computeCenter(*v0, *v1, *v2, *v3);

        v0 = &mesh->getVertices()[tet.f1.a];
        v1 = &mesh->getVertices()[tet.f1.b];
        v2 = &mesh->getVertices()[tet.f1.c];
        tet.f1.computeCenter(*v0, *v1, *v2);
        set.insert(tet.f1);

        v0 = &mesh->getVertices()[tet.f2.a];
        v1 = &mesh->getVertices()[tet.f2.b];
        v2 = &mesh->getVertices()[tet.f2.c];
        tet.f2.computeCenter(*v0, *v1, *v2);
        set.insert(tet.f2);

        v0 = &mesh->getVertices()[tet.f3.a];
        v1 = &mesh->getVertices()[tet.f3.b];
        v2 = &mesh->getVertices()[tet.f3.c];
        tet.f3.computeCenter(*v0, *v1, *v2);
        set.insert(tet.f3);

        v0 = &mesh->getVertices()[tet.f4.a];
        v1 = &mesh->getVertices()[tet.f4.b];
        v2 = &mesh->getVertices()[tet.f4.c];
        tet.f4.computeCenter(*v0, *v1, *v2);
        set.insert(tet.f4);
        mesh->addTetrahedron(tet);

        Edge e1(tet.a, tet.b, -1.0);
        v0 = &mesh->getVertices()[e1.a];
        v1 = &mesh->getVertices()[e1.b];
        e1.computeCenter(*v0, *v1);
        mesh->addEdge(e1);

        Edge e2(tet.a, tet.c, -1.0);
        v0 = &mesh->getVertices()[e2.a];
        v1 = &mesh->getVertices()[e2.b];
        e2.computeCenter(*v0, *v1);
        mesh->addEdge(e2);

        Edge e3(tet.a, tet.d, -1.0);
        v0 = &mesh->getVertices()[e3.a];
        v1 = &mesh->getVertices()[e3.b];
        e3.computeCenter(*v0, *v1);
        mesh->addEdge(e3);

        Edge e4(tet.b, tet.c, -1.0);
        v0 = &mesh->getVertices()[e4.a];
        v1 = &mesh->getVertices()[e4.b];
        e4.computeCenter(*v0, *v1);
        mesh->addEdge(e4);

        Edge e5(tet.b, tet.d, -1.0);
        v0 = &mesh->getVertices()[e5.a];
        v1 = &mesh->getVertices()[e5.b];
        e5.computeCenter(*v0, *v1);
        mesh->addEdge(e5);

        Edge e6(tet.c, tet.d, -1.0);
        v0 = &mesh->getVertices()[e6.a];
        v1 = &mesh->getVertices()[e6.b];
        e6.computeCenter(*v0, *v1);
        mesh->addEdge(e6);
    }

    for (auto &f: set) {
        mesh->addTriangleFace(f);
    }

    mesh->setNumTriangleFaces(mesh->getTriangleFaces().size());
    mesh->setNumEdges(mesh->getEdges().size());
    fclose(fin);
    return std::unique_ptr<Mesh>(mesh);
}


// HexMeshFactory implementation
HexMeshFactory::~HexMeshFactory() {}

std::unique_ptr<Mesh> HexMeshFactory::createMeshFromFile(const char *filename) {
    printf("Reading mesh file %s\n", filename);
    unsigned a, b, c, d, e, f, g, h;
    unsigned dim;
    double x, y, z;
    char line[256];
    FILE *fin = nullptr;

    if ((fin = fopen(filename, "r")) == nullptr) {
        fprintf(stderr, "Error: Reading mesh file %s \n", filename);
        exit(EXIT_FAILURE);
    }

    // HEX format
    while (fgets(line, 256, fin) != nullptr) {
        if (line[0] == '#') {
            continue;
        }

        if (line[0] == 'H' && line[1] == 'E' && line[2] == 'X') {
            break;
        }
    }

    fscanf(fin, "%d\n", &dim);

    // support 3D TET file
    if (dim != 3) {
        fprintf(stderr, "ERROR: corrupted HEX format ...\n");
        exit(EXIT_FAILURE);
    }

    HexMesh* mesh = new HexMesh();

    fscanf(fin, "%d %d\n", &a, &b); // # of vertices, # of hexahedrons
    mesh->setNumVertices(a);
    mesh->setNumHexahedrons(b);

    // read vertices
    for (int v = 0; v < mesh->getNumVertices(); v++) {
        fscanf(fin, "%lf %lf %lf\n", &x, &y, &z);
        mesh->addVertex(Vertex(x, y, z, 1.0));

        std::array<double, 2> minMaxX = mesh->getMinMaxX();
        std::array<double, 2> minMaxY = mesh->getMinMaxY();
        std::array<double, 2> minMaxZ = mesh->getMinMaxZ();
        if (x < minMaxX[0]) {
            minMaxX[0] = x;
        }

        if (x > minMaxX[1]) {
            minMaxX[1] = x;
        }

        if (y < minMaxY[0]) {
            minMaxY[0] = y;
        }

        if (y > minMaxY[1]) {
            minMaxY[1] = y;
        }

        if (z < minMaxZ[0]) {
            minMaxZ[0] = z;
        }

        if (z > minMaxZ[1]) {
            minMaxZ[1] = z;
        }
        mesh->setMinMaxX(minMaxX);
        mesh->setMinMaxY(minMaxY);
        mesh->setMinMaxZ(minMaxZ);
    }

    // read hexahedrons
    std::unordered_set<QuadrangleFace> set;
    for (int t = 0; t < mesh->getNumHexahedrons(); t++) {
        fscanf(fin,"8 %d %d %d %d %d %d %d %d\n",&a, &b, &c, &d, &e, &f, &g, &h);
        Hexahedron hex(a, b, c, d, e, f, g, h, -1.0);
        const Vertex *v0 = &mesh->getVertices()[hex.a];
        const Vertex *v1 = &mesh->getVertices()[hex.b];
        const Vertex *v2 = &mesh->getVertices()[hex.c];
        const Vertex *v3 = &mesh->getVertices()[hex.d];
        const Vertex *v4 = &mesh->getVertices()[hex.e];
        const Vertex *v5 = &mesh->getVertices()[hex.f];
        const Vertex *v6 = &mesh->getVertices()[hex.g];
        const Vertex *v7 = &mesh->getVertices()[hex.h];
        hex.computeCenter(*v0, *v1, *v2, *v3, *v4, *v5, *v6, *v7);

        v0 = &mesh->getVertices()[hex.f1.a];
        v1 = &mesh->getVertices()[hex.f1.b];
        v2 = &mesh->getVertices()[hex.f1.c];
        v3 = &mesh->getVertices()[hex.f1.d];
        hex.f1.computeCenter(*v0, *v1, *v2, *v3);
        set.insert(hex.f1);

        v0 = &mesh->getVertices()[hex.f2.a];
        v1 = &mesh->getVertices()[hex.f2.b];
        v2 = &mesh->getVertices()[hex.f2.c];
        v3 = &mesh->getVertices()[hex.f2.d];
        hex.f2.computeCenter(*v0, *v1, *v2, *v3);
        set.insert(hex.f2);

        v0 = &mesh->getVertices()[hex.f3.a];
        v1 = &mesh->getVertices()[hex.f3.b];
        v2 = &mesh->getVertices()[hex.f3.c];
        v3 = &mesh->getVertices()[hex.f3.d];
        hex.f3.computeCenter(*v0, *v1, *v2, *v3);
        set.insert(hex.f3);

        v0 = &mesh->getVertices()[hex.f4.a];
        v1 = &mesh->getVertices()[hex.f4.b];
        v2 = &mesh->getVertices()[hex.f4.c];
        v3 = &mesh->getVertices()[hex.f4.d];
        hex.f4.computeCenter(*v0, *v1, *v2, *v3);
        set.insert(hex.f4);

        v0 = &mesh->getVertices()[hex.f5.a];
        v1 = &mesh->getVertices()[hex.f5.b];
        v2 = &mesh->getVertices()[hex.f5.c];
        v3 = &mesh->getVertices()[hex.f5.d];
        hex.f5.computeCenter(*v0, *v1, *v2, *v3);
        set.insert(hex.f5);

        v0 = &mesh->getVertices()[hex.f6.a];
        v1 = &mesh->getVertices()[hex.f6.b];
        v2 = &mesh->getVertices()[hex.f6.c];
        v3 = &mesh->getVertices()[hex.f6.d];
        hex.f6.computeCenter(*v0, *v1, *v2, *v3);
        set.insert(hex.f6);
        mesh->addHexahedron(hex);

        Edge e1(hex.a, hex.b, -1.0);
        v0 = &mesh->getVertices()[e1.a];
        v1 = &mesh->getVertices()[e1.b];
        e1.computeCenter(*v0, *v1);
        mesh->addEdge(e1);

        Edge e2(hex.a, hex.d, -1.0);
        v0 = &mesh->getVertices()[e2.a];
        v1 = &mesh->getVertices()[e2.b];
        e2.computeCenter(*v0, *v1);
        mesh->addEdge(e2);

        Edge e3(hex.a, hex.e, -1.0);
        v0 = &mesh->getVertices()[e3.a];
        v1 = &mesh->getVertices()[e3.b];
        e3.computeCenter(*v0, *v1);
        mesh->addEdge(e3);

        Edge e4(hex.b, hex.c, -1.0);
        v0 = &mesh->getVertices()[e4.a];
        v1 = &mesh->getVertices()[e4.b];
        e4.computeCenter(*v0, *v1);
        mesh->addEdge(e4);

        Edge e5(hex.b, hex.f, -1.0);
        v0 = &mesh->getVertices()[e5.a];
        v1 = &mesh->getVertices()[e5.b];
        e5.computeCenter(*v0, *v1);
        mesh->addEdge(e5);

        Edge e6(hex.c, hex.d, -1.0);
        v0 = &mesh->getVertices()[e6.a];
        v1 = &mesh->getVertices()[e6.b];
        e6.computeCenter(*v0, *v1);
        mesh->addEdge(e6);

        Edge e7(hex.c, hex.g, -1.0);
        v0 = &mesh->getVertices()[e7.a];
        v1 = &mesh->getVertices()[e7.b];
        e7.computeCenter(*v0, *v1);
        mesh->addEdge(e7);

        Edge e8(hex.d, hex.h, -1.0);
        v0 = &mesh->getVertices()[e8.a];
        v1 = &mesh->getVertices()[e8.b];
        e8.computeCenter(*v0, *v1);
        mesh->addEdge(e8);

        Edge e9(hex.e, hex.f, -1.0);
        v0 = &mesh->getVertices()[e9.a];
        v1 = &mesh->getVertices()[e9.b];
        e9.computeCenter(*v0, *v1);
        mesh->addEdge(e9);

        Edge e10(hex.e, hex.h, -1.0);
        v0 = &mesh->getVertices()[e10.a];
        v1 = &mesh->getVertices()[e10.b];
        e10.computeCenter(*v0, *v1);
        mesh->addEdge(e10);

        Edge e11(hex.f, hex.g, -1.0);
        v0 = &mesh->getVertices()[e11.a];
        v1 = &mesh->getVertices()[e11.b];
        e11.computeCenter(*v0, *v1);
        mesh->addEdge(e11);

        Edge e12(hex.g, hex.h, -1.0);
        v0 = &mesh->getVertices()[e12.a];
        v1 = &mesh->getVertices()[e12.b];
        e12.computeCenter(*v0, *v1);
        mesh->addEdge(e12);
    }

    for (auto &f: set) {
        mesh->addQuadrangleFace(f);
    }

    mesh->setNumQuadrangleFaces(mesh->getQuadrangleFaces().size());
    mesh->setNumEdges(mesh->getEdges().size());
    fclose(fin);
    return std::unique_ptr<Mesh>(mesh);
}


// MeshFactoryWithExperimentalFeatures implementation
MeshFactoryWithExperimentalFeatures::~MeshFactoryWithExperimentalFeatures() {}

std::unique_ptr<Mesh> MeshFactoryWithExperimentalFeatures::createMeshFromFile(const char* filename) {
    if (m_factory) {
        return m_factory->createMeshFromFile(filename);
    } else {
        return nullptr;
    }
}

void MeshFactoryWithExperimentalFeatures::setMeshFactory(std::unique_ptr<MeshFactory> factory) {
    m_factory = std::move(factory);
}


// MeshFactoryWithEnlargedBoundingBox implementation
std::unique_ptr<Mesh> MeshFactoryWithEnlargedBoundingBox::createMeshFromFile(const char* filename) {
    std::unique_ptr<Mesh> mesh = MeshFactoryWithExperimentalFeatures::createMeshFromFile(filename);
    enlargeBoundingBox(mesh.get());
    return std::move(mesh);
}

void MeshFactoryWithEnlargedBoundingBox::enlargeBoundingBox(Mesh *mesh) {
    std::array<double, 2> minMaxX = mesh->getMinMaxX();
    std::array<double, 2> minMaxY = mesh->getMinMaxY();
    mesh->setMinMaxX({ minMaxX[0]-0.1, minMaxX[1]+0.1 });
    mesh->setMinMaxY({ minMaxY[0]-0.1, minMaxY[1]+0.1 });
}
