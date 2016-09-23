//
// Created by Ke Liu on 8/3/16.
//

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <array>
#include <unordered_set>
#include "../include/Mesh.h"
#include "../include/TriMesh.h"
#include "../include/TetMesh.h"
#include "../include/MeshFactory.h"

using namespace std;

// MeshFactory implementation
MeshFactory::~MeshFactory() {}


// TriMeshFactory implementation
TriMeshFactory::~TriMeshFactory() {}

unique_ptr<Mesh> TriMeshFactory::createMeshFromFile(const char* filename) {
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

        array<double, 2> minMaxX = mesh->getMinMaxX();
        array<double, 2> minMaxY = mesh->getMinMaxY();
        array<double, 2> minMaxZ = mesh->getMinMaxZ();
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
    return unique_ptr<Mesh>(mesh);
}


// TetMeshFactory implementation
TetMeshFactory::~TetMeshFactory() {}

unique_ptr<Mesh> TetMeshFactory::createMeshFromFile(const char *filename) {
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

        if (line[0] == 'T' && line[1] == 'E' && line[2] == 'T') {
            break;
        }
    }

    fscanf(fin, "%d\n", &dim);

    // support 3D TET file
    if (dim != 3) {
        fprintf(stderr, "ERROR: corrupted TRI format ...\n");
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

        array<double, 2> minMaxX = mesh->getMinMaxX();
        array<double, 2> minMaxY = mesh->getMinMaxY();
        array<double, 2> minMaxZ = mesh->getMinMaxZ();
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
    return unique_ptr<Mesh>(mesh);
}


// MeshFactoryWithExperimentalFeatures implementation
MeshFactoryWithExperimentalFeatures::~MeshFactoryWithExperimentalFeatures() {}

unique_ptr<Mesh> MeshFactoryWithExperimentalFeatures::createMeshFromFile(const char* filename) {
    if (m_factory) {
        return m_factory->createMeshFromFile(filename);
    } else {
        return nullptr;
    }
}

void MeshFactoryWithExperimentalFeatures::setMeshFactory(unique_ptr<MeshFactory> factory) {
    m_factory = std::move(factory);
}


// MeshFactoryWithEnlargedBoundingBox implementation
unique_ptr<Mesh> MeshFactoryWithEnlargedBoundingBox::createMeshFromFile(const char* filename) {
    unique_ptr<Mesh> mesh = MeshFactoryWithExperimentalFeatures::createMeshFromFile(filename);
    enlargeBoundingBox(mesh.get());
    return std::move(mesh);
}

void MeshFactoryWithEnlargedBoundingBox::enlargeBoundingBox(Mesh *mesh) {
    array<double, 2> minMaxX = mesh->getMinMaxX();
    array<double, 2> minMaxY = mesh->getMinMaxY();
    mesh->setMinMaxX({ minMaxX[0]-0.1, minMaxX[1]+0.1 });
    mesh->setMinMaxY({ minMaxY[0]-0.1, minMaxY[1]+0.1 });
}
