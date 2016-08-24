//
// Created by Ke Liu on 8/3/16.
//

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <array>
#include "../include/Mesh.h"
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
    mesh->setNumFaces(b);
    mesh->setNumEdges(c);

    // read vertices
    for (int v = 0; v < mesh->getNumVertices(); v++) {
        fscanf(fin, "%lf %lf %lf\n", &x, &y, &z);
        Vertex vertex;
        vertex.x = x;
        vertex.y = y;
        vertex.z = z;
        vertex.intensity = 1;
        vertex.eig1 = 0;
        vertex.eig2 = 0;
        vertex.eig3 = 0;
        vertex.eigVec1 = {1, 0, 0};
        vertex.eigVec2 = {0, 1, 0};
        vertex.eigVec3 = {0, 0, 1};
        mesh->addVertex(vertex);

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

//    // enlarge bounding box temporarily
//    m_minX -= 0.1;
//    m_minY -= 0.1;
//    m_maxX += 0.1;
//    m_maxY += 0.1;

    // read faces
    for (int f = 0; f < mesh->getNumFaces(); f++) {
        fscanf(fin,"%d %d %d %d\n",&a, &b, &c, &d);
        assert(a == 3);
        Face face;
        face.a = b;
        face.b = c;
        face.c = d;
        face.intensity = -1;
        mesh->addFace(face);

        Edge e1, e2, e3;
        e1.a = face.a;
        e1.b = face.b;
        e1.intensity = -1;
        e2.a = face.b;
        e2.b = face.c;
        e2.intensity = -1;
        e3.a = face.c;
        e3.b = face.a;
        e3.intensity = -1;
        mesh->addEdge(e1);
        mesh->addEdge(e2);
        mesh->addEdge(e3);
    }

    fclose(fin);
    assert(mesh->getNumEdges() == mesh->getEdges().size());
//    string outputFileName = string(path) + ".off";
//    deform();
//    write(outputFileName.c_str());
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
