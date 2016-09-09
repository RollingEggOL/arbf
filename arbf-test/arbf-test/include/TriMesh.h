//
// Created by Ke Liu on 9/4/16.
//

#ifndef ARBF_TEST_TRIMESH_H
#define ARBF_TEST_TRIMESH_H

#include "Mesh.h"

class TriMesh: public Mesh {
public:
    TriMesh();
    TriMesh(const TriMesh &other);
    virtual ~TriMesh();
    TriMesh &operator=(const TriMesh &other);
    virtual void exportToFile(const char *filename);
    void deform(); // randomize input to make deformed mesh
};

#endif //ARBF_TEST_TRIMESH_H
