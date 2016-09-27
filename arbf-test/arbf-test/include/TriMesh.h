//
// Created by Ke Liu on 9/4/16.
//

#ifndef ARBF_TEST_TRIMESH_H
#define ARBF_TEST_TRIMESH_H

#include <vector>
#include <unordered_set>
#include "Common.h"
#include "Mesh.h"

class TriMesh: public Mesh {
public:
    TriMesh();
    TriMesh(const TriMesh &other);
    virtual ~TriMesh();
    TriMesh &operator=(const TriMesh &other);
    virtual void exportToFile(const char *filename);
    virtual unsigned getNumQuadrangleFaces() const;
    virtual const std::unordered_set<QuadrangleFace> &getQuadrangleFaces() const;
    virtual const std::vector<QuadrangleFace> getQuadrangleFacesAsList() const;
    void deform(); // randomize input to make deformed mesh
};

#endif //ARBF_TEST_TRIMESH_H
