//
// Created by Ke Liu on 9/4/16.
//

#ifndef ARBF_TEST_TETMESH_H
#define ARBF_TEST_TETMESH_H

#include <vector>
#include <unordered_set>
#include "Common.h"
#include "Mesh.h"

class TetMesh: public Mesh {
public:
    TetMesh();
    TetMesh(const TetMesh &other);
    virtual ~TetMesh();
    TetMesh &operator=(const TetMesh &other);
    virtual void exportToFile(const char *filename);
    virtual void update();
    unsigned getNumTetrahedrons() const;
    void setNumTetrahedrons(unsigned value);
    const std::vector<Tetrahedron> &getTetrahedrons() const;
    void addTetrahedron(const Tetrahedron &tetrahedron);

private:
    std::vector<Tetrahedron> m_tetrahedrons; // tetrahedrons
    unsigned m_nt; // # of tetrahedrons
};

#endif //ARBF_TEST_TETMESH_H
