//
// Created by Ke Liu on 9/22/16.
//

#ifndef ARBF_TEST_HEXMESH_H
#define ARBF_TEST_HEXMESH_H

#include <vector>
#include "Common.h"
#include "Mesh.h"

class HexMesh: public Mesh {
public:
    HexMesh();
    HexMesh(const HexMesh &other);
    virtual ~HexMesh();
    HexMesh &operator=(const HexMesh &other);
    unsigned getNumHexahedrons() const;
    void setNumHexahedrons(unsigned value);
    const std::vector<Hexahedron> &getHexahedrons() const;
    void addHexahedron(const Hexahedron &hexahedron);
    virtual void exportToFile(const char *filename);
    virtual void update();

private:
    std::vector<Hexahedron> m_hexahedrons; // hexahedrons
    unsigned m_nh; // # of hexahedrons
};

#endif //ARBF_TEST_HEXMESH_H
