//
// Created by Ke Liu on 8/4/16.
//

#ifndef ARBF_TEST_MESHNEIGHBORHOOD_H
#define ARBF_TEST_MESHNEIGHBORHOOD_H

#include <vector>
#include <set>
#include "Mesh.h"

/*
 * Abstract mesh neighborhood strategy
 */
class MeshNeighborhood {
public:
    virtual ~MeshNeighborhood();
    virtual void computeVertexNeighbors(const Mesh *mesh) = 0;
    virtual void computeFaceNeighbors(const Mesh *mesh) = 0;
    virtual const std::vector<std::vector<int>>& getVertexToFaceNeighborhood() const = 0;
    virtual const std::vector<std::set<int>>& getVertexToVertexNeighborhood() const = 0;
    virtual const std::vector<std::set<int>>& getFaceToFaceNeighborhood() const = 0;
};

/*
 * Concrete 1-ring neighborhood
 */
class OneRingMeshNeighborhood: public MeshNeighborhood {
public:
    virtual ~OneRingMeshNeighborhood();
    virtual void computeVertexNeighbors(const Mesh *mesh);
    virtual void computeFaceNeighbors(const Mesh *mesh);
    virtual const std::vector<std::vector<int>>& getVertexToFaceNeighborhood() const;
    virtual const std::vector<std::set<int>>& getVertexToVertexNeighborhood() const;
    virtual const std::vector<std::set<int>>& getFaceToFaceNeighborhood() const;

protected:
    std::vector<std::vector<int>> m_vf_1ring; // VERTEX <-> FACE
    std::vector<std::set<int>> m_vv_1ring; // VERTEX <-> VERTEX
    std::vector<std::set<int>> m_ff_1ring; // FACE <-> FACE
};

/*
 * Concrete 2-ring neighborhood
 */
class TwoRingMeshNeighborhood: public OneRingMeshNeighborhood {
public:
    TwoRingMeshNeighborhood(const MeshNeighborhood *smallerNeighborhood);
    virtual ~TwoRingMeshNeighborhood();
    virtual void computeVertexNeighbors(const Mesh *mesh);
    virtual void computeFaceNeighbors(const Mesh *mesh);
    virtual const std::vector<std::vector<int>>& getVertexToFaceNeighborhood() const;
    virtual const std::vector<std::set<int>>& getVertexToVertexNeighborhood() const;
    virtual const std::vector<std::set<int>>& getFaceToFaceNeighborhood() const;

protected:
    std::vector<std::vector<int>> m_vf_2ring; // VERTEX <-> FACE
    std::vector<std::set<int>> m_vv_2ring; // VERTEX <-> VERTEX
    std::vector<std::set<int>> m_ff_2ring; // FACE <-> FACE

private:
    bool m_hasVertexToFace1RingComputed = false;
    bool m_hasVertexToVertex1RingComputed = false;
    bool m_hasFaceToFace1RingComputed = false;
};

#endif //ARBF_TEST_MESHNEIGHBORHOOD_H
