//
//  ARBFInterpolator.cpp
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <cfenv> // for std::fesetround()
#include <tuple>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <Eigen/Eigenvalues>
#include "../include/Config.h"
#include "../include/TriMesh.h"
#include "../include/TetMesh.h"
#include "../include/HexMesh.h"
#include "../include/ARBFInterpolator.h"
#include "../include/MeshNeighborhood.h"

ARBFInterpolator::ARBFInterpolator(): m_basis(nullptr) {}

ARBFInterpolator::~ARBFInterpolator() {
//    delete [] m_centers;
//    m_centers = nullptr;
//    delete [] m_edgeCenters;
//    m_edgeCenters = nullptr;
//    delete [] m_bodyCenters;
//    m_bodyCenters = nullptr;
}

//const std::vector<Vertex>& ARBFInterpolator::getTriangleCenters() const {
//    return m_centers;
//}
//
//const std::vector<Vertex>& ARBFInterpolator::getEdgeCenters() const {
//    return m_edgeCenters;
//}
//
//const std::vector<Vertex>& ARBFInterpolator::getTetrahedronCenters() const {
//    return m_tetrahedronCenters;
//}

ARBFInterpolator::InterpolateResult ARBFInterpolator::getResult() const {
    return m_result;
}

void ARBFInterpolator::setMesh(Mesh * mesh) {
    m_mesh = mesh;
}

void ARBFInterpolator::setBasisFunction(BasisFunction *basisFunction) {
    m_basis = basisFunction;
}

void ARBFInterpolator::computeEdgeMetrics() {
//    m_edgeT.resize(m_mesh->getNumEdges());

    for (auto it = m_mesh->getEdges().begin(); it != m_mesh->getEdges().end(); it++) {
        const Vertex &edgeCenter = it->center;
        //             |lambda1 0      | |v1T|
        // T = [v1 v2] |               | |   |
        //             |0       lambda2| |v2T|
        Eigen::Matrix3d t1, t2;
        t1(0, 0) = edgeCenter.eigVec1[0];
        t1(1, 0) = edgeCenter.eigVec1[1];
        t1(2, 0) = edgeCenter.eigVec1[2];
        t1(0, 1) = edgeCenter.eigVec2[0];
        t1(1, 1) = edgeCenter.eigVec2[1];
        t1(2, 1) = edgeCenter.eigVec2[2];
        t1(0, 2) = edgeCenter.eigVec3[0];
        t1(1, 2) = edgeCenter.eigVec3[1];
        t1(2, 2) = edgeCenter.eigVec3[2];

        t2(0, 0) = edgeCenter.eig1;
        t2(0, 1) = 0.0;
        t2(0, 2) = 0.0;
        t2(1, 0) = 0.0;
        t2(1, 1) = edgeCenter.eig2;
        t2(1, 2) = 0.0;
        t2(2, 0) = 0.0;
        t2(2, 1) = 0.0;
        t2(2, 2) = edgeCenter.eig3;

        // Use definition: tensor_result = t1 * t2 * t1.transpose()
        Eigen::Matrix3d tmp = t1 * t2;
        m_edgeT.push_back(tmp * t1.transpose());
    }

    assert(m_edgeT.size() == m_mesh->getNumEdges());
}

void ARBFInterpolator::computeTetrahedronMetrics() {
    TetMesh *mesh = dynamic_cast<TetMesh*>(m_mesh);

    for (auto &tetrahedron: mesh->getTetrahedrons()) {
        const Vertex &center = tetrahedron.center;
        //             |lambda1 0      | |v1T|
        // T = [v1 v2] |               | |   |
        //             |0       lambda2| |v2T|

        Eigen::Matrix3d t1, t2;
        t1(0, 0) = center.eigVec1[0];
        t1(1, 0) = center.eigVec1[1];
        t1(2, 0) = center.eigVec1[2];
        t1(0, 1) = center.eigVec2[0];
        t1(1, 1) = center.eigVec2[1];
        t1(2, 1) = center.eigVec2[2];
        t1(0, 2) = center.eigVec3[0];
        t1(1, 2) = center.eigVec3[1];
        t1(2, 2) = center.eigVec3[2];

        t2(0, 0) = center.eig1;
        t2(0, 1) = 0.0;
        t2(0, 2) = 0.0;
        t2(1, 0) = 0.0;
        t2(1, 1) = center.eig2;
        t2(1, 2) = 0.0;
        t2(2, 0) = 0.0;
        t2(2, 1) = 0.0;
        t2(2, 2) = center.eig3;

        // Use definition: tensor_result = t1 * t2 * t1.transpose()
        auto temp = t1 * t2;
        m_tetrahedronT.push_back(temp * t1.transpose());
    }
}

inline double ARBFInterpolator::m_computeWeight(double r) {
    return exp(-SQ(r) / SQ(30.0));
}

ARBFInterpolator::InterpolateResult ARBFInterpolator::m_interpolate2dGlobal(unsigned numEvalPoints) {
    // solve coefficients
    unsigned nv = m_mesh->getNumVertices();
    unsigned nf = m_mesh->getNumTriangleFaces();
    unsigned ne = m_mesh->getNumEdges();
    auto faces = m_mesh->getTriangleFacesAsList();
    auto edges = m_mesh->getEdgesAsList();
    unsigned dim = nv + nf + ne;

    m_distanceMatrix.resize(dim, dim);
    m_u.resize(dim);
    const Vertex *vi = nullptr, *vj = nullptr;
    Eigen::Matrix3d tensor;

    for (unsigned i = 0; i < dim; i++) {
        if (i < nv) {
            vi = &m_mesh->getVertices()[i];
        } else if (i < (nv + nf)) {
            vi = &faces[i - nv].center;
        } else {
            vi = &edges[i - nv - nf].center;
        }

        m_u(i) = vi->intensity;
        double p0[] = { vi->x, vi->y, vi->z };

        for (int j = 0; j < dim; j++) {
            double c1 = 0.01;
            if (j < nv) {
                vj = &m_mesh->getVertices()[j];
                tensor = Eigen::Matrix3d::Identity();
            } else if (j < (nv + nf)) {
                vj = &faces[j - nv].center;
                tensor = Eigen::Matrix3d::Identity();
            } else {
                vj = &edges[j - nv - nf].center;
                tensor = m_edgeT[j-nv-nf];
            }

            double p1[] = { vj->x, vj->y, vj->z };
            double r = m_computeDistance(p0, p1, tensor);
            m_distanceMatrix(i, j) = m_basis->phi(r, c1);
        };
    }
    m_coeff = m_distanceMatrix.fullPivLu().solve(m_u); // call Eigen library to solve

    // interpolate
    unsigned dimX = (unsigned) std::sqrt(numEvalPoints);
    unsigned dimY = (unsigned) std::sqrt(numEvalPoints);
    std::vector<double> x = linspace(m_mesh->getMinX(), m_mesh->getMaxX(), dimX);
    std::vector<double> y = linspace(m_mesh->getMinY(), m_mesh->getMaxY(), dimY);
    std::vector<double> intensities(numEvalPoints, -1.0); // solution

    for (unsigned j = 0; j < dimY; j++) {
        for (unsigned i = 0; i < dimX; i++) {
            double p0[] = { x[i], y[j], 0.0 };
            double intensity = 0, phi = 0;

            for (unsigned l = 0; l < dim; l++) {
                double r = 0.0, c1 = 0.01;

                if (l < nv) {
                    double p1[] = { m_mesh->getVertices()[l].x,
                                    m_mesh->getVertices()[l].y,
                                    m_mesh->getVertices()[l].z };
                    r = m_computeDistance(p0, p1, Eigen::Matrix3d::Identity());
                } else if (l < (nv + nf)) {
                    double p1[] = { faces[l - nv].center.x,
                                    faces[l - nv].center.y,
                                    faces[l - nv].center.z };
                    r = m_computeDistance(p0, p1, Eigen::Matrix3d::Identity());
                } else {
                    double p1[] = { edges[l - nv - nf].center.x,
                                    edges[l - nv - nf].center.y,
                                    edges[l - nv - nf].center.z };
                    r = m_computeDistance(p0, p1, m_edgeT[l - nv - nf]);
                }

                intensity += (m_basis->phi(r, c1) * m_coeff(l));
            }

            intensities[j*dimX+i] = intensity;
        }
    }

    return InterpolateResult(intensities, 255.0, { dimX, dimY, 0 });
}

ARBFInterpolator::InterpolateResult ARBFInterpolator::m_interpolateHexahedron3dGlobal(unsigned numEvalPoints) {
    // solve coefficients
    HexMesh *mesh = dynamic_cast<HexMesh*>(m_mesh);
    unsigned nv = mesh->getNumVertices();
    unsigned nf = mesh->getNumQuadrangleFaces();
    unsigned nh = mesh->getNumHexahedrons();
    unsigned rows = nv + nf + nh + (nh * 6);
    unsigned cols = rows;
    m_distanceMatrix.resize(rows, cols);
    m_u.resize(rows);
    auto faces = mesh->getQuadrangleFacesAsList();
    double c1 = 0.1;
    double p0[] = {0, 0, 0}, p1[] = {0, 0, 0};
    double v0[3], v1[3], v2[3], v3[3], v4[3], v5[3], v6[3], w0[3], w1[3], w2[3], w3[3], w4[3], w5[3], w6[3];

    for (unsigned i = 0; i < (nv + nf + nh); i++) {
        m_u(i) = -1;

        if (i < nv) {
            const Vertex &vi = mesh->getVertices()[i];
            p0[0] = vi.x;
            p0[1] = vi.y;
            p0[2] = vi.z;
            m_u(i) = vi.intensity;
        } else if (i < (nv + nf)) {
            const Vertex &vi = faces[i - nv].center;
            p0[0] = vi.x;
            p0[1] = vi.y;
            p0[2] = vi.z;
            m_u(i) = vi.intensity;
        } else {
            const Vertex &vi = mesh->getHexahedrons()[i - nv - nf].center;
            p0[0] = vi.x;
            p0[1] = vi.y;
            p0[2] = vi.z;
            m_u(i) = vi.intensity;
        }

        for (unsigned j = 0; j < cols;) {
            if (j < nv) {
                const Vertex &vj = mesh->getVertices()[j];
                p1[0] = vj.x;
                p1[1] = vj.y;
                p1[2] = vj.z;
                double r = m_computeDistance(p0, p1);
                m_distanceMatrix(i, j) = m_basis->phi(r, c1);
                j++;
            } else if (j < (nv + nf)) {
                const Vertex &vj = faces[j - nv].center;
                p1[0] = vj.x;
                p1[1] = vj.y;
                p1[2] = vj.z;
                double r = m_computeDistance(p0, p1);
                m_distanceMatrix(i, j) = m_basis->phi(r, c1);
                j++;
            } else if (j < (nv + nf + nh)) {
                const Vertex &vj = mesh->getHexahedrons()[j - nv - nf].center;
                p1[0] = vj.x;
                p1[1] = vj.y;
                p1[2] = vj.z;
                double r = m_computeDistance(p0, p1);
                m_distanceMatrix(i, j) = m_basis->phi(r, c1);
                j++;
            } else {
                unsigned n = j - nv - nf - nh;
                const Hexahedron &hj = mesh->getHexahedrons()[n / 6];
                w0[0] = hj.center.x, w0[1] = hj.center.y, w0[2] = hj.center.z;
                w1[0] = hj.f1.center.x, w1[1] = hj.f1.center.y, w1[2] = hj.f1.center.z;
                w2[0] = hj.f2.center.x, w2[1] = hj.f2.center.y, w2[2] = hj.f2.center.z;
                w3[0] = hj.f3.center.x, w3[1] = hj.f3.center.y, w3[2] = hj.f3.center.z;
                w4[0] = hj.f4.center.x, w4[1] = hj.f4.center.y, w4[2] = hj.f4.center.z;
                w5[0] = hj.f5.center.x, w5[1] = hj.f5.center.y, w5[2] = hj.f5.center.z;
                w6[0] = hj.f6.center.x, w6[1] = hj.f6.center.y, w6[2] = hj.f6.center.z;
                double r1 = m_computeDistance(p0, w0, w1);
                double r2 = m_computeDistance(p0, w0, w2);
                double r3 = m_computeDistance(p0, w0, w3);
                double r4 = m_computeDistance(p0, w0, w4);
                double r5 = m_computeDistance(p0, w0, w5);
                double r6 = m_computeDistance(p0, w0, w6);
                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i, j + 1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i, j + 2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i, j + 3) = m_basis->phi(r4, c1);
                m_distanceMatrix(i, j + 4) = m_basis->phi(r5, c1);
                m_distanceMatrix(i, j + 5) = m_basis->phi(r6, c1);
                j += 6;
            }
        }
    }

    for (unsigned i = (nv + nf + nh); i < rows; i += 6) {
        m_u(i) = -1;
        m_u(i + 1) = -1;
        m_u(i + 2) = -1;
        m_u(i + 3) = -1;
        m_u(i + 4) = -1;
        m_u(i + 5) = -1;

        unsigned m = i - nv - nf - nh;
        const Hexahedron &hi = mesh->getHexahedrons()[m / 6];
        v0[0] = hi.center.x, v0[1] = hi.center.y, v0[2] = hi.center.z;
        v1[0] = hi.f1.center.x, v1[1] = hi.f1.center.y, v1[2] = hi.f1.center.z;
        v2[0] = hi.f2.center.x, v2[1] = hi.f2.center.y, v2[2] = hi.f2.center.z;
        v3[0] = hi.f3.center.x, v3[1] = hi.f3.center.y, v3[2] = hi.f3.center.z;
        v4[0] = hi.f4.center.x, v4[1] = hi.f4.center.y, v4[2] = hi.f4.center.z;
        v5[0] = hi.f5.center.x, v5[1] = hi.f5.center.y, v5[2] = hi.f5.center.z;
        v6[0] = hi.f6.center.x, v6[1] = hi.f6.center.y, v6[2] = hi.f6.center.z;

        for (int j = 0; j < cols;) {
            if (j < nv) {
                const Vertex &vj = m_mesh->getVertices()[j];
                p1[0] = vj.x;
                p1[1] = vj.y;
                p1[2] = vj.z;
                double r1 = m_computeDistance(p1, v0, v1);
                double r2 = m_computeDistance(p1, v0, v2);
                double r3 = m_computeDistance(p1, v0, v3);
                double r4 = m_computeDistance(p1, v0, v4);
                double r5 = m_computeDistance(p1, v0, v5);
                double r6 = m_computeDistance(p1, v0, v6);
                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 1, j) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 2, j) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 3, j) = m_basis->phi(r4, c1);
                m_distanceMatrix(i + 4, j) = m_basis->phi(r5, c1);
                m_distanceMatrix(i + 5, j) = m_basis->phi(r6, c1);
                j++;
            } else if (j < (nv + nf)) {
                const Vertex &vj = faces[j - nv].center;
                p1[0] = vj.x;
                p1[1] = vj.y;
                p1[2] = vj.z;
                double r1 = m_computeDistance(p1, v0, v1);
                double r2 = m_computeDistance(p1, v0, v2);
                double r3 = m_computeDistance(p1, v0, v3);
                double r4 = m_computeDistance(p1, v0, v4);
                double r5 = m_computeDistance(p1, v0, v5);
                double r6 = m_computeDistance(p1, v0, v6);
                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 1, j) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 2, j) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 3, j) = m_basis->phi(r4, c1);
                m_distanceMatrix(i + 4, j) = m_basis->phi(r5, c1);
                m_distanceMatrix(i + 5, j) = m_basis->phi(r6, c1);
                j++;
            } else if (j < (nv + nf + nh)) {
                const Vertex &vj = mesh->getHexahedrons()[j - nv - nf].center;
                p1[0] = vj.x;
                p1[1] = vj.y;
                p1[2] = vj.z;
                double r1 = m_computeDistance(p1, v0, v1);
                double r2 = m_computeDistance(p1, v0, v2);
                double r3 = m_computeDistance(p1, v0, v3);
                double r4 = m_computeDistance(p1, v0, v4);
                double r5 = m_computeDistance(p1, v0, v5);
                double r6 = m_computeDistance(p1, v0, v6);
                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 1, j) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 2, j) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 3, j) = m_basis->phi(r4, c1);
                m_distanceMatrix(i + 4, j) = m_basis->phi(r5, c1);
                m_distanceMatrix(i + 5, j) = m_basis->phi(r6, c1);
                j++;
            } else {
                unsigned n = j - nv - nf - nh;
                const Hexahedron &hj = mesh->getHexahedrons()[n / 6];
                double r1, r2, r3, r4, r5, r6;
                w0[0] = hj.center.x, w0[1] = hj.center.y, w0[2] = hj.center.z;
                w1[0] = hj.f1.center.x, w1[1] = hj.f1.center.y, w1[2] = hj.f1.center.z;
                w2[0] = hj.f2.center.x, w2[1] = hj.f2.center.y, w2[2] = hj.f2.center.z;
                w3[0] = hj.f3.center.x, w3[1] = hj.f3.center.y, w3[2] = hj.f3.center.z;
                w4[0] = hj.f4.center.x, w4[1] = hj.f4.center.y, w4[2] = hj.f4.center.z;
                w5[0] = hj.f5.center.x, w5[1] = hj.f5.center.y, w5[2] = hj.f5.center.z;

                r1 = m_computeDistance(v0, v1, w0, w1);
                r2 = m_computeDistance(v0, v1, w0, w2);
                r3 = m_computeDistance(v0, v1, w0, w3);
                r4 = m_computeDistance(v0, v1, w0, w4);
                r5 = m_computeDistance(v0, v1, w0, w5);
                r6 = m_computeDistance(v0, v1, w0, w6);
                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i, j + 1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i, j + 2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i, j + 3) = m_basis->phi(r4, c1);
                m_distanceMatrix(i, j + 4) = m_basis->phi(r5, c1);
                m_distanceMatrix(i, j + 5) = m_basis->phi(r6, c1);

                r1 = m_computeDistance(v0, v2, w0, w1);
                r2 = m_computeDistance(v0, v2, w0, w2);
                r3 = m_computeDistance(v0, v2, w0, w3);
                r4 = m_computeDistance(v0, v2, w0, w4);
                r5 = m_computeDistance(v0, v2, w0, w5);
                r6 = m_computeDistance(v0, v2, w0, w6);
                m_distanceMatrix(i + 1, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 1, j + 1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 1, j + 2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 1, j + 3) = m_basis->phi(r4, c1);
                m_distanceMatrix(i + 1, j + 4) = m_basis->phi(r5, c1);
                m_distanceMatrix(i + 1, j + 5) = m_basis->phi(r6, c1);

                r1 = m_computeDistance(v0, v3, w0, w1);
                r2 = m_computeDistance(v0, v3, w0, w2);
                r3 = m_computeDistance(v0, v3, w0, w3);
                r4 = m_computeDistance(v0, v3, w0, w4);
                r5 = m_computeDistance(v0, v3, w0, w5);
                r6 = m_computeDistance(v0, v3, w0, w6);
                m_distanceMatrix(i + 2, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 2, j + 1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 2, j + 2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 2, j + 3) = m_basis->phi(r4, c1);
                m_distanceMatrix(i + 2, j + 4) = m_basis->phi(r5, c1);
                m_distanceMatrix(i + 2, j + 5) = m_basis->phi(r6, c1);

                r1 = m_computeDistance(v0, v4, w0, w1);
                r2 = m_computeDistance(v0, v4, w0, w2);
                r3 = m_computeDistance(v0, v4, w0, w3);
                r4 = m_computeDistance(v0, v4, w0, w4);
                r5 = m_computeDistance(v0, v4, w0, w5);
                r6 = m_computeDistance(v0, v4, w0, w6);
                m_distanceMatrix(i + 3, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 3, j + 1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 3, j + 2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 3, j + 3) = m_basis->phi(r4, c1);
                m_distanceMatrix(i + 3, j + 4) = m_basis->phi(r5, c1);
                m_distanceMatrix(i + 3, j + 5) = m_basis->phi(r6, c1);

                r1 = m_computeDistance(v0, v5, w0, w1);
                r2 = m_computeDistance(v0, v5, w0, w2);
                r3 = m_computeDistance(v0, v5, w0, w3);
                r4 = m_computeDistance(v0, v5, w0, w4);
                r5 = m_computeDistance(v0, v5, w0, w5);
                r6 = m_computeDistance(v0, v5, w0, w6);
                m_distanceMatrix(i + 4, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 4, j + 1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 4, j + 2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 4, j + 3) = m_basis->phi(r4, c1);
                m_distanceMatrix(i + 4, j + 4) = m_basis->phi(r5, c1);
                m_distanceMatrix(i + 4, j + 5) = m_basis->phi(r6, c1);

                r1 = m_computeDistance(v0, v6, w0, w1);
                r2 = m_computeDistance(v0, v6, w0, w2);
                r3 = m_computeDistance(v0, v6, w0, w3);
                r4 = m_computeDistance(v0, v6, w0, w4);
                r5 = m_computeDistance(v0, v6, w0, w5);
                r6 = m_computeDistance(v0, v6, w0, w6);
                m_distanceMatrix(i + 5, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 5, j + 1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 5, j + 2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 5, j + 3) = m_basis->phi(r4, c1);
                m_distanceMatrix(i + 5, j + 4) = m_basis->phi(r5, c1);
                m_distanceMatrix(i + 5, j + 5) = m_basis->phi(r6, c1);

                j += 6;
            }
        }
    }

    m_coeff = m_distanceMatrix.fullPivLu().solve(m_u); // call Eigen library to solve

    // interpolate
    unsigned dimX = (unsigned) std::round(std::pow(numEvalPoints, 1.0 / 3.0));
    unsigned dimY = (unsigned) std::round(std::pow(numEvalPoints, 1.0 / 3.0));
    unsigned dimZ = (unsigned) std::round(std::pow(numEvalPoints, 1.0 / 3.0));
    std::vector<double> x = linspace(m_mesh->getMinX(), m_mesh->getMaxX(), dimX);
    std::vector<double> y = linspace(m_mesh->getMinY(), m_mesh->getMaxY(), dimY);
    std::vector<double> z = linspace(m_mesh->getMinZ(), m_mesh->getMaxZ(), dimZ);
    std::vector<double> intensities(numEvalPoints, -1.0); // solution

    //    cout << "\nx = ";
    //    copy(x.begin(), x.end(), ostream_iterator<double>(cout, ", "));
    //    cout << "\ny = ";
    //    copy(y.begin(), y.end(), ostream_iterator<double>(cout, ", "));
    //    cout << "\nz = ";
    //    copy(z.begin(), z.end(), ostream_iterator<double>(cout, ", "));

    for (int k = 0; k < dimZ; k++) {
        for (int j = 0; j < dimY; j++) {
            for (int i = 0; i < dimX; i++) {
                p0[0] = x[i], p0[1] = y[j], p0[2] = z[k];
                if (!m_isInAnyHexahedron(p0)) {
                    continue;
                }

                double intensity = 0;
                for (int l = 0; l < cols;) {
                    if (l < nv) {
                        const Vertex &vj = m_mesh->getVertices()[l];
                        p1[0] = vj.x;
                        p1[1] = vj.y;
                        p1[2] = vj.z;
                        double r = m_computeDistance(p0, p1);
                        intensity += (m_basis->phi(r, c1) * m_coeff(l++));
                    } else if (l < (nv + nf)) {
                        const Vertex &vj = faces[l - nv].center;
                        p1[0] = vj.x;
                        p1[1] = vj.y;
                        p1[2] = vj.z;
                        double r = m_computeDistance(p0, p1);
                        intensity += (m_basis->phi(r, c1) * m_coeff(l++));
                    } else if (l < (nv + nf + nh)) {
                        const Vertex &vj = mesh->getHexahedrons()[l - nv - nf].center;
                        p1[0] = vj.x;
                        p1[1] = vj.y;
                        p1[2] = vj.z;
                        double r = m_computeDistance(p0, p1);
                        intensity += (m_basis->phi(r, c1) * m_coeff(l++));
                    } else {
                        unsigned n = l - nv - nf - nh;
                        const Hexahedron &hj = mesh->getHexahedrons()[n / 6];
                        w0[0] = hj.center.x, w0[1] = hj.center.y, w0[2] = hj.center.z;
                        w1[0] = hj.f1.center.x, w1[1] = hj.f1.center.y, w1[2] = hj.f1.center.z;
                        w2[0] = hj.f2.center.x, w2[1] = hj.f2.center.y, w2[2] = hj.f2.center.z;
                        w3[0] = hj.f3.center.x, w3[1] = hj.f3.center.y, w3[2] = hj.f3.center.z;
                        w4[0] = hj.f4.center.x, w4[1] = hj.f4.center.y, w4[2] = hj.f4.center.z;
                        w5[0] = hj.f5.center.x, w5[1] = hj.f5.center.y, w5[2] = hj.f5.center.z;
                        w6[0] = hj.f6.center.x, w6[1] = hj.f6.center.y, w6[2] = hj.f6.center.z;
                        double r1 = m_computeDistance(p0, w0, w1);
                        double r2 = m_computeDistance(p0, w0, w2);
                        double r3 = m_computeDistance(p0, w0, w3);
                        double r4 = m_computeDistance(p0, w0, w4);
                        double r5 = m_computeDistance(p0, w0, w5);
                        double r6 = m_computeDistance(p0, w0, w6);
                        intensity += (m_basis->phi(r1, c1) * m_coeff(l++));
                        intensity += (m_basis->phi(r2, c1) * m_coeff(l++));
                        intensity += (m_basis->phi(r3, c1) * m_coeff(l++));
                        intensity += (m_basis->phi(r4, c1) * m_coeff(l++));
                        intensity += (m_basis->phi(r5, c1) * m_coeff(l++));
                        intensity += (m_basis->phi(r6, c1) * m_coeff(l++));
                    }

//                            if (l < nv) {
//                                p1[0] = m_mesh->getVertices()[l].x;
//                                p1[1] = m_mesh->getVertices()[l].y;
//                                p1[2] = m_mesh->getVertices()[l].z;
//                            } else if (l < (nv + nf)) {
//                                p1[0] = m_centers[l - nv].x;
//                                p1[1] = m_centers[l - nv].y;
//                                p1[2] = m_centers[l - nv].z;
//                            } else {
//                                p1[0] = m_tetrahedronCenters[l - nv - nf].x;
//                                p1[1] = m_tetrahedronCenters[l - nv - nf].y;
//                                p1[2] = m_tetrahedronCenters[l - nv - nf].z;
//                            }

//                            double r = m_computeDistance(p0, p1, tensor);
//                            intensity += (m_basis->phi(r, c1) * m_coeff(l));
                }
                intensities[k * dimX * dimY + j * dimX + i] = intensity;
            }
        }
    }

    return InterpolateResult(intensities, 255.0, { dimX, dimY, dimZ });
}

ARBFInterpolator::InterpolateResult ARBFInterpolator::m_interpolate3dGlobal(unsigned numEvalPoints) {
    // solve coefficients
    TetMesh *mesh = dynamic_cast<TetMesh*>(m_mesh);
    unsigned nv = mesh->getNumVertices();
    unsigned nf = mesh->getNumTriangleFaces();
    unsigned nt = mesh->getNumTetrahedrons();
    unsigned rows = nv + nf + nt + (nt * 4);
    unsigned cols = rows;
    m_distanceMatrix.resize(rows, cols);
    m_u.resize(rows);
    auto faces = mesh->getTriangleFacesAsList();
    double c1 = 0.1;
    double p0[] = {0, 0, 0}, p1[] = {0, 0, 0};
    double v0[3], v1[3], v2[3], v3[3], v4[3], w0[3], w1[3], w2[3], w3[3], w4[3];

    for (unsigned i = 0; i < (nv + nf + nt); i++) {
        m_u(i) = -1;

        if (i < nv) {
            const Vertex &vi = mesh->getVertices()[i];
            p0[0] = vi.x;
            p0[1] = vi.y;
            p0[2] = vi.z;
            m_u(i) = vi.intensity;
        } else if (i < (nv + nf)) {
            const Vertex &vi = faces[i - nv].center;
            p0[0] = vi.x;
            p0[1] = vi.y;
            p0[2] = vi.z;
            m_u(i) = vi.intensity;
        } else {
            const Vertex &vi = mesh->getTetrahedrons()[i - nv - nf].center;
            p0[0] = vi.x;
            p0[1] = vi.y;
            p0[2] = vi.z;
            m_u(i) = vi.intensity;
        }

        for (int j = 0; j < cols;) {
            if (j < nv) {
                const Vertex &vj = m_mesh->getVertices()[j];
                p1[0] = vj.x;
                p1[1] = vj.y;
                p1[2] = vj.z;
                double r = m_computeDistance(p0, p1);
                m_distanceMatrix(i, j) = m_basis->phi(r, c1);
                j++;
            } else if (j < (nv + nf)) {
                const Vertex &vj = faces[j - nv].center;
                p1[0] = vj.x;
                p1[1] = vj.y;
                p1[2] = vj.z;
                double r = m_computeDistance(p0, p1);
                m_distanceMatrix(i, j) = m_basis->phi(r, c1);
                j++;
            } else if (j < (nv + nf + nt)) {
                const Vertex &vj = mesh->getTetrahedrons()[j - nv - nf].center;
                p1[0] = vj.x;
                p1[1] = vj.y;
                p1[2] = vj.z;
                double r = m_computeDistance(p0, p1);
                m_distanceMatrix(i, j) = m_basis->phi(r, c1);
                j++;
            } else {
                unsigned n = j - nv - nf - nt;
                const Tetrahedron &tj = mesh->getTetrahedrons()[n / 4];
                w0[0] = tj.center.x, w0[1] = tj.center.y, w0[2] = tj.center.z;
                w1[0] = tj.f1.center.x, w1[1] = tj.f1.center.y, w1[2] = tj.f1.center.z;
                w2[0] = tj.f2.center.x, w2[1] = tj.f2.center.y, w2[2] = tj.f2.center.z;
                w3[0] = tj.f3.center.x, w3[1] = tj.f3.center.y, w3[2] = tj.f3.center.z;
                w4[0] = tj.f4.center.x, w4[1] = tj.f4.center.y, w4[2] = tj.f4.center.z;
                double r1 = m_computeDistance(p0, w0, w1);
                double r2 = m_computeDistance(p0, w0, w2);
                double r3 = m_computeDistance(p0, w0, w3);
                double r4 = m_computeDistance(p0, w0, w4);
                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i, j + 1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i, j + 2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i, j + 3) = m_basis->phi(r4, c1);
                j += 4;
            }
        }
    }

    for (unsigned i = (nv + nf + nt); i < rows; i += 4) {
        m_u(i) = -1;
        m_u(i + 1) = -1;
        m_u(i + 2) = -1;
        m_u(i + 3) = -1;

        unsigned m = i - nv - nf - nt;
        const Tetrahedron &ti = mesh->getTetrahedrons()[m / 4];
        v0[0] = ti.center.x, v0[1] = ti.center.y, v0[2] = ti.center.z;
        v1[0] = ti.f1.center.x, v1[1] = ti.f1.center.y, v1[2] = ti.f1.center.z;
        v2[0] = ti.f2.center.x, v2[1] = ti.f2.center.y, v2[2] = ti.f2.center.z;
        v3[0] = ti.f3.center.x, v3[1] = ti.f3.center.y, v3[2] = ti.f3.center.z;
        v4[0] = ti.f4.center.x, v4[1] = ti.f4.center.y, v4[2] = ti.f4.center.z;

        for (int j = 0; j < cols;) {
            if (j < nv) {
                const Vertex &vj = m_mesh->getVertices()[j];
                p1[0] = vj.x;
                p1[1] = vj.y;
                p1[2] = vj.z;
                double r1 = m_computeDistance(p1, v0, v1);
                double r2 = m_computeDistance(p1, v0, v2);
                double r3 = m_computeDistance(p1, v0, v3);
                double r4 = m_computeDistance(p1, v0, v4);
                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 1, j) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 2, j) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 3, j) = m_basis->phi(r4, c1);
                j++;
            } else if (j < (nv + nf)) {
                const Vertex &vj = faces[j - nv].center;
                p1[0] = vj.x;
                p1[1] = vj.y;
                p1[2] = vj.z;
                double r1 = m_computeDistance(p1, v0, v1);
                double r2 = m_computeDistance(p1, v0, v2);
                double r3 = m_computeDistance(p1, v0, v3);
                double r4 = m_computeDistance(p1, v0, v4);
                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 1, j) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 2, j) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 3, j) = m_basis->phi(r4, c1);
                j++;
            } else if (j < (nv + nf + nt)) {
                const Vertex &vj = mesh->getTetrahedrons()[j - nv - nf].center;
                p1[0] = vj.x;
                p1[1] = vj.y;
                p1[2] = vj.z;
                double r1 = m_computeDistance(p1, v0, v1);
                double r2 = m_computeDistance(p1, v0, v2);
                double r3 = m_computeDistance(p1, v0, v3);
                double r4 = m_computeDistance(p1, v0, v4);
                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 1, j) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 2, j) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 3, j) = m_basis->phi(r4, c1);
                j++;
            } else {
                unsigned n = j - nv - nf - nt;
                const Tetrahedron &tj = mesh->getTetrahedrons()[n / 4];
                double r1, r2, r3, r4;
                w0[0] = tj.center.x, w0[1] = tj.center.y, w0[2] = tj.center.z;
                w1[0] = tj.f1.center.x, w1[1] = tj.f1.center.y, w1[2] = tj.f1.center.z;
                w2[0] = tj.f2.center.x, w2[1] = tj.f2.center.y, w2[2] = tj.f2.center.z;
                w3[0] = tj.f3.center.x, w3[1] = tj.f3.center.y, w3[2] = tj.f3.center.z;
                w4[0] = tj.f4.center.x, w4[1] = tj.f4.center.y, w4[2] = tj.f4.center.z;

                r1 = m_computeDistance(v0, v1, w0, w1);
                r2 = m_computeDistance(v0, v1, w0, w2);
                r3 = m_computeDistance(v0, v1, w0, w3);
                r4 = m_computeDistance(v0, v1, w0, w4);
                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i, j + 1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i, j + 2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i, j + 3) = m_basis->phi(r4, c1);

                r1 = m_computeDistance(v0, v2, w0, w1);
                r2 = m_computeDistance(v0, v2, w0, w2);
                r3 = m_computeDistance(v0, v2, w0, w3);
                r4 = m_computeDistance(v0, v2, w0, w4);
                m_distanceMatrix(i + 1, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 1, j + 1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 1, j + 2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 1, j + 3) = m_basis->phi(r4, c1);

                r1 = m_computeDistance(v0, v3, w0, w1);
                r2 = m_computeDistance(v0, v3, w0, w2);
                r3 = m_computeDistance(v0, v3, w0, w3);
                r4 = m_computeDistance(v0, v3, w0, w4);
                m_distanceMatrix(i + 2, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 2, j + 1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 2, j + 2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 2, j + 3) = m_basis->phi(r4, c1);

                r1 = m_computeDistance(v0, v4, w0, w1);
                r2 = m_computeDistance(v0, v4, w0, w2);
                r3 = m_computeDistance(v0, v4, w0, w3);
                r4 = m_computeDistance(v0, v4, w0, w4);
                m_distanceMatrix(i + 3, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i + 3, j + 1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i + 3, j + 2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i + 3, j + 3) = m_basis->phi(r4, c1);

                j += 4;
            }
        }
    }

    m_coeff = m_distanceMatrix.fullPivLu().solve(m_u); // call Eigen library to solve

    // interpolate
    unsigned dimX = (unsigned) std::round(std::pow(numEvalPoints, 1.0 / 3.0));
    unsigned dimY = (unsigned) std::round(std::pow(numEvalPoints, 1.0 / 3.0));
    unsigned dimZ = (unsigned) std::round(std::pow(numEvalPoints, 1.0 / 3.0));
    std::vector<double> x = linspace(m_mesh->getMinX(), m_mesh->getMaxX(), dimX);
    std::vector<double> y = linspace(m_mesh->getMinY(), m_mesh->getMaxY(), dimY);
    std::vector<double> z = linspace(m_mesh->getMinZ(), m_mesh->getMaxZ(), dimZ);
    std::vector<double> intensities(numEvalPoints, -1.0); // solution

    //    cout << "\nx = ";
    //    copy(x.begin(), x.end(), ostream_iterator<double>(cout, ", "));
    //    cout << "\ny = ";
    //    copy(y.begin(), y.end(), ostream_iterator<double>(cout, ", "));
    //    cout << "\nz = ";
    //    copy(z.begin(), z.end(), ostream_iterator<double>(cout, ", "));

    for (int k = 0; k < dimZ; k++) {
        for (int j = 0; j < dimY; j++) {
            for (int i = 0; i < dimX; i++) {
                p0[0] = x[i], p0[1] = y[j], p0[2] = z[k];
                if (!m_isInAnyTetrahedron(p0)) {
                    continue;
                }

                double intensity = 0;
                for (int l = 0; l < cols;) {
                    if (l < nv) {
                        const Vertex &vj = m_mesh->getVertices()[l];
                        p1[0] = vj.x;
                        p1[1] = vj.y;
                        p1[2] = vj.z;
                        double r = m_computeDistance(p0, p1);
                        intensity += (m_basis->phi(r, c1) * m_coeff(l++));
                    } else if (l < (nv + nf)) {
                        const Vertex &vj = faces[l - nv].center;
                        p1[0] = vj.x;
                        p1[1] = vj.y;
                        p1[2] = vj.z;
                        double r = m_computeDistance(p0, p1);
                        intensity += (m_basis->phi(r, c1) * m_coeff(l++));
                    } else if (l < (nv + nf + nt)) {
                        const Vertex &vj = mesh->getTetrahedrons()[l - nv - nf].center;
                        p1[0] = vj.x;
                        p1[1] = vj.y;
                        p1[2] = vj.z;
                        double r = m_computeDistance(p0, p1);
                        intensity += (m_basis->phi(r, c1) * m_coeff(l++));
                    } else {
                        unsigned n = l - nv - nf - nt;
                        const Tetrahedron &tj = mesh->getTetrahedrons()[n / 4];
                        w0[0] = tj.center.x, w0[1] = tj.center.y, w0[2] = tj.center.z;
                        w1[0] = tj.f1.center.x, w1[1] = tj.f1.center.y, w1[2] = tj.f1.center.z;
                        w2[0] = tj.f2.center.x, w2[1] = tj.f2.center.y, w2[2] = tj.f2.center.z;
                        w3[0] = tj.f3.center.x, w3[1] = tj.f3.center.y, w3[2] = tj.f3.center.z;
                        w4[0] = tj.f4.center.x, w4[1] = tj.f4.center.y, w4[2] = tj.f4.center.z;
                        double r1 = m_computeDistance(p0, w0, w1);
                        double r2 = m_computeDistance(p0, w0, w2);
                        double r3 = m_computeDistance(p0, w0, w3);
                        double r4 = m_computeDistance(p0, w0, w4);
                        intensity += (m_basis->phi(r1, c1) * m_coeff(l++));
                        intensity += (m_basis->phi(r2, c1) * m_coeff(l++));
                        intensity += (m_basis->phi(r3, c1) * m_coeff(l++));
                        intensity += (m_basis->phi(r4, c1) * m_coeff(l++));
                    }

//                            if (l < nv) {
//                                p1[0] = m_mesh->getVertices()[l].x;
//                                p1[1] = m_mesh->getVertices()[l].y;
//                                p1[2] = m_mesh->getVertices()[l].z;
//                            } else if (l < (nv + nf)) {
//                                p1[0] = m_centers[l - nv].x;
//                                p1[1] = m_centers[l - nv].y;
//                                p1[2] = m_centers[l - nv].z;
//                            } else {
//                                p1[0] = m_tetrahedronCenters[l - nv - nf].x;
//                                p1[1] = m_tetrahedronCenters[l - nv - nf].y;
//                                p1[2] = m_tetrahedronCenters[l - nv - nf].z;
//                            }

//                            double r = m_computeDistance(p0, p1, tensor);
//                            intensity += (m_basis->phi(r, c1) * m_coeff(l));
                }
                intensities[k * dimX * dimY + j * dimX + i] = intensity;
            }
        }
    }

    return InterpolateResult(intensities, 255.0, { dimX, dimY, dimZ });
}

void ARBFInterpolator::interpolate_global(unsigned numEvalPoints) {
    switch (Config::problemDim) {
        case 2:
            m_result = m_interpolate2dGlobal(numEvalPoints);
            break;
        case 3:
        default:
            switch (Config::meshType) {
                case Config::MeshType::Tetrahedron:
                default:
                    m_result = m_interpolate3dGlobal(numEvalPoints);
                    break;
                case Config::MeshType::Hexahedron:
                    m_result = m_interpolateHexahedron3dGlobal(numEvalPoints);
                    break;
            }
            break;
    }

    rescaleResultData(0, 255);
    thresholdResultData(std::get<1>(m_result));
}

//ARBFInterpolator::InterpolateResult ARBFInterpolator::m_interpolate3d_local(unsigned numEvalPoints) {
//    TetMesh *mesh = dynamic_cast<TetMesh*>(m_mesh);
//    double p0[] = {0, 0, 0}, p1[] = {0, 0, 0};
//    unsigned dimX = (unsigned) std::round(std::pow(numEvalPoints, 1.0 / 3.0));
//    unsigned dimY = (unsigned) std::round(std::pow(numEvalPoints, 1.0 / 3.0));
//    unsigned dimZ = (unsigned) std::round(std::pow(numEvalPoints, 1.0 / 3.0));
//    std::vector<double> x = linspace(m_mesh->getMinX(), m_mesh->getMaxX(), dimX);
//    std::vector<double> y = linspace(m_mesh->getMinY(), m_mesh->getMaxY(), dimY);
//    std::vector<double> z = linspace(m_mesh->getMinZ(), m_mesh->getMaxZ(), dimZ);
//    std::vector<double> intensities(numEvalPoints, -1.0); // solution
//
//    for (const Tetrahedron &tet: mesh->getTetrahedrons()) {
//        OneRingMeshNeighborhood neighborhood
//    }
//
//
//
//
//    for (int k = 0; k < dimZ; k++) {
//        for (int j = 0; j < dimY; j++) {
//            for (int i = 0; i < dimX; i++) {
//                p0[0] = x[i], p0[1] = y[j], p0[2] = z[k];
//                for (int t = 0; t < )
//
//
//
//                if (!m_isInAnyTetrahedron(p0)) {
//                    continue;
//                }
//
//                double intensity = 0;
//                for (int l = 0; l < cols;) {
//                    if (l < nv) {
//                        const Vertex &vj = m_mesh->getVertices()[l];
//                        p1[0] = vj.x;
//                        p1[1] = vj.y;
//                        p1[2] = vj.z;
//                        double r = m_computeDistance(p0, p1);
//                        intensity += (m_basis->phi(r, c1) * m_coeff(l++));
//                    } else if (l < (nv + nf)) {
//                        const Vertex &vj = faces[l - nv].center;
//                        p1[0] = vj.x;
//                        p1[1] = vj.y;
//                        p1[2] = vj.z;
//                        double r = m_computeDistance(p0, p1);
//                        intensity += (m_basis->phi(r, c1) * m_coeff(l++));
//                    } else if (l < (nv + nf + nt)) {
//                        const Vertex &vj = mesh->getTetrahedrons()[l - nv - nf].center;
//                        p1[0] = vj.x;
//                        p1[1] = vj.y;
//                        p1[2] = vj.z;
//                        double r = m_computeDistance(p0, p1);
//                        intensity += (m_basis->phi(r, c1) * m_coeff(l++));
//                    } else {
//                        unsigned n = l - nv - nf - nt;
//                        const Tetrahedron &tj = mesh->getTetrahedrons()[n / 4];
//                        w0[0] = tj.center.x, w0[1] = tj.center.y, w0[2] = tj.center.z;
//                        w1[0] = tj.f1.center.x, w1[1] = tj.f1.center.y, w1[2] = tj.f1.center.z;
//                        w2[0] = tj.f2.center.x, w2[1] = tj.f2.center.y, w2[2] = tj.f2.center.z;
//                        w3[0] = tj.f3.center.x, w3[1] = tj.f3.center.y, w3[2] = tj.f3.center.z;
//                        w4[0] = tj.f4.center.x, w4[1] = tj.f4.center.y, w4[2] = tj.f4.center.z;
//                        double r1 = m_computeDistance(p0, w0, w1);
//                        double r2 = m_computeDistance(p0, w0, w2);
//                        double r3 = m_computeDistance(p0, w0, w3);
//                        double r4 = m_computeDistance(p0, w0, w4);
//                        intensity += (m_basis->phi(r1, c1) * m_coeff(l++));
//                        intensity += (m_basis->phi(r2, c1) * m_coeff(l++));
//                        intensity += (m_basis->phi(r3, c1) * m_coeff(l++));
//                        intensity += (m_basis->phi(r4, c1) * m_coeff(l++));
//                    }
//
////                            if (l < nv) {
////                                p1[0] = m_mesh->getVertices()[l].x;
////                                p1[1] = m_mesh->getVertices()[l].y;
////                                p1[2] = m_mesh->getVertices()[l].z;
////                            } else if (l < (nv + nf)) {
////                                p1[0] = m_centers[l - nv].x;
////                                p1[1] = m_centers[l - nv].y;
////                                p1[2] = m_centers[l - nv].z;
////                            } else {
////                                p1[0] = m_tetrahedronCenters[l - nv - nf].x;
////                                p1[1] = m_tetrahedronCenters[l - nv - nf].y;
////                                p1[2] = m_tetrahedronCenters[l - nv - nf].z;
////                            }
//
////                            double r = m_computeDistance(p0, p1, tensor);
////                            intensity += (m_basis->phi(r, c1) * m_coeff(l));
//                }
//                intensities[k * dimX * dimY + j * dimX + i] = intensity;
//            }
//        }
//    }
//
//
//
//
//
//
//    // solve coefficients
//
//    unsigned nv = mesh->getNumVertices();
//    unsigned nf = mesh->getNumFaces();
//    unsigned nt = mesh->getNumTetrahedrons();
//    unsigned rows = nv + nf + nt + (nt * 4);
//    unsigned cols = rows;
//    m_distanceMatrix.resize(rows, cols);
//    m_u.resize(rows);
//    auto faces = mesh->getTriangleFacesAsList();
//    double c1 = 0.1;
//
//    double v0[3], v1[3], v2[3], v3[3], v4[3], w0[3], w1[3], w2[3], w3[3], w4[3];
//
//    for (unsigned i = 0; i < (nv + nf + nt); i++) {
//        m_u(i) = -1;
//
//        if (i < nv) {
//            const Vertex &vi = mesh->getVertices()[i];
//            p0[0] = vi.x;
//            p0[1] = vi.y;
//            p0[2] = vi.z;
//            m_u(i) = vi.intensity;
//        } else if (i < (nv + nf)) {
//            const Vertex &vi = faces[i - nv].center;
//            p0[0] = vi.x;
//            p0[1] = vi.y;
//            p0[2] = vi.z;
//            m_u(i) = vi.intensity;
//        } else {
//            const Vertex &vi = mesh->getTetrahedrons()[i - nv - nf].center;
//            p0[0] = vi.x;
//            p0[1] = vi.y;
//            p0[2] = vi.z;
//            m_u(i) = vi.intensity;
//        }
//
//        for (int j = 0; j < cols;) {
//            if (j < nv) {
//                const Vertex &vj = m_mesh->getVertices()[j];
//                p1[0] = vj.x;
//                p1[1] = vj.y;
//                p1[2] = vj.z;
//                double r = m_computeDistance(p0, p1);
//                m_distanceMatrix(i, j) = m_basis->phi(r, c1);
//                j++;
//            } else if (j < (nv + nf)) {
//                const Vertex &vj = faces[j - nv].center;
//                p1[0] = vj.x;
//                p1[1] = vj.y;
//                p1[2] = vj.z;
//                double r = m_computeDistance(p0, p1);
//                m_distanceMatrix(i, j) = m_basis->phi(r, c1);
//                j++;
//            } else if (j < (nv + nf + nt)) {
//                const Vertex &vj = mesh->getTetrahedrons()[j - nv - nf].center;
//                p1[0] = vj.x;
//                p1[1] = vj.y;
//                p1[2] = vj.z;
//                double r = m_computeDistance(p0, p1);
//                m_distanceMatrix(i, j) = m_basis->phi(r, c1);
//                j++;
//            } else {
//                unsigned n = j - nv - nf - nt;
//                const Tetrahedron &tj = mesh->getTetrahedrons()[n / 4];
//                w0[0] = tj.center.x, w0[1] = tj.center.y, w0[2] = tj.center.z;
//                w1[0] = tj.f1.center.x, w1[1] = tj.f1.center.y, w1[2] = tj.f1.center.z;
//                w2[0] = tj.f2.center.x, w2[1] = tj.f2.center.y, w2[2] = tj.f2.center.z;
//                w3[0] = tj.f3.center.x, w3[1] = tj.f3.center.y, w3[2] = tj.f3.center.z;
//                w4[0] = tj.f4.center.x, w4[1] = tj.f4.center.y, w4[2] = tj.f4.center.z;
//                double r1 = m_computeDistance(p0, w0, w1);
//                double r2 = m_computeDistance(p0, w0, w2);
//                double r3 = m_computeDistance(p0, w0, w3);
//                double r4 = m_computeDistance(p0, w0, w4);
//                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
//                m_distanceMatrix(i, j + 1) = m_basis->phi(r2, c1);
//                m_distanceMatrix(i, j + 2) = m_basis->phi(r3, c1);
//                m_distanceMatrix(i, j + 3) = m_basis->phi(r4, c1);
//                j += 4;
//            }
//        }
//    }
//
//    for (unsigned i = (nv + nf + nt); i < rows; i += 4) {
//        m_u(i) = -1;
//        m_u(i + 1) = -1;
//        m_u(i + 2) = -1;
//        m_u(i + 3) = -1;
//
//        unsigned m = i - nv - nf - nt;
//        const Tetrahedron &ti = mesh->getTetrahedrons()[m / 4];
//        v0[0] = ti.center.x, v0[1] = ti.center.y, v0[2] = ti.center.z;
//        v1[0] = ti.f1.center.x, v1[1] = ti.f1.center.y, v1[2] = ti.f1.center.z;
//        v2[0] = ti.f2.center.x, v2[1] = ti.f2.center.y, v2[2] = ti.f2.center.z;
//        v3[0] = ti.f3.center.x, v3[1] = ti.f3.center.y, v3[2] = ti.f3.center.z;
//        v4[0] = ti.f4.center.x, v4[1] = ti.f4.center.y, v4[2] = ti.f4.center.z;
//
//        for (int j = 0; j < cols;) {
//            if (j < nv) {
//                const Vertex &vj = m_mesh->getVertices()[j];
//                p1[0] = vj.x;
//                p1[1] = vj.y;
//                p1[2] = vj.z;
//                double r1 = m_computeDistance(p1, v0, v1);
//                double r2 = m_computeDistance(p1, v0, v2);
//                double r3 = m_computeDistance(p1, v0, v3);
//                double r4 = m_computeDistance(p1, v0, v4);
//                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
//                m_distanceMatrix(i + 1, j) = m_basis->phi(r2, c1);
//                m_distanceMatrix(i + 2, j) = m_basis->phi(r3, c1);
//                m_distanceMatrix(i + 3, j) = m_basis->phi(r4, c1);
//                j++;
//            } else if (j < (nv + nf)) {
//                const Vertex &vj = faces[j - nv].center;
//                p1[0] = vj.x;
//                p1[1] = vj.y;
//                p1[2] = vj.z;
//                double r1 = m_computeDistance(p1, v0, v1);
//                double r2 = m_computeDistance(p1, v0, v2);
//                double r3 = m_computeDistance(p1, v0, v3);
//                double r4 = m_computeDistance(p1, v0, v4);
//                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
//                m_distanceMatrix(i + 1, j) = m_basis->phi(r2, c1);
//                m_distanceMatrix(i + 2, j) = m_basis->phi(r3, c1);
//                m_distanceMatrix(i + 3, j) = m_basis->phi(r4, c1);
//                j++;
//            } else if (j < (nv + nf + nt)) {
//                const Vertex &vj = mesh->getTetrahedrons()[j - nv - nf].center;
//                p1[0] = vj.x;
//                p1[1] = vj.y;
//                p1[2] = vj.z;
//                double r1 = m_computeDistance(p1, v0, v1);
//                double r2 = m_computeDistance(p1, v0, v2);
//                double r3 = m_computeDistance(p1, v0, v3);
//                double r4 = m_computeDistance(p1, v0, v4);
//                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
//                m_distanceMatrix(i + 1, j) = m_basis->phi(r2, c1);
//                m_distanceMatrix(i + 2, j) = m_basis->phi(r3, c1);
//                m_distanceMatrix(i + 3, j) = m_basis->phi(r4, c1);
//                j++;
//            } else {
//                unsigned n = j - nv - nf - nt;
//                const Tetrahedron &tj = mesh->getTetrahedrons()[n / 4];
//                double r1, r2, r3, r4;
//                w0[0] = tj.center.x, w0[1] = tj.center.y, w0[2] = tj.center.z;
//                w1[0] = tj.f1.center.x, w1[1] = tj.f1.center.y, w1[2] = tj.f1.center.z;
//                w2[0] = tj.f2.center.x, w2[1] = tj.f2.center.y, w2[2] = tj.f2.center.z;
//                w3[0] = tj.f3.center.x, w3[1] = tj.f3.center.y, w3[2] = tj.f3.center.z;
//                w4[0] = tj.f4.center.x, w4[1] = tj.f4.center.y, w4[2] = tj.f4.center.z;
//
//                r1 = m_computeDistance(v0, v1, w0, w1);
//                r2 = m_computeDistance(v0, v1, w0, w2);
//                r3 = m_computeDistance(v0, v1, w0, w3);
//                r4 = m_computeDistance(v0, v1, w0, w4);
//                m_distanceMatrix(i, j) = m_basis->phi(r1, c1);
//                m_distanceMatrix(i, j + 1) = m_basis->phi(r2, c1);
//                m_distanceMatrix(i, j + 2) = m_basis->phi(r3, c1);
//                m_distanceMatrix(i, j + 3) = m_basis->phi(r4, c1);
//
//                r1 = m_computeDistance(v0, v2, w0, w1);
//                r2 = m_computeDistance(v0, v2, w0, w2);
//                r3 = m_computeDistance(v0, v2, w0, w3);
//                r4 = m_computeDistance(v0, v2, w0, w4);
//                m_distanceMatrix(i + 1, j) = m_basis->phi(r1, c1);
//                m_distanceMatrix(i + 1, j + 1) = m_basis->phi(r2, c1);
//                m_distanceMatrix(i + 1, j + 2) = m_basis->phi(r3, c1);
//                m_distanceMatrix(i + 1, j + 3) = m_basis->phi(r4, c1);
//
//                r1 = m_computeDistance(v0, v3, w0, w1);
//                r2 = m_computeDistance(v0, v3, w0, w2);
//                r3 = m_computeDistance(v0, v3, w0, w3);
//                r4 = m_computeDistance(v0, v3, w0, w4);
//                m_distanceMatrix(i + 2, j) = m_basis->phi(r1, c1);
//                m_distanceMatrix(i + 2, j + 1) = m_basis->phi(r2, c1);
//                m_distanceMatrix(i + 2, j + 2) = m_basis->phi(r3, c1);
//                m_distanceMatrix(i + 2, j + 3) = m_basis->phi(r4, c1);
//
//                r1 = m_computeDistance(v0, v4, w0, w1);
//                r2 = m_computeDistance(v0, v4, w0, w2);
//                r3 = m_computeDistance(v0, v4, w0, w3);
//                r4 = m_computeDistance(v0, v4, w0, w4);
//                m_distanceMatrix(i + 3, j) = m_basis->phi(r1, c1);
//                m_distanceMatrix(i + 3, j + 1) = m_basis->phi(r2, c1);
//                m_distanceMatrix(i + 3, j + 2) = m_basis->phi(r3, c1);
//                m_distanceMatrix(i + 3, j + 3) = m_basis->phi(r4, c1);
//
//                j += 4;
//            }
//        }
//    }
//
//    m_coeff = m_distanceMatrix.fullPivLu().solve(m_u); // call Eigen library to solve
//
//    // interpolate
//
//
//    //    cout << "\nx = ";
//    //    copy(x.begin(), x.end(), ostream_iterator<double>(cout, ", "));
//    //    cout << "\ny = ";
//    //    copy(y.begin(), y.end(), ostream_iterator<double>(cout, ", "));
//    //    cout << "\nz = ";
//    //    copy(z.begin(), z.end(), ostream_iterator<double>(cout, ", "));
//
//
//
//    return InterpolateResult(intensities, 255.0, { dimX, dimY, dimZ });
//}
ARBFInterpolator::InterpolateResult ARBFInterpolator::m_interpolate2d_local(unsigned numEvalPoints) {
    return InterpolateResult();
}

ARBFInterpolator::InterpolateResult ARBFInterpolator::m_interpolate3d_local(unsigned numEvalPoints) {
    return InterpolateResult();
}

void ARBFInterpolator::interpolate_local(unsigned numEvalPoints) {
    switch (Config::problemDim) {
        case 2:
            m_result = m_interpolate2d_local(numEvalPoints);
            break;
        case 3:
        default:
            m_result = m_interpolate3d_local(numEvalPoints);
            break;
    }

    rescaleResultData(0, 255);
    thresholdResultData(std::get<1>(m_result));
}

double ARBFInterpolator::computePSNR(const PPMImage &img1, const PPMImage &img2) {
    return 20.0 * log10(255.0 / m_computeRMSE(img1, img2));
}

inline double ARBFInterpolator::m_computeRMSE(const PPMImage &img1, const PPMImage &img2) {
    double sum = 0.0;
    assert(img1.getX() == img2.getX() && img1.getY() == img2.getY());
    
    for (int j = 0; j < img1.getY(); ++j) {
        for (int i = 0; i < img1.getX(); ++i) {
            sum += SQ(img1.getImageData()[j * img1.getX() + i] - img2.getImageData()[j * img2.getX() + i]);
            // sum += SQ(img1.getImageData[j*N1+i] - piecewise[j*N1+i]);
        }
    }
    
    return sqrt(sum / (img1.getX() * img1.getY()));
}

//std::tuple<double, double, std::vector<double>, std::vector<double>> ARBFInterpolator::m_computeEigen(const Eigen::Matrix2d &m) {
ARBFInterpolator::Eigens2d ARBFInterpolator::computeEigens(const Eigen::Matrix2d &matrix) {
    assert(std::abs(matrix(0,1) - matrix(1,0)) < Config::epsilon); // matrix should be symmetric

    // compute eigenvalues
    // matrix looks like
    //      | a    c|
    // A =  |       |
    //      | c    b|

    double a = matrix(0,0), b = matrix(1,1), c = matrix(0,1);

    // compute eigenvalues
    double lambda1 = 0.5 * (a + b) + std::sqrt(SQ(c) + 0.25 * SQ(a - b));
    double lambda2 = 0.5 * (a + b) - std::sqrt(SQ(c) + 0.25 * SQ(a - b));
    Eigenvalues2d eigenvalues = { lambda1, lambda2 };

    // compute eigenvectors
    double v1[2], v2[2];
    double u1 = 0.0, w1 = 0.0;

    if (std::abs(lambda1 - a) > Config::epsilon) {
        // lambda1 != a
        u1 = c;
        w1 = lambda1 - a;
    } else if (std::abs(lambda1 - b) > Config::epsilon) {
        // lambda1 != b
        u1 = lambda1 - b;
        w1 = c;
    } else {
        v1[0] = 1.0;
        v1[1] = 0.0;
        v2[0] = v1[1];
        v2[1] = -v1[0];
    }

    double p = std::sqrt(SQ(u1) + SQ(w1));
    v1[0] = u1 / p;
    v1[1] = w1 / p;
    v2[0] = v1[1];
    v2[1] = -v1[0];
    std::array<double, 2> vec1 = { v1[0], v1[1] };
    std::array<double, 2> vec2 = { v2[0], v2[1] };
    Eigenvectors2d eigenvectors = { vec1, vec2 };

    return std::make_tuple(eigenvalues, eigenvectors);
}

ARBFInterpolator::Eigens3d ARBFInterpolator::computeEigens(const Eigen::Matrix3d &matrix) {
    Eigen::EigenSolver<Eigen::Matrix3d> solver(matrix);
    auto values = solver.eigenvalues();

    // make sure all eigenvalues are real
    assert(values[0].imag() < std::abs(Config::epsilon) &&
           values[1].imag() < std::abs(Config::epsilon) &&
           values[2].imag() < std::abs(Config::epsilon));

    Eigenvalues3d eigenvalues = { values(0).real(), values(1).real(), values(2).real() };
//    auto minmax = std::minmax_element(eigenvalues.begin(), eigenvalues.end());
//    auto minPos = minmax.first - eigenvalues.begin();
//    auto maxPos = minmax.second - eigenvalues.begin();
//    unsigned middlePos = 0;
//
//    for (unsigned i = 0; i < eigenvalues.size(); i++) {
//        if (eigenvalues[i] != *minmax.first && eigenvalues[i] != *minmax.second) {
//            middlePos = i;
//            break;
//        }
//    }

    // put the minimum eigenvalues first
//    Eigenvalues3d sortedEigenvalues = { eigenvalues[minPos], eigenvalues[middlePos], eigenvalues[maxPos] };
    auto vectors = solver.eigenvectors();
//    std::array<double, 3> vec1 = { vectors.col(minPos)[0].real(), vectors.col(minPos)[1].real(), vectors.col(minPos)[2].real() };
//    std::array<double, 3> vec2 = { vectors.col(middlePos)[0].real(), vectors.col(middlePos)[1].real(), vectors.col(middlePos)[2].real() };
//    std::array<double, 3> vec3 = { vectors.col(maxPos)[0].real(), vectors.col(maxPos)[1].real(), vectors.col(maxPos)[2].real() };
    std::array<double, 3> vec1 = {vectors.col(0)[0].real(), vectors.col(0)[1].real(), vectors.col(0)[2].real()};
    std::array<double, 3> vec2 = {vectors.col(1)[0].real(), vectors.col(1)[1].real(), vectors.col(1)[2].real()};
    std::array<double, 3> vec3 = {vectors.col(2)[0].real(), vectors.col(2)[1].real(), vectors.col(2)[2].real()};
    Eigenvectors3d eigenvectors = { vec1, vec2, vec3 };
//    return std::make_tuple(sortedEigenvalues, eigenvectors);
    return std::make_tuple(eigenvalues, eigenvectors);
}

void ARBFInterpolator::rescaleEigenvalues(Eigenvalues2d &values) {
    double mu1 = values[0], mu2 = values[1];
    double sumMu = mu1 + mu2;
    double c = 1.0;

    if (sumMu > Config::epsilon || sumMu < -Config::epsilon) {
        c = std::abs(mu1 - mu2) / sumMu;
    }

    if (c > Config::epsilon || c < -Config::epsilon) {
        values[1] = 1.0 - std::exp(-3.8 * SQ(Config::c0 / c));
    } else {
        values[1] = 1.0;
    }

    values[0] = 1.0;
}

void ARBFInterpolator::rescaleEigenvalues(Eigenvalues3d &values) {
    double mu1 = values[0], mu2 = values[1], mu3 = values[2];
    double sum1 = mu1 + mu2;
    double sum2 = mu1 + mu3;
    double c1 = 1.0, c2 = 1.0;

    if (sum1 > Config::epsilon || sum1 < -Config::epsilon) {
        c1 = std::abs(mu1 - mu2) / sum1;
    }

    if (sum2 > Config::epsilon || sum2 < -Config::epsilon) {
        c2 = std::abs(mu1 - mu3) / sum2;
    }

    if (c1 > Config::epsilon || c1 < -Config::epsilon) {
        values[1] = 1.0 - exp(-3.8 * SQ(Config::c0 / c1));
        values[2] = 1.0 - exp(-3.8 * SQ(Config::c0 / c2));
    } else {
        values[1] = 1.0;
        values[2] = 1.0;
    }

    values[0] = 1.0;
}

bool ARBFInterpolator::m_isInTriangle(const double* x, int ta, int tb, int tc) {
    const Vertex &a = m_mesh->getVertices()[ta];
    const Vertex &b = m_mesh->getVertices()[tb];
    const Vertex &c = m_mesh->getVertices()[tc];

    Eigen::Vector2d ab(b.x - a.x, b.y - a.y);
    Eigen::Vector2d bc(c.x - b.x, c.y - b.y);
    Eigen::Vector2d ca(a.x - c.x, a.y - c.y);
    Eigen::Vector2d ax(x[0] - a.x, x[1] - a.y);
    Eigen::Vector2d bx(x[0] - b.x, x[1] - b.y);
    Eigen::Vector2d cx(x[0] - c.x, x[1] - c.y);

    double area = 0, area1 = 0, area2 = 0, area3 = 0;
    double cosTheta = ab.dot(-ca) / (ab.norm() * ca.norm());
    area = 0.5 * ab.norm() * ca.norm() * sqrt(1.0 - SQ(cosTheta));
    double cosTheta1 = ab.dot(ax) / (ab.norm() * ax.norm());
    area1 = 0.5 * ab.norm() * ax.norm() * sqrt(1.0 - SQ(cosTheta1));
    double cosTheta2 = cx.dot(ca) / (cx.norm() * ca.norm());
    area2 = 0.5 * cx.norm() * ca.norm() * sqrt(1.0 - SQ(cosTheta2));
    double cosTheta3 = bx.dot(bc) / (bx.norm() * bc.norm());
    area3 = 0.5 * bx.norm() * bc.norm() * sqrt(1.0 - SQ(cosTheta3));

    if (std::isnan(area)) area = 0.0;
    if (std::isnan(area1)) area1 = 0.0;
    if (std::isnan(area2)) area2 = 0.0;
    if (std::isnan(area3)) area3 = 0.0;

    return (std::abs(area1+area2+area3-area) < Config::epsilon);
}

double ARBFInterpolator::m_computeDistance(const double *v0, const double *v1) {
    return std::sqrt(SQ(v1[0]-v0[0]) + SQ(v1[1]-v0[1]) + SQ(v1[2]-v0[2]));
}

double ARBFInterpolator::m_computeDistance(const double *v0, const double *v1, const Eigen::Matrix3d & T) {
    Eigen::RowVector3d dis(v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]);
    Eigen::RowVector3d tmp = dis * T;

    // use definition: distance_T = dis * T * dis
    return tmp(0)*dis(0) + tmp(1)*dis(1) + tmp(2)*dis(2);
}

double ARBFInterpolator::m_computeDistance(const double *x, const double *v0, const double *v1) {
    double a, b, c;
    a = m_computeDistance(v0, v1);
    b = m_computeDistance(x, v0);
    c = m_computeDistance(x, v1);

    // [v0, v1] is not a line segment but a point
    if (a < Config::epsilon) {
        return b;
    }

    // x is on line segment [v0, v1]
    if (std::abs(b + c - a) < Config::epsilon) {
        return 0.0;
    }

    // [x, v0, v1] is a perpendicular or obtuse triangle, v0 is the perpendicular angle
    if (SQ(c) >= (SQ(a) + SQ(b))) {
        return b;
    }

    // [x, v0, v1] is a perpendicular or obtuse triangle, v1 is the perpendicular angle
    if (SQ(b) >= (SQ(a) + SQ(c))) {
        return c;
    }

    // [x, v0, v1] is a acute triangle, need to calculate triangle height
    double hc = (a + b + c) / 2.0; // half circumference
    double s = std::sqrt(hc * (hc - a) * (hc - b) * (hc - c)); // calculate area by Heron's formula
    return (s * 2.0) / a; // height of triangle
}

double ARBFInterpolator::m_computeDistance(const double *v0, const double *v1, const double *w0, const double *w1) {
    double dist = m_computeDistance(v0, w0, w1);
    dist = std::min(dist, m_computeDistance(v1, w0, w1));
    dist = std::min(dist, m_computeDistance(w0, v0, v1));
    dist = std::min(dist, m_computeDistance(w1, v0, v1));
    return dist;
}

std::vector<double> ARBFInterpolator::linspace(double a, double b, int n) {
    assert(b > a && n > 0);
    std::vector<double> res;
    double step = (b-a) / (n-1);

    while (a < b) {
        res.push_back(a);
        a += step;
    }

    // get last point
    if ((a-b) < Config::epsilon) {
        res.push_back(b);
    }

    return res;
}

void ARBFInterpolator::rescaleResultData(int low, int high) {
    assert((high - low) > Config::epsilon);

    // original data is in range [-1.0, 1.0]
    for (double &intensity: std::get<0>(m_result)) {
        intensity = (intensity - (-1)) / (1.0 - (-1.0)) * (high - low);
    }

//    for (int i = 0; i < NP; i++) {
//        intensities[i] = (intensities[i] + 1) / 2 * 255;
//    }
}

void ARBFInterpolator::thresholdResultData(double maxIntensity) {
    assert(maxIntensity > Config::epsilon);

    for (unsigned i = 0; i < std::get<0>(m_result).size(); i++) {
        if (std::get<0>(m_result)[i] < 0.0) {
            std::get<0>(m_result)[i] = 0.0;
        }

        if (std::get<0>(m_result)[i] > maxIntensity) {
            std::get<0>(m_result)[i] = maxIntensity;
        }
    }

//    for (auto &intensity: std::get<0>(m_result)) {
//        if (intensity < 0.0) {
//            intensity = 0.0;
//        }
//
//        if (intensity > maxIntensity) {
//            intensity = maxIntensity;
//        }
//    }

//        for (int j = 0; j < dimY; ++j) {
//            for (int i = 0; i < dimX; ++i) {
//                if (intensities[j*dimX + i] < 0.0) {
//                    intensities[j*dimX + i] = 0.0;
//                }
//
//                if (intensities[j*dimX + i] > MAX_INTENSITY) {
//                    intensities[j*dimX + i] = MAX_INTENSITY;
//                }
//            }
//        }
//    } else {
//        for (int k = 0; k < dimZ; ++k) {
//            for (int j = 0; j < dimY; ++j) {
//                for (int i = 0; i < dimX; ++i) {
//                    if (intensities[k*dimX*dimY + j*dimX + i] < 0.0) {
//                        intensities[k*dimX*dimY + j*dimX + i] = 0.0;
//                    }
//
//                    if (intensities[k*dimX*dimY + j*dimX + i] > MAX_INTENSITY) {
//                        intensities[k*dimX*dimY + j*dimX + i] = MAX_INTENSITY;
//                    }
//                }
//            }
//        }
//    }
}

bool ARBFInterpolator::m_isInTetrahedron(const double *x, int tetrahedronId) {
    TetMesh *mesh = dynamic_cast<TetMesh*>(m_mesh);
    const Tetrahedron &tet = mesh->getTetrahedrons()[tetrahedronId];
    int a = tet.a;
    int b = tet.b;
    int c = tet.c;
    int d = tet.d;
    const Vertex &v1 = mesh->getVertices()[a];
    const Vertex &v2 = mesh->getVertices()[b];
    const Vertex &v3 = mesh->getVertices()[c];
    const Vertex &v4 = mesh->getVertices()[d];

    Eigen::Matrix4d m0, m1, m2, m3, m4;
    m0 << v1.x, v1.y, v1.z, 1.0,
          v2.x, v2.y, v2.z, 1.0,
          v3.x, v3.y, v3.z, 1.0,
          v4.x, v4.y, v4.z, 1.0;

    m1 << x[0], x[1], x[2], 1.0,
          v2.x, v2.y, v2.z, 1.0,
          v3.x, v3.y, v3.z, 1.0,
          v4.x, v4.y, v4.z, 1.0;

    m2 << v1.x, v1.y, v1.z, 1.0,
          x[0], x[1], x[2], 1.0,
          v3.x, v3.y, v3.z, 1.0,
          v4.x, v4.y, v4.z, 1.0;

    m3 << v1.x, v1.y, v1.z, 1.0,
          v2.x, v2.y, v2.z, 1.0,
          x[0], x[1], x[2], 1.0,
          v4.x, v4.y, v4.z, 1.0;

    m4 << v1.x, v1.y, v1.z, 1.0,
          v2.x, v2.y, v2.z, 1.0,
          v3.x, v3.y, v3.z, 1.0,
          x[0], x[1], x[2], 1.0;

    double det0 = m0.determinant();
    double det1 = m1.determinant();
    double det2 = m2.determinant();
    double det3 = m3.determinant();
    double det4 = m4.determinant();
    return (det0 * det1 > Config::epsilon) && (det0 * det2 > Config::epsilon) && (det0 * det3 > Config::epsilon) &&
            (det0 * det4 > Config::epsilon);
}

bool ARBFInterpolator::m_isInHexahedron(const double *x, int hexahedronId) {
    HexMesh *mesh = dynamic_cast<HexMesh*>(m_mesh);
    const Hexahedron &hex = mesh->getHexahedrons()[hexahedronId];
    int a = hex.a;
    int b = hex.b;
    int c = hex.c;
    int d = hex.d;
    int e = hex.e;
    int f = hex.f;
    int g = hex.g;
    int h = hex.h;
    const Vertex &v1 = mesh->getVertices()[a];
    const Vertex &v2 = mesh->getVertices()[b];
    const Vertex &v3 = mesh->getVertices()[c];
    const Vertex &v4 = mesh->getVertices()[d];
    const Vertex &v5 = mesh->getVertices()[e];
    const Vertex &v6 = mesh->getVertices()[f];
    const Vertex &v7 = mesh->getVertices()[g];
    const Vertex &v8 = mesh->getVertices()[h];
    auto mmx = std::minmax({ v1.x, v2.x, v3.x, v4.x, v5.x, v6.x, v7.x, v8.x });
    auto mmy = std::minmax({ v1.y, v2.y, v3.y, v4.y, v5.y, v6.y, v7.y, v8.y });
    auto mmz = std::minmax({ v1.z, v2.z, v3.z, v4.z, v5.z, v6.z, v7.z, v8.z });
    return x[0] >= mmx.first && x[0] <= mmx.second &&
            x[1] >= mmy.first && x[1] <= mmy.second &&
            x[2] >= mmz.first && x[2] <= mmz.second;
}

bool ARBFInterpolator::m_isInAnyTetrahedron(const double *x) {
    TetMesh *mesh = dynamic_cast<TetMesh*>(m_mesh);

    for (unsigned i = 0; i < mesh->getNumTetrahedrons(); i++) {
        if (m_isInTetrahedron(x, i)) {
            return true;
        }
    }

    return false;
}

bool ARBFInterpolator::m_isInAnyHexahedron(const double *x) {
    HexMesh *mesh = dynamic_cast<HexMesh*>(m_mesh);

    for (unsigned i = 0; i < mesh->getNumHexahedrons(); i++) {
        if (m_isInHexahedron(x, i)) {
            return true;
        }
    }

    return false;
}

//Vertex ARBFInterpolator::createFalseCenter(const Vertex &trueCenter, const Vertex &vertex, double offsetRatio) {
//    int idx = std::rand() % m_samplingDirections.size();
//    double falseX = trueCenter.x + offsetRatio * m_samplingDirections[idx][0];
//    double falseY = trueCenter.y + offsetRatio * m_samplingDirections[idx][1];
//    double falseZ = trueCenter.z + offsetRatio * m_samplingDirections[idx][2];
//    Vertex falseCenter(falseX, falseY, falseZ, trueCenter.intensity);
//    double dx = vertex.x - falseCenter.x;
//    double dy = vertex.y - falseCenter.y;
//    double dz = vertex.z - falseCenter.z;
//    double norm = std::sqrt(SQ(dx) + SQ(dy) + SQ(dz));
//    double gradX = dx / norm;
//    double gradY = dy / norm;
//    double gradZ = dz / norm;
//    Eigen::Matrix3d mEdge;
//    mEdge(0, 0) = gradX * gradX;
//    mEdge(0, 1) = gradX * gradY;
//    mEdge(0, 2) = gradX * gradZ;
//    mEdge(1, 0) = gradY * gradX;
//    mEdge(1, 1) = gradY * gradY;
//    mEdge(1, 2) = gradY * gradZ;
//    mEdge(2, 0) = gradZ * gradX;
//    mEdge(2, 1) = gradZ * gradY;
//    mEdge(2, 2) = gradZ * gradZ;
//    std::cout << mEdge << std::endl;
//    Eigens3d eigens = computeEigens(mEdge); // compute eigenvalues and eigenvectors
//    Eigenvalues3d &eigenvalues = std::get<0>(eigens);
//    printf("Eigenvalues: %lf, %lf, %lf\n", eigenvalues[0], eigenvalues[1], eigenvalues[2]);
//    Eigenvectors3d &eigenvectors = std::get<1>(eigens);
//    rescaleEigenvalues(eigenvalues);
//    falseCenter.eig1 = eigenvalues[0];
//    falseCenter.eig2 = eigenvalues[1];
//    falseCenter.eig3 = eigenvalues[2];
//    falseCenter.eigVec1 = { eigenvectors[0][0], eigenvectors[0][1], eigenvectors[0][2] };
//    falseCenter.eigVec2 = { eigenvectors[1][0], eigenvectors[1][1], eigenvectors[1][2] };
//    falseCenter.eigVec3 = { eigenvectors[2][0], eigenvectors[2][1], eigenvectors[2][2] };
//    return falseCenter;
//}
