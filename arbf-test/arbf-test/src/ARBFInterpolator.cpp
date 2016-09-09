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
#include "../include/ARBFInterpolator.h"

ARBFInterpolator::ARBFInterpolator(): m_basis(nullptr), m_hasCoeffSolved(false) {}

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

//void ARBFInterpolator::calculateTriangleCenters() {
//    m_centers.resize(m_mesh->getNumFaces());
//
//    for (unsigned f = 0; f < m_mesh->getNumFaces(); ++f) {
//        int a = m_mesh->getFaces()[f].a;
//        int b = m_mesh->getFaces()[f].b;
//        int c = m_mesh->getFaces()[f].c;
//
//        // compute coordinates
//        m_centers[f].x = (m_mesh->getVertices()[a].x + m_mesh->getVertices()[b].x + m_mesh->getVertices()[c].x) / 3.0;
//        m_centers[f].y = (m_mesh->getVertices()[a].y + m_mesh->getVertices()[b].y + m_mesh->getVertices()[c].y) / 3.0;
//        m_centers[f].z = (m_mesh->getVertices()[a].z + m_mesh->getVertices()[b].z + m_mesh->getVertices()[c].z) / 3.0;
//
//        // compute intensity
//        m_centers[f].intensity = m_mesh->getFaces()[f].intensity;
//
//        // assign eigenvalues and eigenvectors
//        m_centers[f].eig1 = 0;
//        m_centers[f].eig2 = 0;
//        m_centers[f].eig3 = 0;
//        m_centers[f].eigVec1 = {1, 0, 0};
//        m_centers[f].eigVec2 = {0, 1, 0};
//        m_centers[f].eigVec3 = {0, 0, 1};
//    }
//
//    if (Config::isDebugEnabled) {
//        std::string filename = "DEBUG.centers.txt";
//        FILE *fp = nullptr;
//
//        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
//            fprintf(stderr, "WARNING: open DEBUG.centers.txt failed.\n");
//            return;
//        }
//
//        fprintf(fp, "X Y Z INTENSITY\n");
//
//        for (int f = 0; f < m_mesh->getNumFaces(); ++f) {
//            fprintf(fp, "%lf %lf %lf %lf\n", m_centers[f].x, m_centers[f].y, m_centers[f].z, m_centers[f].intensity);
//        }
//
//        fclose(fp);
//    }
//}

//void ARBFInterpolator::calculateEdgeCenters() {
//    m_edgeCenters.resize(m_mesh->getNumEdges());
//
//    int e = 0;
//    std::unordered_set<Edge>::const_iterator it;
//    for (it = m_mesh->getEdges().begin(); it != m_mesh->getEdges().end(); it++, e++) {
//        int a = it->a;
//        int b = it->b;
//        m_edgeCenters[e].x = m_mesh->getVertices()[a].x + 0.5 * (m_mesh->getVertices()[b].x - m_mesh->getVertices()[a].x);
//        m_edgeCenters[e].y = m_mesh->getVertices()[a].y + 0.5 * (m_mesh->getVertices()[b].y - m_mesh->getVertices()[a].y);
//        m_edgeCenters[e].z = m_mesh->getVertices()[a].z + 0.5 * (m_mesh->getVertices()[b].z - m_mesh->getVertices()[a].z);
//        m_edgeCenters[e].intensity = it->intensity;
//
//        double dx = m_mesh->getVertices()[b].x - m_mesh->getVertices()[a].x;
//        double dy = m_mesh->getVertices()[b].y - m_mesh->getVertices()[a].y;
//        double dz = m_mesh->getVertices()[b].z - m_mesh->getVertices()[a].z;
//        double norm = std::sqrt(SQ(dx) + SQ(dy) + SQ(dz));
//        double gradX = dx / norm;
//        double gradY = dy / norm;
//        double gradZ = dz / norm;
//        Eigen::Matrix3d mEdge;
////        Eigen::Matrix2d mEdge;
//        mEdge(0, 0) = gradX * gradX;
//        mEdge(0, 1) = gradX * gradY;
//        mEdge(0, 2) = gradX * gradZ;
//        mEdge(1, 0) = gradY * gradX;
//        mEdge(1, 1) = gradY * gradY;
//        mEdge(1, 2) = gradY * gradZ;
//        mEdge(2, 0) = gradZ * gradX;
//        mEdge(2, 1) = gradZ * gradY;
//        mEdge(2, 2) = gradZ * gradZ;
//
//        Eigens3d eigens = computeEigens(mEdge); // compute eigenvalues and eigenvectors
//        Eigenvalues3d &eigenvalues = std::get<0>(eigens);
//        Eigenvectors3d &eigenvectors = std::get<1>(eigens);
//        rescaleEigenvalues(eigenvalues);
//
//        m_edgeCenters[e].eig1 = eigenvalues[0];
//        m_edgeCenters[e].eig2 = eigenvalues[1];
//        m_edgeCenters[e].eig3 = eigenvalues[2];
//        m_edgeCenters[e].eigVec1 = { eigenvectors[0][0], eigenvectors[0][1], eigenvectors[0][2] };
//        m_edgeCenters[e].eigVec2 = { eigenvectors[1][0], eigenvectors[1][1], eigenvectors[1][2] };
//        m_edgeCenters[e].eigVec3 = { eigenvectors[2][0], eigenvectors[2][1], eigenvectors[2][2] };
//    }
//
//    //    PPMImage im(301, 301);
//    //    im.setMaxIntensity(255);
//    //    im.setPath("/Users/keliu/Documents/projects/arbf/arbf-test/arbf-test/edge.ppm");
//    //    double* data = new double[301*301];
//    //    memset(data, 0, sizeof(double) * 301 * 301);
//    //    for (int j = 0; j < 301; j++) {
//    //        for (int i = 0; i < 301; i++) {
//    //            for (int e = 0;  e < m_mesh->getNumEdges(); e++) {
//    //                if (m_edgeCenters[e].x * 100 == i && m_edgeCenters[e].y * 100 == j) {
//    //                    data[j*301+i] = 255;
//    //                }
//    //            }
//    //        }
//    //    }
//    //    im.setImageData(data);
//    //    im.write();
//
//    if (Config::isDebugEnabled) {
//        std::string filename = "DEBUG.edge_centers.txt";
//        FILE *fp = nullptr;
//
//        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
//            fprintf(stderr, "WARNING: open DEBUG.edge_centers.txt failed.\n");
//            return;
//        }
//
//        fprintf(fp, "X Y INTENSITY\n");
//
//        for (auto center: m_edgeCenters) {
//            fprintf(fp, "%lf %lf %lf %lf\n", center.x, center.y, center.z, center.intensity);
//        }
//
//        fclose(fp);
//    }
//}

//void ARBFInterpolator::calculateTetrahedronCenters() {
//    double tmpx = 0, tmpy = 0, tmpz = 0;
//    TetMesh *mesh = dynamic_cast<TetMesh*>(m_mesh);
//
//    for (int t = 0; t < mesh->getNumTetrahedrons(); t++) {
//        int a = mesh->getTetrahedrons()[t].a;
//        int b = mesh->getTetrahedrons()[t].b;
//        int c = mesh->getTetrahedrons()[t].c;
//        int d = mesh->getTetrahedrons()[t].d;
//        tmpx += (mesh->getVertices()[a].x + mesh->getVertices()[b].x + mesh->getVertices()[c].x + mesh->getVertices()[d].x);
//        tmpy += (mesh->getVertices()[a].y + mesh->getVertices()[b].y + mesh->getVertices()[c].y + mesh->getVertices()[d].y);
//        tmpz += (mesh->getVertices()[a].z + mesh->getVertices()[b].z + mesh->getVertices()[c].z + mesh->getVertices()[d].z);
//        Vertex center(tmpx / 4.0, tmpy / 4.0, tmpz / 4.0, -1.0);
////        double ratio = 0;
////        std::srand((unsigned) time(NULL));
////        Vertex falseCenter1 = createFalseCenter(center, m_centers[0], ratio);
////    Vertex falseCenter2 = createFalseCenter(center, m_centers[1], ratio);
////    Vertex falseCenter3 = createFalseCenter(center, m_centers[2], ratio);
////    Vertex falseCenter4 = createFalseCenter(center, m_centers[3], ratio);
//        m_tetrahedronCenters.push_back(center);
////    m_tetrahedronCenters.push_back(falseCenter2);
////    m_tetrahedronCenters.push_back(falseCenter3);
////    m_tetrahedronCenters.push_back(falseCenter4);
//        tmpx = tmpy = tmpz = 0;
//    }
//}

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

void ARBFInterpolator::interpolate(unsigned numEvalPoints) {
    if (!m_hasCoeffSolved) {
        fprintf(stderr, "Error: coefficients are not solved out yet.");
        exit(EXIT_FAILURE);
    }

    unsigned dimX = 0, dimY = 0, dimZ = 0;
    std::vector<double> x, y, z;
    unsigned nv = m_mesh->getNumVertices();
    unsigned nf = m_mesh->getNumFaces();
    unsigned ne = m_mesh->getNumEdges();
    auto faces = m_mesh->getFacesAsList();
    auto edges = m_mesh->getEdgesAsList();
    unsigned dim = 0; // dimension of distance matrix
    std::vector<double> intensities(numEvalPoints, -1.0); // solution

    switch (Config::problemDim) {
        case 2:
            dimX = (unsigned) std::sqrt(numEvalPoints);
            dimY = (unsigned) std::sqrt(numEvalPoints);
            x = linspace(m_mesh->getMinX(), m_mesh->getMaxX(), dimX);
            y = linspace(m_mesh->getMinY(), m_mesh->getMaxY(), dimY);
            dim = nv + nf + ne;

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

//            std::cout << "\n----- Intensity at face centers -----\n";
//            for (int j = 0; j < dimY; j++) {
//                for (int i = 0; i < dimX; i++) {
//                    for (int v = 0; v < m_mesh->getNumVertices(); v++) {
//                        if (std::abs(x[i] - m_mesh->getVertices()[v].x) < 1e-6 && std::abs(y[j] - m_mesh->getVertices()[v].y) < 1e-6) {
//                            printf("\tAt vert[%d][%lf, %lf], intensity[%d, %d][%lf, %lf] = %lf\n",
//                                   v, m_mesh->getVertices()[v].x, m_mesh->getVertices()[v].y, i, j,
//                                   x[i], y[j], intensities[j*dimX+i]);
//                            break;
//                        }
//                    }
//
//                    for (int f = 0; f < m_mesh->getNumFaces(); f++) {
//                        if (std::abs(x[i] - getTriangleCenters()[f].x) < 1e-6 && std::abs(y[j] - getTriangleCenters()[f].y) < 1e-6) {
//                            printf("\tAt face[%d][%lf, %lf], intensity[%d, %d][%lf, %lf] = %lf\n",
//                                   f, getTriangleCenters()[f].x, getTriangleCenters()[f].y, i, j,
//                                   x[i], y[j], intensities[j*dimX+i]);
//                            break;
//                        }
//                    }
//
//                    for (int e = 0; e < m_mesh->getNumEdges(); e++) {
//                        if (std::abs(x[i] - getEdgeCenters()[e].x) < 1e-6 && std::abs(y[j] - getEdgeCenters()[e].y) < 1e-6) {
//                            printf("\tAt edgeCenter[%d][%lf, %lf], intensity[%d, %d][%lf, %lf] = %lf\n",
//                                   e, getEdgeCenters()[e].x, getEdgeCenters()[e].y, i, j,
//                                   x[i], y[j], intensities[j*dimX+i]);
//                            break;
//                        }
//                    }
//                }
//            }
//            std::cout << "\n" << "-----" << std::endl;
            break;
        case 3:
        default:
            TetMesh *mesh = dynamic_cast<TetMesh*>(m_mesh);
            unsigned nt = mesh->getNumTetrahedrons();
            dimX = (unsigned) std::round(std::pow(numEvalPoints, 1.0/3.0));
            dimY = (unsigned) std::round(std::pow(numEvalPoints, 1.0/3.0));
            dimZ = (unsigned) std::round(std::pow(numEvalPoints, 1.0/3.0));
            x = linspace(m_mesh->getMinX(), m_mesh->getMaxX(), dimX);
            y = linspace(m_mesh->getMinY(), m_mesh->getMaxY(), dimY);
            z = linspace(m_mesh->getMinZ(), m_mesh->getMaxZ(), dimZ);
            dim = nv + nf + nt + (nt * 4);
            double c1 = 0.1;
            double p0[] = { 0, 0, 0 }, p1[] = { 0, 0, 0 };
            double v0[3], v1[3], v2[3], v3[3], v4[3], w0[3], w1[3] ,w2[3], w3[3], w4[3];

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
                        if (! m_isInAnyTetrahedron(p0)) {
                            continue;
                        }

                        double intensity = 0;
                        for (int l = 0; l < dim;) {
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
                        intensities[k*dimX*dimY + j*dimX + i] = intensity;
                    }
                }
            }

            std::cout << "\n----- Intensity at face centers -----\n";
            for (int k = 0; k < dimZ; k++) {
                for (int j = 0; j < dimY; j++) {
                    for (int i = 0; i < dimX; i++) {
                        for (int v = 0; v < m_mesh->getNumVertices(); v++) {
                            if (std::abs(x[i] - m_mesh->getVertices()[v].x) < 1e-3 &&
                                std::abs(y[j] - m_mesh->getVertices()[v].y) < 1e-3 &&
                                std::abs(z[k] - m_mesh->getVertices()[v].z) < 1e-3)
                            {
                                printf("\tAt vertex[%d][%lf, %lf, %lf], intensity[%d, %d, %d][%lf, %lf, %lf] = %lf\n",
                                       v, m_mesh->getVertices()[v].x, m_mesh->getVertices()[v].y, m_mesh->getVertices()[v].z,
                                       i, j, k, x[i], y[j], z[k], intensities[k*dimY*dimX+j*dimX+i]);
                                break;
                            }
                        }

                        for (int f = 0; f < m_mesh->getNumFaces(); f++) {
                            if (std::abs(x[i] - faces[f].center.x) < 1e-3 &&
                                    std::abs(y[j] - faces[f].center.y) < 1e-3 &&
                                    std::abs(z[k] - faces[f].center.z) < 1e-3)
                            {
                                printf("\tAt face[%d][%lf, %lf, %lf], intensity[%d, %d, %d][%lf, %lf, %lf] = %lf\n",
                                       f, faces[f].center.x, faces[f].center.y, faces[f].center.z,
                                       i, j, k, x[i], y[j], z[k], intensities[k*dimY*dimX+j*dimX+i]);
                                break;
                            }
                        }

                        for (unsigned c = 0; c < mesh->getNumTetrahedrons(); c++) {
                            if (std::abs(x[i] - mesh->getTetrahedrons()[c].center.x) < 1e-3 &&
                                std::abs(y[j] - mesh->getTetrahedrons()[c].center.y) < 1e-3 &&
                                std::abs(z[k] - mesh->getTetrahedrons()[c].center.z) < 1e-3)
                            {
                                printf("\tAt false center[%d][%lf, %lf, %lf], intensity[%d, %d, %d][%lf, %lf, %lf] = %lf\n",
                                       c, mesh->getTetrahedrons()[c].center.x, mesh->getTetrahedrons()[c].center.y, mesh->getTetrahedrons()[c].center.z,
                                       i, j, k, x[i], y[j], z[k], intensities[k*dimY*dimX+j*dimX+i]);
                                break;
                            }
                        }
                    }
                }
            }
            std::cout << "\n" << "-----" << std::endl;
            break;
    }


    //    cout << "\nma = \n";
    //    cout << ma << "\n" << "-----" << endl;
    //    cout << ma - ma.transpose() << "-----" << endl;
    //    cout << u << "\n" << "-----" << endl;

    // coeff = ma.inverse() * u;
    // cout << "\ncoeff = \n" << coeff << "\n" << "-----" << endl;
    // cout << "\ncoeff.sum() = " << coeff.sum() << endl;

    std::get<0>(m_result) = intensities;
    rescaleResultData(0, 255);
    double maxIntensity = 255.0;
    thresholdResultData(maxIntensity);
    
    //    cout << "\nintensity after rescale = \n";
    //    for (int j = 0; j < dimY; j++) {
    //        for (int i = 0; i < dimX; i++) {
    //            cout << (int) round(intensities[j*dimX+i]) << " ";
    //        }
    //        cout << "\n";
    //    }
    //    cout << "\n" << "-----" << endl;

    std::get<1>(m_result) = maxIntensity;
    std::get<2>(m_result) = { dimX, dimY, dimZ };
//    return std::make_tuple(intensities, maxIntensity, { dimX, dimY, dimZ });
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

void ARBFInterpolator::solveDistanceMatrix2d() {
    unsigned nv = m_mesh->getNumVertices();
    unsigned nf = m_mesh->getNumFaces();
    unsigned ne = m_mesh->getNumEdges();
    unsigned matrixDimension = nv + nf + ne;

    m_distanceMatrix.resize(matrixDimension, matrixDimension);
    m_u.resize(matrixDimension);
    const Vertex *vi = nullptr, *vj = nullptr;
    Eigen::Matrix3d tensor;
    auto faces = m_mesh->getFacesAsList();
    auto edges = m_mesh->getEdgesAsList();

    for (unsigned i = 0; i < matrixDimension; i++) {
        if (i < nv) {
            vi = &m_mesh->getVertices()[i];
        } else if (i < (nv + nf)) {
            vi = &faces[i - nv].center;
        } else {
            vi = &edges[i - nv - nf].center;
        }

        m_u(i) = vi->intensity;
        double p0[] = { vi->x, vi->y, vi->z };

        for (int j = 0; j < matrixDimension; j++) {
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
    m_hasCoeffSolved = true;
}

void ARBFInterpolator::solveDistanceMatrix3d() {
    TetMesh *mesh = dynamic_cast<TetMesh*>(m_mesh);
    unsigned nv = mesh->getNumVertices();
    unsigned nf = mesh->getNumFaces();
    unsigned nt = mesh->getNumTetrahedrons();
    unsigned rows = nv + nf + nt + (nt * 4);
    unsigned cols = rows;
    m_distanceMatrix.resize(rows, cols);
    m_u.resize(rows);
    auto faces = mesh->getFacesAsList();
    double c1 = 0.1;
    double p0[] = { 0, 0, 0 }, p1[] = { 0, 0, 0 };
    double v0[3], v1[3], v2[3], v3[3], v4[3], w0[3], w1[3] ,w2[3], w3[3], w4[3];

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
                m_distanceMatrix(i, j+1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i, j+2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i, j+3) = m_basis->phi(r4, c1);
                j += 4;
            }
        }
    }

    for (unsigned i = (nv + nf + nt); i < rows; i+=4) {
        m_u(i) = -1;
        m_u(i+1) = -1;
        m_u(i+2) = -1;
        m_u(i+3) = -1;

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
                m_distanceMatrix(i+1, j) = m_basis->phi(r2, c1);
                m_distanceMatrix(i+2, j) = m_basis->phi(r3, c1);
                m_distanceMatrix(i+3, j) = m_basis->phi(r4, c1);
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
                m_distanceMatrix(i+1, j) = m_basis->phi(r2, c1);
                m_distanceMatrix(i+2, j) = m_basis->phi(r3, c1);
                m_distanceMatrix(i+3, j) = m_basis->phi(r4, c1);
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
                m_distanceMatrix(i+1, j) = m_basis->phi(r2, c1);
                m_distanceMatrix(i+2, j) = m_basis->phi(r3, c1);
                m_distanceMatrix(i+3, j) = m_basis->phi(r4, c1);
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
                m_distanceMatrix(i, j+1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i, j+2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i, j+3) = m_basis->phi(r4, c1);

                r1 = m_computeDistance(v0, v2, w0, w1);
                r2 = m_computeDistance(v0, v2, w0, w2);
                r3 = m_computeDistance(v0, v2, w0, w3);
                r4 = m_computeDistance(v0, v2, w0, w4);
                m_distanceMatrix(i+1, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i+1, j+1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i+1, j+2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i+1, j+3) = m_basis->phi(r4, c1);

                r1 = m_computeDistance(v0, v3, w0, w1);
                r2 = m_computeDistance(v0, v3, w0, w2);
                r3 = m_computeDistance(v0, v3, w0, w3);
                r4 = m_computeDistance(v0, v3, w0, w4);
                m_distanceMatrix(i+2, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i+2, j+1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i+2, j+2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i+2, j+3) = m_basis->phi(r4, c1);

                r1 = m_computeDistance(v0, v4, w0, w1);
                r2 = m_computeDistance(v0, v4, w0, w2);
                r3 = m_computeDistance(v0, v4, w0, w3);
                r4 = m_computeDistance(v0, v4, w0, w4);
                m_distanceMatrix(i+3, j) = m_basis->phi(r1, c1);
                m_distanceMatrix(i+3, j+1) = m_basis->phi(r2, c1);
                m_distanceMatrix(i+3, j+2) = m_basis->phi(r3, c1);
                m_distanceMatrix(i+3, j+3) = m_basis->phi(r4, c1);

                j += 4;
            }
        }
    }

    m_coeff = m_distanceMatrix.fullPivLu().solve(m_u); // call Eigen library to solve
//    std::cout << m_distanceMatrix << std::endl;

    std::cout << m_coeff << std::endl;

    m_hasCoeffSolved = true;
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

bool ARBFInterpolator::m_isInAnyTetrahedron(const double *x) {
    TetMesh *mesh = dynamic_cast<TetMesh*>(m_mesh);

    for (unsigned i = 0; i < mesh->getNumTetrahedrons(); i++) {
        if (m_isInTetrahedron(x, i)) {
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
