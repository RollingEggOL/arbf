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
#include <iostream>
#include <cfenv> // for std::fesetround()
#include <tuple>
#include <vector>
#include <unordered_set>
#include <Eigen/Eigenvalues>
#include "../include/Config.h"
#include "../include/ARBFInterpolator.h"

// ARBFInterpolator implementations
ARBFInterpolator::ARBFInterpolator(): m_basis(nullptr), m_hasCoeffSolved(false) {}

ARBFInterpolator::~ARBFInterpolator() {
//    delete [] m_centers;
//    m_centers = nullptr;
//    delete [] m_edgeCenters;
//    m_edgeCenters = nullptr;
//    delete [] m_bodyCenters;
//    m_bodyCenters = nullptr;
}

const std::vector<Vertex>& ARBFInterpolator::getTriangleCenters() const {
    return m_centers;
}

const std::vector<Vertex>& ARBFInterpolator::getEdgeCenters() const {
    return m_edgeCenters;
}

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

void ARBFInterpolator::calculateTriangleCenters() {
    m_centers.resize(m_mesh->getNumFaces());

    for (unsigned f = 0; f < m_mesh->getNumFaces(); ++f) {
        int a = m_mesh->getFaces()[f].a;
        int b = m_mesh->getFaces()[f].b;
        int c = m_mesh->getFaces()[f].c;

        // compute coordinates
        m_centers[f].x = (m_mesh->getVertices()[a].x + m_mesh->getVertices()[b].x + m_mesh->getVertices()[c].x) / 3.0;
        m_centers[f].y = (m_mesh->getVertices()[a].y + m_mesh->getVertices()[b].y + m_mesh->getVertices()[c].y) / 3.0;
        m_centers[f].z = (m_mesh->getVertices()[a].z + m_mesh->getVertices()[b].z + m_mesh->getVertices()[c].z) / 3.0;

        // compute intensity
        m_centers[f].intensity = m_mesh->getFaces()[f].intensity;

        // assign eigenvalues and eigenvectors
        m_centers[f].eig1 = 0;
        m_centers[f].eig2 = 0;
        m_centers[f].eig3 = 0;
        m_centers[f].eigVec1 = {1, 0, 0};
        m_centers[f].eigVec2 = {0, 1, 0};
        m_centers[f].eigVec3 = {0, 0, 1};
    }

    if (Config::isDebugEnabled) {
        std::string filename = "DEBUG.centers.txt";
        FILE *fp = nullptr;

        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
            fprintf(stderr, "WARNING: open DEBUG.centers.txt failed.\n");
            return;
        }

        fprintf(fp, "X Y Z INTENSITY\n");

        for (int f = 0; f < m_mesh->getNumFaces(); ++f) {
            fprintf(fp, "%lf %lf %lf %lf\n", m_centers[f].x, m_centers[f].y, m_centers[f].z, m_centers[f].intensity);
        }

        fclose(fp);
    }
}

void ARBFInterpolator::calculateEdgeCenters() {
    m_edgeCenters.resize(m_mesh->getNumEdges());

    int e = 0;
    std::unordered_set<Edge, edgeHasher, edgeComparator>::const_iterator it;
    for (it = m_mesh->getEdges().begin(); it != m_mesh->getEdges().end(); it++, e++) {
        int a = it->a;
        int b = it->b;
        m_edgeCenters[e].x = m_mesh->getVertices()[a].x + 0.5 * (m_mesh->getVertices()[b].x - m_mesh->getVertices()[a].x);
        m_edgeCenters[e].y = m_mesh->getVertices()[a].y + 0.5 * (m_mesh->getVertices()[b].y - m_mesh->getVertices()[a].y);
        m_edgeCenters[e].z = m_mesh->getVertices()[a].z + 0.5 * (m_mesh->getVertices()[b].z - m_mesh->getVertices()[a].z);
        m_edgeCenters[e].intensity = it->intensity;

        double dx = m_mesh->getVertices()[b].x - m_mesh->getVertices()[a].x;
        double dy = m_mesh->getVertices()[b].y - m_mesh->getVertices()[a].y;
        double dz = m_mesh->getVertices()[b].z - m_mesh->getVertices()[a].z;
        double norm = std::sqrt(SQ(dx) + SQ(dy) + SQ(dz));
        double gradX = dx / norm;
        double gradY = dy / norm;
        double gradZ = dz / norm;
        Eigen::Matrix3d mEdge;
//        Eigen::Matrix2d mEdge;
        mEdge(0, 0) = gradX * gradX;
        mEdge(0, 1) = gradX * gradY;
        mEdge(0, 2) = gradX * gradZ;
        mEdge(1, 0) = gradY * gradX;
        mEdge(1, 1) = gradY * gradY;
        mEdge(1, 2) = gradY * gradZ;
        mEdge(2, 0) = gradZ * gradX;
        mEdge(2, 1) = gradZ * gradY;
        mEdge(2, 2) = gradZ * gradZ;

        Eigens3d eigens = computeEigens(mEdge); // compute eigenvalues and eigenvectors
        Eigenvalues3d &eigenvalues = std::get<0>(eigens);
        Eigenvectors3d &eigenvectors = std::get<1>(eigens);
        rescaleEigenvalues(eigenvalues);

        m_edgeCenters[e].eig1 = eigenvalues[0];
        m_edgeCenters[e].eig2 = eigenvalues[1];
        m_edgeCenters[e].eig3 = eigenvalues[2];
        m_edgeCenters[e].eigVec1 = { eigenvectors[0][0], eigenvectors[0][1], eigenvectors[0][2] };
        m_edgeCenters[e].eigVec2 = { eigenvectors[1][0], eigenvectors[1][1], eigenvectors[1][2] };
        m_edgeCenters[e].eigVec3 = { eigenvectors[2][0], eigenvectors[2][1], eigenvectors[2][2] };
    }

    //    PPMImage im(301, 301);
    //    im.setMaxIntensity(255);
    //    im.setPath("/Users/keliu/Documents/projects/arbf/arbf-test/arbf-test/edge.ppm");
    //    double* data = new double[301*301];
    //    memset(data, 0, sizeof(double) * 301 * 301);
    //    for (int j = 0; j < 301; j++) {
    //        for (int i = 0; i < 301; i++) {
    //            for (int e = 0;  e < m_mesh->getNumEdges(); e++) {
    //                if (m_edgeCenters[e].x * 100 == i && m_edgeCenters[e].y * 100 == j) {
    //                    data[j*301+i] = 255;
    //                }
    //            }
    //        }
    //    }
    //    im.setImageData(data);
    //    im.write();

    if (Config::isDebugEnabled) {
        std::string filename = "DEBUG.edge_centers.txt";
        FILE *fp = nullptr;

        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
            fprintf(stderr, "WARNING: open DEBUG.edge_centers.txt failed.\n");
            return;
        }

        fprintf(fp, "X Y INTENSITY\n");

        for (auto center: m_edgeCenters) {
            fprintf(fp, "%lf %lf %lf %lf\n", center.x, center.y, center.z, center.intensity);
        }

        fclose(fp);
    }
}

void ARBFInterpolator::computeEdgeMetrics() {
    m_T.resize(m_mesh->getNumEdges());

    for (int e = 0; e < m_mesh->getNumEdges(); e++) {
        const Vertex& edgeCenter = m_edgeCenters[e];
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
        m_T[e] = t1 * t2;
        m_T[e] = m_T[e] * t1.transpose();
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

//    const unsigned NP = 160000; // number of evaluation points
    unsigned dimX = 0, dimY = 0, dimZ = 0;
    std::vector<double> x, y, z;
    unsigned nv = m_mesh->getNumVertices(), nf = m_mesh->getNumFaces(), ne = m_mesh->getNumEdges();
    unsigned dim = 0; // dimension of distance matrix
    std::vector<double> intensities(numEvalPoints, 0.0); // solution

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
                    double intensity = 0;

                    for (unsigned l = 0; l < dim; l++) {
                        double r = 0.0;

                        if (l < nv) {
                            double p1[] = { m_mesh->getVertices()[l].x, m_mesh->getVertices()[l].y, m_mesh->getVertices()[l].z };
                            r = m_computeDistance(p0, p1, Eigen::Matrix3d::Identity());
                        } else if (l < (nv + nf)) {
                            double p1[] = { m_centers[l - nv].x, m_centers[l - nv].y, m_centers[l - nv].z };
                            r = m_computeDistance(p0, p1, Eigen::Matrix3d::Identity());
                        } else {
                            double p1[] = { m_edgeCenters[l - nv - nf].x, m_edgeCenters[l - nv - nf].y, m_edgeCenters[l - nv - nf].z };
                            r = m_computeDistance(p0, p1, m_T[l - nv - nf]);
                        }

                        intensity += (m_basis->phi(r) * m_coeff(l));
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
            dimX = (unsigned) std::round(std::pow(numEvalPoints, 1.0/3.0));
            dimY = (unsigned) std::round(std::pow(numEvalPoints, 1.0/3.0));
            dimZ = (unsigned) std::round(std::pow(numEvalPoints, 1.0/3.0));
            x = linspace(m_mesh->getMinX(), m_mesh->getMaxX(), dimX);
            y = linspace(m_mesh->getMinY(), m_mesh->getMaxY(), dimY);
            z = linspace(m_mesh->getMinZ(), m_mesh->getMaxZ(), dimZ);
            dim = nv + nf + 1;

            //    cout << "\nx = ";
            //    copy(x.begin(), x.end(), ostream_iterator<double>(cout, ", "));
            //    cout << "\ny = ";
            //    copy(y.begin(), y.end(), ostream_iterator<double>(cout, ", "));
            //    cout << "\nz = ";
            //    copy(z.begin(), z.end(), ostream_iterator<double>(cout, ", "));

            double centerX = (m_mesh->getMinX() + m_mesh->getMaxX()) / 2.0,
                    centerY = (m_mesh->getMinY() + m_mesh->getMaxY()) / 2.0,
                    centerZ = (m_mesh->getMinZ() + m_mesh->getMaxZ()) / 2.0;

            for (int k = 0; k < dimZ; k++) {
                for (int j = 0; j < dimY; j++) {
                    for (int i = 0; i < dimX; i++) {
                        double p0[] = { x[i], y[j], z[k] };
                        double intensity = 0;
                        for (int l = 0; l < dim; l++) {
                            double p1[] = { centerX, centerY, centerZ };

                            if (l < nv) {
                                p1[0] = m_mesh->getVertices()[l].x;
                                p1[1] = m_mesh->getVertices()[l].y;
                                p1[2] = m_mesh->getVertices()[l].z;
                            } else if (l < (nv + nf)) {
                                p1[0] = m_centers[l - nv].x;
                                p1[1] = m_centers[l - nv].y;
                                p1[2] = m_centers[l - nv].z;
                            }

                            double r = m_computeDistance(p0, p1, Eigen::Matrix3d::Identity());
                            intensity += (m_basis->phi(r) * m_coeff(l));
                        }
                        intensities[k*dimX*dimY + j*dimX + i] = intensity;
                    }
                }
            }

//            std::cout << "\n----- Intensity at face centers -----\n";
//            for (int k = 0; k < dimZ; k++) {
//                for (int j = 0; j < dimY; j++) {
//                    for (int i = 0; i < dimX; i++) {
//                        for (int v = 0; v < m_mesh->getNumVertices(); v++) {
//                            if (std::abs(x[i] - m_mesh->getVertices()[v].x) < 1e-6 &&
//                                std::abs(y[j] - m_mesh->getVertices()[v].y) < 1e-6 &&
//                                std::abs(z[k] - m_mesh->getVertices()[v].z) < 1e-6)
//                            {
//                                printf("\tAt vert[%d][%lf, %lf, %lf], intensity[%d, %d, %d][%lf, %lf, %lf] = %lf\n",
//                                       v, m_mesh->getVertices()[v].x, m_mesh->getVertices()[v].y, m_mesh->getVertices()[v].z, i, j, k,
//                                       x[i], y[j], z[k], intensities[k*dimY*dimX+j*dimX+i]);
//                                break;
//                            }
//                        }
//
//                        for (int f = 0; f < m_mesh->getNumFaces(); f++) {
//                            if (std::abs(x[i] - getTriangleCenters()[f].x) < 1e-6 &&
//                                std::abs(y[j] - getTriangleCenters()[f].y) < 1e-6 &&
//                                std::abs(z[k] - getTriangleCenters()[f].z) < 1e-6)
//                            {
//                                printf("\tAt face[%d][%lf, %lf, %lf], intensity[%d, %d, %d][%lf, %lf, %lf] = %lf\n",
//                                       f, getTriangleCenters()[f].x, getTriangleCenters()[f].y, getTriangleCenters()[f].z, i, j, k,
//                                       x[i], y[j], z[k], intensities[k*dimY*dimX+j*dimX+i]);
//                                break;
//                            }
//                        }
//
//                        if (std::abs(x[i] - (m_mesh->getMinX()+m_mesh->getMaxX())/2.0) < 1e-6 &&
//                            std::abs(y[j] - (m_mesh->getMinY()+m_mesh->getMaxY())/2.0) < 1e-6 &&
//                            std::abs(z[k] - (m_mesh->getMinZ()+m_mesh->getMaxZ())/2.0) < 1e-6)
//                        {std::
//                            printf("\tAt center[%lf, %lf, %lf], intensity[%d, %d, %d][%lf, %lf, %lf] = %lf\n",
//                                   (m_mesh->getMinX()+m_mesh->getMaxX())/2.0,
//                                   (m_mesh->getMinY()+m_mesh->getMaxY())/2.0,
//                                   (m_mesh->getMinZ()+m_mesh->getMaxZ())/2.0, i, j, k,
//                                   x[i], y[j], z[k], intensities[k*dimY*dimX+j*dimX+i]);
//                        }
//                    }
//                }
//            }
//            std::cout << "\n" << "-----" << std::endl;
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
    auto vectors = solver.eigenvectors();
    std::array<double, 3> vec1 = { vectors(0, 0).real(), vectors(0, 1).real(), vectors(0, 2).real() };
    std::array<double, 3> vec2 = { vectors(1, 0).real(), vectors(1, 1).real(), vectors(1, 2).real() };
    std::array<double, 3> vec3 = { vectors(2, 0).real(), vectors(2, 1).real(), vectors(2, 2).real() };
    Eigenvectors3d eigenvectors = { vec1, vec2, vec3 };
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

bool ARBFInterpolator::m_isInTriangle(const double* x, int triangleId) {
    const Vertex &a = m_mesh->getVertices()[m_mesh->getFaces()[triangleId].a];
    const Vertex &b = m_mesh->getVertices()[m_mesh->getFaces()[triangleId].b];
    const Vertex &c = m_mesh->getVertices()[m_mesh->getFaces()[triangleId].c];

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

    for (unsigned i = 0; i < matrixDimension; i++) {
        if (i < nv) {
            vi = &m_mesh->getVertices()[i];
        } else if (i < (nv + nf)) {
            vi = &getTriangleCenters()[i - nv];
        } else {
            vi = &getEdgeCenters()[i - nv - nf];
        }

        m_u(i) = vi->intensity;
        double p0[] = { vi->x, vi->y, vi->z };

        for (int j = 0; j < matrixDimension; j++) {
            if (j < nv) {
                vj = &m_mesh->getVertices()[j];
                tensor = Eigen::Matrix3d::Identity();
            } else if (j < (nv + nf)) {
                vj = &getTriangleCenters()[j - nv];
                tensor = Eigen::Matrix3d::Identity();
            } else {
                vj = &getEdgeCenters()[j - nv - nf];
                tensor = m_T[j-nv-nf];
            }

            double p1[] = { vj->x, vj->y, vj->z };
            double r = m_computeDistance(p0, p1, tensor);
            m_distanceMatrix(i, j) = m_basis->phi(r);
        };
    }

    m_coeff = m_distanceMatrix.fullPivLu().solve(m_u); // call Eigen library to solve
    m_hasCoeffSolved = true;
}

void ARBFInterpolator::solveDistanceMatrix3d() {
    unsigned nv = m_mesh->getNumVertices();
    unsigned nf = m_mesh->getNumFaces();

    unsigned matrixDimension = nv + nf + 1; // 1 means tetrahedron center

    m_distanceMatrix.resize(matrixDimension, matrixDimension);
    m_u.resize(matrixDimension);
    const Vertex *vi = nullptr, *vj = nullptr;
    Eigen::Matrix3d tensor = Eigen::Matrix3d::Identity(); // use RBF for 3D

    // tetrahedron center
    double centerX = (m_mesh->getMinX() + m_mesh->getMaxX()) / 2.0;
    double centerY = (m_mesh->getMinY() + m_mesh->getMaxY()) / 2.0;
    double centerZ = (m_mesh->getMinZ() + m_mesh->getMaxZ()) / 2.0;

    for (int i = 0; i < matrixDimension; i++) {
        double p0[] = { centerX, centerY, centerZ };
        m_u(i) = -1;

        if (i < nv) {
            vi = &m_mesh->getVertices()[i];
            p0[0] = vi->x;
            p0[1] = vi->y;
            p0[2] = vi->z;
            m_u(i) = vi->intensity;
        } else if (i < (nv + nf)) {
            vi = &getTriangleCenters()[i - nv];
            p0[0] = vi->x;
            p0[1] = vi->y;
            p0[2] = vi->z;
            m_u(i) = vi->intensity;
        }

        for (int j = 0; j < matrixDimension; j++) {
            double p1[] = { centerX, centerY, centerZ };

            if (j < nv) {
                vj = &m_mesh->getVertices()[j];
                p1[0] = vj->x;
                p1[1] = vj->y;
                p1[2] = vj->z;
            } else if (j < (nv + nf)) {
                vj = &getTriangleCenters()[j - nv];
                p1[0] = vj->x;
                p1[1] = vj->y;
                p1[2] = vj->z;
            }

            double r = m_computeDistance(p0, p1, tensor);
            m_distanceMatrix(i, j) = m_basis->phi(r);
        }
    }

    m_coeff = m_distanceMatrix.fullPivLu().solve(m_u); // call Eigen library to solve
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
