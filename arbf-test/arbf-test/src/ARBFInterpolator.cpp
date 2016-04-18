//
//  ARBFInterpolator.cpp
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#include <cassert>
#include <cstdio>
#include <iostream>
#include <cfenv> // for std::fesetround()
#include <unordered_set>
#include <Eigen/Eigenvalues>
#include "../include/Config.h"
#include "../include/ARBFInterpolator.h"

using namespace std;
using namespace Eigen;

Vertex* ARBFInterpolator::s_centers = nullptr;
Vertex* ARBFInterpolator::s_edgeCenters = nullptr;
const TriMesh* ARBFInterpolator::s_mesh = nullptr;
BasisType ARBFInterpolator::s_basis = BasisType::MQ;
vector<Matrix3d> ARBFInterpolator::s_T;
vector<vector<int>> ARBFInterpolator::s_vf_1ring;
vector<set<int>> ARBFInterpolator::s_vv_1ring;
vector<set<int>> ARBFInterpolator::s_ff_1ring;
vector<vector<int>> ARBFInterpolator::s_vf_2ring;
vector<set<int>> ARBFInterpolator::s_vv_2ring;
vector<set<int>> ARBFInterpolator::s_ff_2ring;

ARBFInterpolator::ARBFInterpolator() {}

// MARK: Getters/Setters
void ARBFInterpolator::setMesh(const TriMesh *mesh) {
    s_mesh = mesh;
}

inline vector<vector<int> >& ARBFInterpolator::getNeigborVertexFace1Ring() {
    return s_vf_1ring;
}

inline void ARBFInterpolator::setNeighborVertexFace1Ring(const vector<vector<int> > &vf) {
    s_vf_1ring = vf;
}

inline vector<set<int> >& ARBFInterpolator::getNeighborVertexVertex1Ring() {
    return s_vv_1ring;
}

inline void ARBFInterpolator::setNeighborVertexVertex1Ring(const vector<set<int> > &vv) {
    s_vv_1ring = vv;
}

inline vector<set<int> >& ARBFInterpolator::getNeighborFaceFace1Ring() {
    return s_ff_1ring;
}

inline void ARBFInterpolator::setNeighborFaceFace1Ring(vector<set<int> > &ff) {
    s_ff_1ring = ff;
}

inline vector<vector<int> >& ARBFInterpolator::getNeigborVertexFace2Ring() {
    return s_vf_2ring;
}

inline void ARBFInterpolator::setNeighborVertexFace2Ring(const vector<vector<int> > &vf) {
    s_vf_2ring = vf;
}

inline vector<set<int> >& ARBFInterpolator::getNeighborVertexVertex2Ring() {
    return s_vv_2ring;
}

inline void ARBFInterpolator::setNeighborVertexVertex2Ring(const vector<set<int> > &vv) {
    s_vv_2ring = vv;
}

inline vector<set<int> >& ARBFInterpolator::getNeighborFaceFace2Ring() {
    return s_ff_2ring;
}

inline void ARBFInterpolator::setNeighborFaceFace2Ring(vector<set<int> > &ff) {
    s_ff_2ring = ff;
}

// MARK: Shape Functions
inline double ARBFInterpolator::s_basisMQ(double r) {
    return std::sqrt(SQ(r) + SQ(0));
}

inline double ARBFInterpolator::s_basisIMQ(double r) {
    return 1.0 / (sqrt(SQ(r) + SQ(1.8)));
}

inline double ARBFInterpolator::s_basisGaussian(double r) {
    double c = 0.01;
    return exp(-SQ(c*r));
}

inline double ARBFInterpolator::s_basisTPS(double r) {
    if (abs(r) < Config::EPSILON)
        return 0;
    else
        return SQ(r) * log(r);
}

// MARK: ARBF Related
inline double ARBFInterpolator::s_computeC(double mu1, double mu2) {
    double sumMu = mu1 + mu2;
    
    if (sumMu > Config::EPSILON || sumMu < -Config::EPSILON) {
        return abs(mu1 - mu2) / sumMu;
    } else {
        return 1.0;
    }
}

inline double ARBFInterpolator::s_computeWeight(double r) {
    return exp(-SQ(r) / SQ(30.0));
}

inline double ARBFInterpolator::s_phi(double r) {
    switch (s_basis) {
        case MQ:
        default:
            return s_basisMQ(r);
        case IMQ:
            return s_basisIMQ(r);
        case Gaussian:
            return s_basisGaussian(r);
        case TPS:
            return s_basisTPS(r);
    }
}

tuple<double, double, vector<double>, vector<double>> ARBFInterpolator::s_computeEigen(const Eigen::Matrix2d &m) {
    assert(abs(m(0, 1) - m(1, 0)) < Config::EPSILON); // m should be symmetric
    
    // compute eigenvalues
    // matrix looks like
    //      | a    c|
    // A =  |       |
    //      | c    b|
    
    double a = m(0, 0), b = m(1, 1), c = m(0, 1);
    double lambda1 = 0.5 * (a + b) + sqrt(SQ(c) + 0.25 * SQ(a - b));
    double lambda2 = 0.5 * (a + b) - sqrt(SQ(c) + 0.25 * SQ(a - b));
    
    // compute eigenvectors
    double v1[2], v2[2];
    double u1 = 0.0, w1 = 0.0;
    
    if (abs(lambda1 - a) > Config::EPSILON) {
        // lambda1 != a
        u1 = c;
        w1 = lambda1 - a;
    } else if (abs(lambda1 - b) > Config::EPSILON) {
        // lambda1 != b
        u1 = lambda1 - b;
        w1 = c;
    } else {
        v1[0] = 1.0;
        v1[1] = 0.0;
        v2[0] = v1[1];
        v2[1] = -v1[0];
    }
    
    double p = sqrt(SQ(u1) + SQ(w1));
    v1[0] = u1 / p;
    v1[1] = w1 / p;
    v2[0] = v1[1];
    v2[1] = -v1[0];
    
    return make_tuple(lambda1, lambda2,
                      vector<double>(v1, v1 + sizeof(v1) / sizeof(double)),
                      vector<double>(v2, v2 + sizeof(v2) / sizeof(double)));
}

void ARBFInterpolator::computeMetric() {
    s_T.resize(s_mesh->getNumEdges());
    
    for (int e = 0; e < s_mesh->getNumEdges(); e++) {
        const Vertex& edgeCenter = s_edgeCenters[e];
        //             |lambda1 0      | |v1T|
        // T = [v1 v2] |               | |   |
        //             |0       lambda2| |v2T|
        
        Matrix3d t1, t2;
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
        
        s_T[e] = t1 * t2;
        s_T[e] = s_T[e] * t1.transpose();
    }
}

bool ARBFInterpolator::s_isInTriangle(const double* x, int triangleID) {
    const Vertex &a = s_mesh->getVertices()[s_mesh->getFaces()[triangleID].a];
    const Vertex &b = s_mesh->getVertices()[s_mesh->getFaces()[triangleID].b];
    const Vertex &c = s_mesh->getVertices()[s_mesh->getFaces()[triangleID].c];
    
    Eigen::Vector2d ab(b.x - a.x, b.y - a.y);
    Eigen::Vector2d bc(c.x - b.x, c.y - b.y);
    Eigen::Vector2d ca(a.x - c.x, a.y - c.y);
    Eigen::Vector2d ax(x[0] - a.x, x[1] - a.y);
    Eigen::Vector2d bx(x[0] - b.x, x[1] - b.y);
    Eigen::Vector2d cx(x[0] - c.x, x[1] - c.y);
    
    double area = 0.0, area1 = 0.0, area2 = 0.0, area3 = 0.0;
    double cosTheta = ab.dot(-ca) / (ab.norm() * ca.norm());
    area = 0.5 * ab.norm() * ca.norm() * sqrt(1.0 - SQ(cosTheta));
    double cosTheta1 = ab.dot(ax) / (ab.norm() * ax.norm());
    area1 = 0.5 * ab.norm() * ax.norm() * sqrt(1.0 - SQ(cosTheta1));
    double cosTheta2 = cx.dot(ca) / (cx.norm() * ca.norm());
    area2 = 0.5 * cx.norm() * ca.norm() * sqrt(1.0 - SQ(cosTheta2));
    double cosTheta3 = bx.dot(bc) / (bx.norm() * bc.norm());
    area3 = 0.5 * bx.norm() * bc.norm() * sqrt(1.0 - SQ(cosTheta3));
    
    if (isnan(area)) area = 0.0;
    if (isnan(area1)) area1 = 0.0;
    if (isnan(area2)) area2 = 0.0;
    if (isnan(area3)) area3 = 0.0;
    
    if (abs(area1 + area2 + area3 - area) < Config::EPSILON)
        return true;
    else
        return false;
}

inline double ARBFInterpolator::s_computeDistance(const double *v0, const double *v1, const Matrix3d& T) {
    RowVector3d dis(v0[0]-v1[0], v0[1]-v1[1], v0[2]-v1[2]);
    RowVector3d tmp = dis * T;
    return sqrt(tmp(0)*dis(0) + tmp(1)*dis(1) + tmp(2)*dis(2));
//    return sqrt(SQ(v1[0]-v0[0]) + SQ(v1[1]-v0[1]) + SQ(v1[2]-v0[2]));
}

// MARK: RBF Related
tuple<double*, double, int, int, int> ARBFInterpolator::interpolate() {
    const int NP = 250000; // number of evaluation points
    int dimX = 0, dimY = 0, dimZ = 0;
    
    if (Config::dim == 2) {
        dimX = (int) sqrt(NP);
        dimY = (int) sqrt(NP);
    } else {
        dimX = round(pow(NP, 1.0/3.0));
        dimY = round(pow(NP, 1.0/3.0));
        dimZ = round(pow(NP, 1.0/3.0));
    }
    
    vector<double> x, y, z;
    
    if (Config::dim == 2) {
        x = linspace(s_mesh->getMinX(), s_mesh->getMaxX(), sqrt(NP));
        y = linspace(s_mesh->getMinY(), s_mesh->getMaxY(), sqrt(NP));
    } else {
        x = linspace(s_mesh->getMinX(), s_mesh->getMaxX(), round(pow(NP, 1.0/3.0)));
        y = linspace(s_mesh->getMinY(), s_mesh->getMaxY(), round(pow(NP, 1.0/3.0)));
        z = linspace(s_mesh->getMinZ(), s_mesh->getMaxZ(), round(pow(NP, 1.0/3.0)));
    }
    
    const int NV = s_mesh->getNumVertices(), NF = s_mesh->getNumFaces(), NE = s_mesh->getNumEdges();
    int DIM = 0; // distance matrix dimension
    
    if (Config::dim == 2) {
        DIM = NV + NF + NE;
    } else {
        DIM = NV + NF + 1; // with additional tetrahedron center
    }
    
    //    cout << "\nx = ";
    //    copy(x.begin(), x.end(), ostream_iterator<double>(cout, ", "));
    //    cout << "\ny = ";
    //    copy(y.begin(), y.end(), ostream_iterator<double>(cout, ", "));
    //    cout << "\nz = ";
    //    copy(z.begin(), z.end(), ostream_iterator<double>(cout, ", "));
    
    MatrixXd ma(DIM, DIM); // distance matrix
    VectorXd u(DIM); // intensities
    VectorXd coeff(DIM); // weight
    double *intensities = new double[NP]; // solution
    memset(intensities, 0, sizeof(double) * NP);
    buildMatrix(DIM, NV, NF, ma, u); // assemble distance matrix and vector
    
    //    cout << "\nma = \n";
    //    cout << ma << "\n" << "-----" << endl;
    //    cout << ma - ma.transpose() << "-----" << endl;
    //    cout << u << "\n" << "-----" << endl;
    
    coeff = ma.fullPivLu().solve(u); // solve interpolation coefficents
                                     //    coeff = ma.inverse() * u;
                                     //    cout << "\ncoeff = \n" << coeff << "\n" << "-----" << endl;
                                     //    cout << "\ncoeff.sum() = " << coeff.sum() << endl;
    
    if (Config::dim == 2) {
        for (int j = 0; j < dimY; j++) {
            for (int i = 0; i < dimX; i++) {
                double p0[] = { x[i], y[j], 0.0 };
                double intensity = 0;
                for (int l = 0; l < DIM; l++) {
                    double p1[] = { -1, -1, -1 };
                    double r = 0.0;
                    
                    if (l < NV) {
                        p1[0] = s_mesh->getVertices()[l].x;
                        p1[1] = s_mesh->getVertices()[l].y;
                        p1[2] = s_mesh->getVertices()[l].z;
                        r = s_computeDistance(p0, p1, Matrix3d::Identity());
                    } else if (l < (NV + NF)) {
                        p1[0] = getCenters()[l-NV].x;
                        p1[1] = getCenters()[l-NV].y;
                        p1[2] = getCenters()[l-NV].z;
                        r = s_computeDistance(p0, p1, Matrix3d::Identity());
                    } else {
                        p1[0] = getEdgeCenters()[l-NV-NF].x;
                        p1[1] = getEdgeCenters()[l-NV-NF].y;
                        p1[2] = getEdgeCenters()[l-NV-NF].z;
                        r = s_computeDistance(p0, p1, s_T[l-NV-NF]);
                    }
                    
//                    r = s_computeDistance(p0, p1);
                    intensity += (s_phi(r) * coeff(l));
                }
                intensities[j*dimX+i] = intensity;
            }
        }
        
        cout << "\n----- Intensity at face centers -----\n";
        for (int j = 0; j < dimY; j++) {
            for (int i = 0; i < dimX; i++) {
                for (int v = 0; v < s_mesh->getNumVertices(); v++) {
                    if (abs(x[i] - s_mesh->getVertices()[v].x) < 1e-6 && abs(y[j] - s_mesh->getVertices()[v].y) < 1e-6) {
                        printf("\tAt vert[%d][%lf, %lf], intensity[%d, %d][%lf, %lf] = %lf\n",
                               v, s_mesh->getVertices()[v].x, s_mesh->getVertices()[v].y, i, j,
                               x[i], y[j], intensities[j*dimX+i]);
                        break;
                    }
                }
                
                for (int f = 0; f < s_mesh->getNumFaces(); f++) {
                    if (abs(x[i] - getCenters()[f].x) < 1e-6 && abs(y[j] - getCenters()[f].y) < 1e-6) {
                        printf("\tAt face[%d][%lf, %lf], intensity[%d, %d][%lf, %lf] = %lf\n",
                               f, getCenters()[f].x, getCenters()[f].y, i, j,
                               x[i], y[j], intensities[j*dimX+i]);
                        break;
                    }
                }
                
                for (int e = 0; e < s_mesh->getNumEdges(); e++) {
                    if (abs(x[i] - getEdgeCenters()[e].x) < 1e-6 && abs(y[j] - getEdgeCenters()[e].y) < 1e-6) {
                        printf("\tAt edgeCenter[%d][%lf, %lf], intensity[%d, %d][%lf, %lf] = %lf\n",
                               e, getEdgeCenters()[e].x, getEdgeCenters()[e].y, i, j,
                               x[i], y[j], intensities[j*dimX+i]);
                        break;
                    }
                }
            }
        }
        cout << "\n" << "-----" << endl;
    }
    // 3D
    else {
        for (int k = 0; k < dimZ; k++) {
            for (int j = 0; j < dimY; j++) {
                for (int i = 0; i < dimX; i++) {
                    double p0[] = { x[i], y[j], z[k] };
                    double intensity = 0;
                    for (int l = 0; l < DIM; l++) {
                        double p1[] = { -1, -1, -1 };
                        
                        if (l < NV) {
                            p1[0] = s_mesh->getVertices()[l].x;
                            p1[1] = s_mesh->getVertices()[l].y;
                            p1[2] = s_mesh->getVertices()[l].z;
                        } else if (l < (NV + NF)) {
                            p1[0] = getCenters()[l-NV].x;
                            p1[1] = getCenters()[l-NV].y;
                            p1[2] = getCenters()[l-NV].z;
                        } else {
                            p1[0] = (s_mesh->getMinX() + s_mesh->getMaxX()) / 2.0;
                            p1[1] = (s_mesh->getMinY() + s_mesh->getMaxY()) / 2.0;
                            p1[2] = (s_mesh->getMinZ() + s_mesh->getMaxZ()) / 2.0;
                        }
                        
                        double r = s_computeDistance(p0, p1, Matrix3d::Identity());
                        intensity += (s_phi(r) * coeff(l));
                    }
                    intensities[k*dimX*dimY + j*dimX + i] = intensity;
                }
            }
        }
        
        cout << "\n----- Intensity at face centers -----\n";
        for (int k = 0; k < dimZ; k++) {
            for (int j = 0; j < dimY; j++) {
                for (int i = 0; i < dimX; i++) {
                    for (int v = 0; v < s_mesh->getNumVertices(); v++) {
                        if (abs(x[i] - s_mesh->getVertices()[v].x) < 1e-6 &&
                            abs(y[j] - s_mesh->getVertices()[v].y) < 1e-6 &&
                            abs(z[k] - s_mesh->getVertices()[v].z) < 1e-6)
                        {
                            printf("\tAt vert[%d][%lf, %lf, %lf], intensity[%d, %d, %d][%lf, %lf, %lf] = %lf\n",
                                   v, s_mesh->getVertices()[v].x, s_mesh->getVertices()[v].y, s_mesh->getVertices()[v].z, i, j, k,
                                   x[i], y[j], z[k], intensities[k*dimY*dimX+j*dimX+i]);
                            break;
                        }
                    }
                    
                    for (int f = 0; f < s_mesh->getNumFaces(); f++) {
                        if (abs(x[i] - getCenters()[f].x) < 1e-6 &&
                            abs(y[j] - getCenters()[f].y) < 1e-6 &&
                            abs(z[k] - getCenters()[f].z) < 1e-6)
                        {
                            printf("\tAt face[%d][%lf, %lf, %lf], intensity[%d, %d, %d][%lf, %lf, %lf] = %lf\n",
                                   f, getCenters()[f].x, getCenters()[f].y, getCenters()[f].z, i, j, k,
                                   x[i], y[j], z[k], intensities[k*dimY*dimX+j*dimX+i]);
                            break;
                        }
                    }
                    
                    if (abs(x[i] - (s_mesh->getMinX()+s_mesh->getMaxX())/2.0) < 1e-6 &&
                        abs(y[j] - (s_mesh->getMinY()+s_mesh->getMaxY())/2.0) < 1e-6 &&
                        abs(z[k] - (s_mesh->getMinZ()+s_mesh->getMaxZ())/2.0) < 1e-6)
                    {
                        printf("\tAt center[%lf, %lf, %lf], intensity[%d, %d, %d][%lf, %lf, %lf] = %lf\n",
                               (s_mesh->getMinX()+s_mesh->getMaxX())/2.0,
                               (s_mesh->getMinY()+s_mesh->getMaxY())/2.0,
                               (s_mesh->getMinZ()+s_mesh->getMaxZ())/2.0, i, j, k,
                               x[i], y[j], z[k], intensities[k*dimY*dimX+j*dimX+i]);
                    }
                }
            }
        }
        cout << "\n" << "-----" << endl;
    }
    
    // rescale to [0, 255]
    for (int i = 0; i < NP; i++) {
        intensities[i] = (intensities[i] + 1) / 2 * 255;
    }
    
    const double MAX_INTENSITY = 255.0;
    
    // set threshold
    if (Config::dim == 2) {
        for (int j = 0; j < dimY; ++j) {
            for (int i = 0; i < dimX; ++i) {
                if (intensities[j*dimX + i] < 0.0) {
                    intensities[j*dimX + i] = 0.0;
                }
                
                if (intensities[j*dimX + i] > MAX_INTENSITY) {
                    intensities[j*dimX + i] = MAX_INTENSITY;
                }
            }
        }
    } else {
        for (int k = 0; k < dimZ; ++k) {
            for (int j = 0; j < dimY; ++j) {
                for (int i = 0; i < dimX; ++i) {
                    if (intensities[k*dimX*dimY + j*dimX + i] < 0.0) {
                        intensities[k*dimX*dimY + j*dimX + i] = 0.0;
                    }
                    
                    if (intensities[k*dimX*dimY + j*dimX + i] > MAX_INTENSITY) {
                        intensities[k*dimX*dimY + j*dimX + i] = MAX_INTENSITY;
                    }
                }
            }
        }
    }
    
    //    cout << "\nintensity after rescale = \n";
    //    for (int j = 0; j < dimY; j++) {
    //        for (int i = 0; i < dimX; i++) {
    //            cout << (int) round(intensities[j*dimX+i]) << " ";
    //        }
    //        cout << "\n";
    //    }
    //    cout << "\n" << "-----" << endl;
    
    return tuple<double*, double, int, int, int>(intensities, MAX_INTENSITY, dimX, dimY, dimZ);
}

void ARBFInterpolator::buildMatrix(int DIM, int NV, int NF, MatrixXd &ma, VectorXd &u) {
    if (Config::dim == 2) { // 2D
        for (int i = 0; i < DIM; i++) {
            if (i < NV) {
                u(i) = s_mesh->getVertices()[i].intensity;
                double p0[] = { s_mesh->getVertices()[i].x, s_mesh->getVertices()[i].y, s_mesh->getVertices()[i].z };
                
                for (int j = 0; j < DIM; j++) {
                    double p1[] = { -1, -1, -1 };
                    double r = 0.0;
                    
                    if (j < NV) {
                        p1[0] = s_mesh->getVertices()[j].x;
                        p1[1] = s_mesh->getVertices()[j].y;
                        p1[2] = s_mesh->getVertices()[j].z;
                        r = s_computeDistance(p0, p1, Matrix3d::Identity());
                    } else if (j < (NV + NF)) {
                        p1[0] = getCenters()[j-NV].x;
                        p1[1] = getCenters()[j-NV].y;
                        p1[2] = getCenters()[j-NV].z;
                        r = s_computeDistance(p0, p1, Matrix3d::Identity());
                    } else {
                        p1[0] = getEdgeCenters()[j-NV-NF].x;
                        p1[1] = getEdgeCenters()[j-NV-NF].y;
                        p1[2] = getEdgeCenters()[j-NV-NF].z;
                        r = s_computeDistance(p0, p1, s_T[j-NV-NF]);
                    }
                    
//                    r = s_computeDistance(p0, p1);
                    ma(i, j) = s_phi(r);
                }
            } else if (i < (NV + NF)) {
                u(i) = getCenters()[i-NV].intensity;
                double p0[] = { getCenters()[i-NV].x, getCenters()[i-NV].y, getCenters()[i-NV].z };
                
                for (int j = 0; j < DIM; j++) {
                    double p1[] = { -1, -1, -1 };
                    double r = 0.0;
                    
                    if (j < NV) {
                        p1[0] = s_mesh->getVertices()[j].x;
                        p1[1] = s_mesh->getVertices()[j].y;
                        p1[2] = s_mesh->getVertices()[j].z;
                        r = s_computeDistance(p0, p1, Matrix3d::Identity());
                    } else if (j < (NV + NF)) {
                        p1[0] = getCenters()[j-NV].x;
                        p1[1] = getCenters()[j-NV].y;
                        p1[2] = getCenters()[j-NV].z;
                        r = s_computeDistance(p0, p1, Matrix3d::Identity());
                    } else {
                        p1[0] = getEdgeCenters()[j-NV-NF].x;
                        p1[1] = getEdgeCenters()[j-NV-NF].y;
                        p1[2] = getEdgeCenters()[j-NV-NF].z;
                        r = s_computeDistance(p0, p1, s_T[j-NV-NF]);
                    }
                    
//                    r = s_computeDistance(p0, p1);
                    ma(i, j) = s_phi(r);
                }
            } else {
                u(i) = getEdgeCenters()[i-NV-NF].intensity;
                double p0[] = { getEdgeCenters()[i-NV-NF].x, getEdgeCenters()[i-NV-NF].y, getEdgeCenters()[i-NV-NF].z };
                
                for (int j = 0; j < DIM; j++) {
                    double p1[] = { -1, -1, -1 };
                    double r = 0.0;
                    
                    if (j < NV) {
                        p1[0] = s_mesh->getVertices()[j].x;
                        p1[1] = s_mesh->getVertices()[j].y;
                        p1[2] = s_mesh->getVertices()[j].z;
                        r = s_computeDistance(p0, p1, Matrix3d::Identity());
                    } else if (j < (NV + NF)) {
                        p1[0] = getCenters()[j-NV].x;
                        p1[1] = getCenters()[j-NV].y;
                        p1[2] = getCenters()[j-NV].z;
                        r = s_computeDistance(p0, p1, Matrix3d::Identity());
                    } else {
                        p1[0] = getEdgeCenters()[j-NV-NF].x;
                        p1[1] = getEdgeCenters()[j-NV-NF].y;
                        p1[2] = getEdgeCenters()[j-NV-NF].z;
                        r = s_computeDistance(p0, p1, s_T[j-NV-NF]);
                    }
                    
//                    r = s_computeDistance(p0, p1);
                    ma(i, j) = s_phi(r);
                }
            }
        }
    }
    // 3D
    else {
        // tetrahedron center
        double centerX = (s_mesh->getMinX() + s_mesh->getMaxX()) / 2.0;
        double centerY = (s_mesh->getMinY() + s_mesh->getMaxY()) / 2.0;
        double centerZ = (s_mesh->getMinZ() + s_mesh->getMaxZ()) / 2.0;
        
        for (int i = 0; i < DIM; i++) {
            if (i < NV) {
                u(i) = s_mesh->getVertices()[i].intensity;
                double p0[] = { s_mesh->getVertices()[i].x, s_mesh->getVertices()[i].y, s_mesh->getVertices()[i].z };
                
                for (int j = 0; j < DIM; j++) {
                    double p1[] = { -1, -1, -1 };
                    
                    if (j < NV) {
                        p1[0] = s_mesh->getVertices()[j].x;
                        p1[1] = s_mesh->getVertices()[j].y;
                        p1[2] = s_mesh->getVertices()[j].z;
                    } else if (j < (NV + NF)) {
                        p1[0] = getCenters()[j-NV].x;
                        p1[1] = getCenters()[j-NV].y;
                        p1[2] = getCenters()[j-NV].z;
                    } else {
                        p1[0] = centerX;
                        p1[1] = centerY;
                        p1[2] = centerZ;
                    }
                    
                    double r = s_computeDistance(p0, p1, Matrix3d::Identity());
                    ma(i, j) = s_phi(r);
                }
            } else if (i < (NV + NF)) {
                u(i) = getCenters()[i-NV].intensity;
                double p0[] = { getCenters()[i-NV].x, getCenters()[i-NV].y, getCenters()[i-NV].z };
                
                for (int j = 0; j < DIM; j++) {
                    double p1[] = { -1, -1, -1 };
                    
                    if (j < NV) {
                        p1[0] = s_mesh->getVertices()[j].x;
                        p1[1] = s_mesh->getVertices()[j].y;
                        p1[2] = s_mesh->getVertices()[j].z;
                    } else if (j < (NV + NF)) {
                        p1[0] = getCenters()[j-NV].x;
                        p1[1] = getCenters()[j-NV].y;
                        p1[2] = getCenters()[j-NV].z;
                    } else {
                        p1[0] = centerX;
                        p1[1] = centerY;
                        p1[2] = centerZ;
                    }
                    
                    double r = s_computeDistance(p0, p1, Matrix3d::Identity());
                    ma(i, j) = s_phi(r);
                }
            } else {
                u(i) = -1;
                double p0[] = { centerX, centerY, centerZ };
                
                for (int j = 0; j < DIM; j++) {
                    double p1[] = { -1, -1, -1 };
                    
                    if (j < NV) {
                        p1[0] = s_mesh->getVertices()[j].x;
                        p1[1] = s_mesh->getVertices()[j].y;
                        p1[2] = s_mesh->getVertices()[j].z;
                    } else if (j < (NV + NF)) {
                        p1[0] = getCenters()[j-NV].x;
                        p1[1] = getCenters()[j-NV].y;
                        p1[2] = getCenters()[j-NV].z;
                    } else {
                        p1[0] = centerX;
                        p1[1] = centerY;
                        p1[2] = centerZ;
                    }
                    
                    double r = s_computeDistance(p0, p1, Matrix3d::Identity());
                    ma(i, j) = s_phi(r);
                }
            }
        }
    }
}

// MARK: Utility Functions
void ARBFInterpolator::setBasisType(const BasisType type) {
    s_basis = type;
}

double ARBFInterpolator::computePSNR(const PPMImage &img1, const PPMImage &img2) {
    return 20.0 * log10(255.0 / s_computeRMSE(img1, img2));
}

inline double ARBFInterpolator::s_computeRMSE(const PPMImage &img1, const PPMImage &img2) {
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

inline vector<double> ARBFInterpolator::linspace(double a, double b, int n) {
    assert(b > a && n > 0);
    vector<double> res;
    double step = (b-a) / (n-1);
    
    while (a < b) {
        res.push_back(a);
        a += step;
    }
    
    // get last point
    if ((a-b) < Config::EPSILON) {
        res.push_back(b);
    }
    
    return res;
}

void ARBFInterpolator::clean() {
    delete [] s_edgeCenters;
    s_edgeCenters = nullptr;
    
    delete [] s_centers;
    s_centers = nullptr;
}

// MARK: Find Neighboring Triangles and Vertices
void ARBFInterpolator::findNeighVertices1Ring() {
    s_vv_1ring.resize(s_mesh->getNumVertices());
    s_vf_1ring.resize(s_mesh->getNumVertices());
    
    // build 1-ring VERTEX <-> FACE mapping
    for (int f = 0; f < s_mesh->getNumFaces(); ++f) {
        s_vf_1ring[s_mesh->getFaces()[f].a].push_back(f);
        s_vf_1ring[s_mesh->getFaces()[f].b].push_back(f);
        s_vf_1ring[s_mesh->getFaces()[f].c].push_back(f);
    }
    
    // find 1-ring VERTEX <-> VERTEX mapping
    for (int v = 0; v < s_mesh->getNumVertices(); ++v) {
        for (size_t j = 0; j < s_vf_1ring[v].size(); ++j) {
            s_vv_1ring[v].insert(s_mesh->getFaces()[s_vf_1ring[v][j]].a);
            s_vv_1ring[v].insert(s_mesh->getFaces()[s_vf_1ring[v][j]].b);
            s_vv_1ring[v].insert(s_mesh->getFaces()[s_vf_1ring[v][j]].c);
        }
    }
    
    if (Config::isPrintingDebugInfo) {
        string filename = string("DEBUG.vf_1ring.txt");
        FILE *fp = nullptr;
        
        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
            fprintf(stderr, "WARNING: open DEBUG.vf_1ring.txt failed ...\n");
            return;
        }
        
        for (int i = 0; i < s_vf_1ring.size(); ++i) {
            fprintf(fp, "[v <-> f] [%ld] v%d: ", s_vf_1ring[i].size(), i);
            for (int j = 0; j < s_vf_1ring[i].size(); ++j) {
                fprintf(fp, "%d ", s_vf_1ring[i][j]);
            }
            fprintf(fp, "\n");
        }
        
        fclose(fp);
        fp = nullptr;
        filename = string("DEBUG.vv_1ring.txt");
        
        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
            fprintf(stderr, "WARNING: open DEBUG.vv_1ring.txt failed ...\n");
            return;
        }
        
        for (int i = 0; i < s_vv_1ring.size(); ++i) {
            fprintf(fp, "[v <-> v] [%ld] v%d: ", s_vv_1ring[i].size(), i);
            for (set<int>::const_iterator it = s_vv_1ring[i].begin();
                 it != s_vv_1ring[i].end(); ++it) {
                fprintf(fp, "%d ", *it);
            }
            fprintf(fp, "\n");
        }
        
        fclose(fp);
    }
}

void ARBFInterpolator::findNeighFaces1Ring() {
    s_ff_1ring.resize(s_mesh->getNumFaces());
    
    for (int f = 0; f < s_mesh->getNumFaces(); ++f) {
        int a = s_mesh->getFaces()[f].a;
        int b = s_mesh->getFaces()[f].b;
        int c = s_mesh->getFaces()[f].c;
        
        for (size_t j = 0; j < s_vf_1ring[a].size(); ++j) {
            s_ff_1ring[f].insert(s_vf_1ring[a][j]);
        }
        
        for (size_t j = 0; j < s_vf_1ring[b].size(); ++j) {
            s_ff_1ring[f].insert(s_vf_1ring[b][j]);
        }
        
        for (size_t j = 0; j < s_vf_1ring[c].size(); ++j) {
            s_ff_1ring[f].insert(s_vf_1ring[c][j]);
        }
    }
    
    if (Config::isPrintingDebugInfo) {
        string filename = string("DEBUG.ff_1ring.txt");
        FILE *fp = nullptr;
        
        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
            fprintf(stderr, "WARNING: open DEBUG.ff_1ring.txt failed ...\n");
            return;
        }
        
        for (int f = 0; f < s_ff_1ring.size(); ++f) {
            fprintf(fp, "[f <-> f] [%ld] f%d: ", s_ff_1ring[f].size(), f);
            for (set<int>::const_iterator it = s_ff_1ring[f].begin();
                 it != s_ff_1ring[f].end(); ++it) {
                fprintf(fp, "%d ", *it);
            }
            fprintf(fp, "\n");
        }
        
        fclose(fp);
    }
}

void ARBFInterpolator::findNeighVertices2Ring() {
    s_vv_2ring.resize(s_mesh->getNumVertices());
    s_vf_2ring.resize(s_mesh->getNumVertices());
    
    // build 2-ring VERTEX <-> FACE mapping
    s_vf_2ring.assign(s_vf_1ring.begin(), s_vf_1ring.end());
    for (int i = 0; i < s_mesh->getNumFaces(); ++i) {
        // 3 vertices of current face
        int a = s_mesh->getFaces()[i].a;
        int b = s_mesh->getFaces()[i].b;
        int c = s_mesh->getFaces()[i].c;
        
        // neighboring vertices of 3 vertices
        const std::set<int> &neib1 = s_vv_1ring[a];
        const std::set<int> &neib2 = s_vv_1ring[b];
        const std::set<int> &neib3 = s_vv_1ring[c];
        
        // go over each neighboring vertices set, check if current face is in
        // the neighboring faces list
        std::set<int>::const_iterator it;
        
        for (it = neib1.begin(); it != neib1.end(); ++it) {
            std::vector<int> &neib_faces = s_vf_1ring[*it];
            if (find(neib_faces.begin(), neib_faces.end(), i) == neib_faces.end()) {
                // current face is not in the list, add it
                s_vf_2ring[*it].push_back(i);
            }
        }
        
        for (it = neib2.begin(); it != neib2.end(); ++it) {
            std::vector<int> &neib_faces = s_vf_1ring[*it];
            if (find(neib_faces.begin(), neib_faces.end(), i) == neib_faces.end()) {
                // current face is not in the list, add it
                s_vf_2ring[*it].push_back(i);
            }
        }
        
        for (it = neib3.begin(); it != neib3.end(); ++it) {
            std::vector<int> &neib_faces = s_vf_1ring[*it];
            if (find(neib_faces.begin(), neib_faces.end(), i) == neib_faces.end()) {
                // current face is not in the list, add it
                s_vf_2ring[*it].push_back(i);
            }
        }
    }
    
    // find 2-ring neighboring vertices
    s_vv_2ring.assign(s_vv_1ring.begin(), s_vv_1ring.end());
    for (int i = 0; i < s_mesh->getNumVertices(); ++i) {
        for (int j = 0; j < s_vf_2ring[i].size(); ++j) {
            s_vv_2ring[i].insert(s_mesh->getFaces()[s_vf_2ring[i][j]].a);
            s_vv_2ring[i].insert(s_mesh->getFaces()[s_vf_2ring[i][j]].b);
            s_vv_2ring[i].insert(s_mesh->getFaces()[s_vf_2ring[i][j]].c);
        }
    }
    
    if (Config::isPrintingDebugInfo) {
        string filename = "DEBUG.vf_2ring.txt";
        FILE *fp = nullptr;
        
        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
            fprintf(stderr, "WARNING: open DEBUG.vf_2ring.txt failed.\n");
            return;
        }
        
        for (int i = 0; i < s_vf_2ring.size(); ++i) {
            fprintf(fp, "[v <-> f] [%ld] v%d: ", s_vf_2ring[i].size(), i);
            for (int j = 0; j < s_vf_2ring[i].size(); ++j) {
                fprintf(fp, "%d ", s_vf_2ring[i][j]);
            }
            fprintf(fp, "\n");
        }
        
        fclose(fp);
        fp = nullptr;
        filename = "DEBUG.vv_2ring.txt";
        
        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
            fprintf(stderr, "WARNING: open DEBUG.vv_2ring.txt failed ...");
            return;
        }
        
        for (int i = 0; i < s_vv_2ring.size(); ++i) {
            fprintf(fp, "[v <-> v] [%ld] v%d: ", s_vv_2ring[i].size(), i);
            for (set<int>::const_iterator it = s_vv_2ring[i].begin();
                 it != s_vv_2ring[i].end(); ++it) {
                fprintf(fp, "%d ", *it);
            }
            fprintf(fp, "\n");
        }
        
        fclose(fp);
    }
}

void ARBFInterpolator::findNeighFaces2Ring() {
    s_ff_2ring.resize(s_mesh->getNumFaces());
    
    for (int f = 0; f < s_mesh->getNumFaces(); ++f) {
        int a = s_mesh->getFaces()[f].a;
        int b = s_mesh->getFaces()[f].b;
        int c = s_mesh->getFaces()[f].c;
        
        for (int j = 0; j < s_vf_2ring[a].size(); ++j) {
            s_ff_2ring[f].insert(s_vf_2ring[a][j]);
        }
        
        for (int j = 0; j < s_vf_2ring[b].size(); ++j) {
            s_ff_2ring[f].insert(s_vf_2ring[b][j]);
        }
        
        for (int j = 0; j < s_vf_2ring[c].size(); ++j) {
            s_ff_2ring[f].insert(s_vf_2ring[c][j]);
        }
    }
    
    if (Config::isPrintingDebugInfo) {
        string filename = "DEBUG.ff_2ring.txt";
        FILE *fp = nullptr;
        
        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
            fprintf(stderr, "WARNINIG: open DEBUG.ff_2ring.txt failed ...");
            return;
        }
        
        for (int f = 0; f < s_ff_2ring.size(); ++f) {
            fprintf(fp, "[f <-> f] [%ld] f%d: ", s_ff_2ring[f].size(), f);
            for (set<int>::const_iterator it = s_ff_2ring[f].begin();
                 it != s_ff_2ring[f].end(); ++it) {
                fprintf(fp, "%d ", *it);
            }
            fprintf(fp, "\n");
        }
        
        fclose(fp);
    }
}

// MARK: Calculate Triangle Centers and Edge Centers
void ARBFInterpolator::calculateEdgeCenters() {
    s_edgeCenters = new Vertex[s_mesh->getNumEdges()];
    
    if (s_edgeCenters == nullptr) {
        fprintf(stderr, "ERROR: allocate memory for triangle edge centers failed.\n");
        exit(EXIT_FAILURE);
    }
    
    // function to compute C, used to transform eigenvalues
    auto computeC = [](double mu1, double mu2) -> double {
        double sum = mu1 + mu2;
        
        if (sum > Config::EPSILON || sum < -Config::EPSILON) {
            return abs(mu1 - mu2) / sum;
        } else {
            return 1.0;
        }
    };
    
    int e = 0;
    for (unordered_set<Edge, edgeComparator>::const_iterator it = s_mesh->getEdges().begin(); it != s_mesh->getEdges().end(); it++, e++) {
        int a = it->a;
        int b = it->b;
        s_edgeCenters[e].x = s_mesh->getVertices()[a].x + 0.5 * (s_mesh->getVertices()[b].x - s_mesh->getVertices()[a].x);
        s_edgeCenters[e].y = s_mesh->getVertices()[a].y + 0.5 * (s_mesh->getVertices()[b].y - s_mesh->getVertices()[a].y);
        s_edgeCenters[e].z = s_mesh->getVertices()[a].z + 0.5 * (s_mesh->getVertices()[b].z - s_mesh->getVertices()[a].z);
        s_edgeCenters[e].intensity = it->intensity;
        
        // only calculate eigenvalues and eigenvectors for edges of 1st triangle
//        double dx = s_centers[0].x - s_edgeCenters[e].x;
//        double dy = s_centers[0].y - s_edgeCenters[e].y;
        double dx = s_mesh->getVertices()[b].x - s_mesh->getVertices()[a].x;
        double dy = s_mesh->getVertices()[b].y - s_mesh->getVertices()[a].y;
        double normDxy = sqrt(SQ(dx) + SQ(dy));
//        double gradX = dx / normDxy;
//        double gradY = dy / normDxy;
        double gradX = dx / normDxy;
        double gradY = dy / normDxy;
//        double gradZ = s_centers[0].z - s_edgeCenters[e].z;
//        Matrix3d mEdge;
        Matrix2d mEdge;
        mEdge(0, 0) = gradX * gradX;
        mEdge(0, 1) = gradX * gradY;
//        mEdge(0, 2) = gradX * gradZ;
        mEdge(1, 0) = gradY * gradX;
        mEdge(1, 1) = gradY * gradY;
//        mEdge(1, 2) = gradY * gradZ;
//        mEdge(2, 0) = gradZ * gradX;
//        mEdge(2, 1) = gradZ * gradY;
//        mEdge(2, 2) = gradZ * gradZ;
        
        // compute eigenvalues and eigenvectors
//        EigenSolver<Matrix3d> solver(mEdge);
        auto result = s_computeEigen(mEdge);
        double eig2 = 0.0;
        double c = computeC(get<0>(result), get<1>(result));
        
        if (c > Config::EPSILON || c < -Config::EPSILON) {
            eig2 = 1.0 - exp(-3.8 * SQ(Config::C0 / c));
        } else {
            eig2 = 1.0;
        }
        
        s_edgeCenters[e].eig1 = 1.0;
        s_edgeCenters[e].eig2 = eig2;
        s_edgeCenters[e].eig3 = 0.0;
        s_edgeCenters[e].eigVec1 = get<2>(result);
        s_edgeCenters[e].eigVec1.push_back(0);
        s_edgeCenters[e].eigVec2 = get<3>(result);
        s_edgeCenters[e].eigVec2.push_back(0);
        s_edgeCenters[e].eigVec3 = {0, 0, 1};
    }
    
    //    PPMImage im(301, 301);
    //    im.setMaxIntensity(255);
    //    im.setPath("/Users/keliu/Documents/projects/arbf/arbf-test/arbf-test/edge.ppm");
    //    double* data = new double[301*301];
    //    memset(data, 0, sizeof(double) * 301 * 301);
    //    for (int j = 0; j < 301; j++) {
    //        for (int i = 0; i < 301; i++) {
    //            for (int e = 0;  e < s_mesh->getNumEdges(); e++) {
    //                if (s_edgeCenters[e].x * 100 == i && s_edgeCenters[e].y * 100 == j) {
    //                    data[j*301+i] = 255;
    //                }
    //            }
    //        }
    //    }
    //    im.setImageData(data);
    //    im.write();
    
    if (Config::isPrintingDebugInfo) {
        string filename = "DEBUG.edge_centers.txt";
        FILE *fp = nullptr;
        
        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
            fprintf(stderr, "WARNING: open DEBUG.edge_centers.txt failed.\n");
            return;
        }
        
        fprintf(fp, "X Y INTENSITY\n");
        
        for (int e = 0; e < s_mesh->getNumEdges(); ++e) {
            fprintf(fp, "%lf %lf %lf %lf\n", s_edgeCenters[e].x, s_edgeCenters[e].y, s_edgeCenters[e].z, s_edgeCenters[e].intensity);
        }
        
        fclose(fp);
    }
}

void ARBFInterpolator::calculateCenters() {
    s_centers = new Vertex[s_mesh->getNumFaces()];
    
    if (s_centers == nullptr) {
        fprintf(stderr, "ERROR: allocate memory for triangle centers failed.\n");
        exit(EXIT_FAILURE);
    }
    
    for (int f = 0; f < s_mesh->getNumFaces(); ++f) {
        int a = s_mesh->getFaces()[f].a;
        int b = s_mesh->getFaces()[f].b;
        int c = s_mesh->getFaces()[f].c;
        
        // compute coordinates
        s_centers[f].x = (s_mesh->getVertices()[a].x + s_mesh->getVertices()[b].x + s_mesh->getVertices()[c].x) / 3.0;
        s_centers[f].y = (s_mesh->getVertices()[a].y + s_mesh->getVertices()[b].y + s_mesh->getVertices()[c].y) / 3.0;
        s_centers[f].z = (s_mesh->getVertices()[a].z + s_mesh->getVertices()[b].z + s_mesh->getVertices()[c].z) / 3.0;
        
        // compute intensity
        s_centers[f].intensity = s_mesh->getFaces()[f].intensity;
        
        // assign eigenvalues and eigenvectors
        s_centers[f].eig1 = 0;
        s_centers[f].eig2 = 0;
        s_centers[f].eig3 = 0;
        s_centers[f].eigVec1 = {1, 0, 0};
        s_centers[f].eigVec2 = {0, 1, 0};
        s_centers[f].eigVec3 = {0, 0, 1};
    }
    
    if (Config::isPrintingDebugInfo) {
        string filename = "DEBUG.centers.txt";
        FILE *fp = nullptr;
        
        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
            fprintf(stderr, "WARNING: open DEBUG.centers.txt failed.\n");
            return;
        }
        
        fprintf(fp, "X Y Z INTENSITY\n");
        
        for (int f = 0; f < s_mesh->getNumFaces(); ++f) {
            fprintf(fp, "%lf %lf %lf %lf\n", s_centers[f].x, s_centers[f].y, s_centers[f].z, s_centers[f].intensity);
        }
        
        fclose(fp);
    }
}