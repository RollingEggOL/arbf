//
//  ARBFInterpolator.h
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#ifndef __arbf_test_ARBFInterpolator__
#define __arbf_test_ARBFInterpolator__

#include <memory>
#include <vector>
#include <set>
#include <tuple>
#include <array>
#include <Eigen/Dense>
#include "Common.h"
#include "PPMImage.h"
#include "Mesh.h"
#include "BasisFunction.h"

class ARBFInterpolator {
public:
    typedef std::array<double, 2> Eigenvalues2d;
    typedef std::array<std::array<double, 2>, 2> Eigenvectors2d;
    typedef std::tuple<Eigenvalues2d, Eigenvectors2d> Eigens2d; // contains eigenvalue and eigenvectors

    typedef std::array<double, 3> Eigenvalues3d;
    typedef std::array<std::array<double, 3>, 3> Eigenvectors3d;
    typedef std::tuple<Eigenvalues3d, Eigenvectors3d> Eigens3d; // contains eigenvalue and eigenvectors

    // contains result data, maximum intensity of result, dimension of result data
    typedef std::tuple<std::vector<double>, double, std::array<unsigned, 3>> InterpolateResult;

public:
    ARBFInterpolator();
    virtual ~ARBFInterpolator();

    InterpolateResult getResult() const;
    void setMesh(Mesh *mesh);
    void setBasisFunction(BasisFunction *basisFunction);
    void computeEdgeMetrics();
    void computeTetrahedronMetrics();
    void interpolate_global(unsigned numEvalPoints);
    void interpolate_local(unsigned numEvalPoints);
    double computeDistance(const double *v0, const double *v1); // Euclidean distance

private:
    Eigens2d computeEigens(const Eigen::Matrix2d &matrix);
    Eigens3d computeEigens(const Eigen::Matrix3d &matrix);
    void rescaleEigenvalues(Eigenvalues2d &values); // for 2D
    void rescaleEigenvalues(Eigenvalues3d &values); // for 3D
    InterpolateResult m_interpolate3dGlobal(unsigned numEvalPoints);
    InterpolateResult m_interpolateHexahedron3dGlobal(unsigned numEvalPoints);
    InterpolateResult m_interpolate2dGlobal(unsigned numEvalPoints);
    InterpolateResult m_interpolate3d_local(unsigned numEvalPoints);
    InterpolateResult m_interpolate2d_local(unsigned numEvalPoints);
    bool m_isInTriangle(const double *x, int ta, int tb, int tc);
    bool m_isInTetrahedron(const double *x, int tetrahedronId);
    bool m_isInHexahedron(const double *x, int hexahedronId);
    bool m_isInAnyTetrahedron(const double *x); // check if a point x is in any tetrahedron
    bool m_isInAnyHexahedron(const double *x);
    double m_computeDistance(const double *v0, const double *v1, const Eigen::Matrix3d& T); // distance of tensor
    double m_computeDistance(const double *x, const double *v0, const double *v1); // distance between point x and line segment [v0, v1]
    double m_computeDistance(const double *v0, const double *v1, const double *w0, const double *w1); // distance of line segment [v0, v1] and [w0, w1]
    std::vector<double> linspace(double a, double b, int n); // generate n points between a and b inclusively
    void rescaleResultData(int low, int high); // rescale result data to range [low, high]
    void thresholdResultData(double maxIntensity); // threshold result data
    double computePSNR(const PPMImage &img1, const PPMImage &img2);
//    Vertex createFalseCenter(const Vertex &trueCenter, const Vertex &vertex, double offsetRatio);

private:
    double m_computeWeight(double r);
    double m_computeRMSE(const PPMImage &img1, const PPMImage &img2);

private:
    Mesh *m_mesh = nullptr; // Mesh
    BasisFunction *m_basis; // basis function
    std::vector<Eigen::Matrix3d> m_edgeT; // edge metric tensor
    std::vector<Eigen::Matrix3d> m_tetrahedronT; // tetrahedron metric tensor
    Eigen::MatrixXd m_distanceMatrix; // A in Ax = u
    Eigen::VectorXd m_u; // u in Ax = u
    Eigen::VectorXd m_coeff; // x in Ax = u
    InterpolateResult m_result;
};

#endif /* defined(__arbf_test_ARBFInterpolator__) */
