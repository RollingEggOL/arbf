//
//  ARBFInterpolator.h
//  arbf-image-super_resolution
//
//  Created by Ke Liu on 2/27/15.
//  Copyright (c) 2015 Ke Liu. All rights reserved.
//

#ifndef __arbf_image_super_resolution__ARBFInterpolator__
#define __arbf_image_super_resolution__ARBFInterpolator__

#include <vector>
#include <set>
#include <Eigen/Dense>
#include "Common.h"
#include "PPMImage.h"
#include "TriMesh.h"

class ARBFInterpolator
{
public:
    static const DBL7VECT* getCenters();
    
    static void setSRImage(PPMImage *img);
    
    static void setOriginalImage(const PPMImage *img);
    
    static void setMesh(const TriMesh *mesh);
    
    static void setBasisType(const BasisType type = BasisType::MQ);
    
    static std::vector<std::vector<int> >& getNeigborVertexFace1Ring();
    static void setNeighborVertexFace1Ring(const std::vector<std::vector<int> >& vf);
    
    static std::vector<std::set<int> >& getNeighborVertexVertex1Ring();
    static void setNeighborVertexVertex1Ring(const std::vector<std::set<int> >& vv);
    
    static std::vector<std::set<int> >& getNeighborFaceFace1Ring();
    static void setNeighborFaceFace1Ring(std::vector<std::set<int> >& ff);
    
    static std::vector<std::vector<int> >& getNeigborVertexFace2Ring();
    static void setNeighborVertexFace2Ring(const std::vector<std::vector<int> >& vf);
    
    static std::vector<std::set<int> >& getNeighborVertexVertex2Ring();
    static void setNeighborVertexVertex2Ring(const std::vector<std::set<int> >& vv);
    
    static std::vector<std::set<int> >& getNeighborFaceFace2Ring();
    static void setNeighborFaceFace2Ring(std::vector<std::set<int> >& ff);

    /*
     * Calculate coordinates, intensity, eigenvalue and eigenvectors for triangle centers.
     * \return: 1 for success, 0 for failure.
     */
    static int calculateCenters();
    
    /*
     * Find 1-ring VERTEX <-> VERTEX, VERTEX <-> FACE.
     * \return: 1 for success, 0 for failure.
     */
    static int findNeighVertices1Ring();
    
    /*
     * Find 1-ring FACE <-> FACE.
     * \return: 1 for success, 0 for failure.
     */
    static int findNeighFaces1Ring();
    
    /*
     * Find 2-ring VERTEX <-> VERTEX, VERTEX <-> FACE.
     * \return: 1 for success, 0 for failure.
     */
    static int findNeighVertices2Ring();
    
    /*
     * Find 2-ring FACE <-> FACE.
     * \return: 1 for success, 0 for failure.
     */
    static int findNeighFaces2Ring();

    
    /*
     * Rescale eigenvalues for every face.
     * \return: 1 for success, 0 for failure.
     */
    static int transformEigenvalues();
    
    /*
     * Compute anisotropic metrics for every face.
     */
    static void computeMetrics();
    
    /*
     * ARBF interpolation.
     * \return: 1 for success, 0 for failure.
     */
    static int interpolate();
    
    /*
     * ARBF weighted-interpolation (piecewise interpolation).
     * \return: 1 for success, 0 for failure.
     */
    static int weightedInterpolate();
    
    /*
     * Compute Peak-Signal-Noise-Ratio between original image and reconstructed image.
     * \return: PSNR.
     */
    static double computePSNR();
    
    /*
     * Release resources.
     */
    static void clean();
    
private:
    ARBFInterpolator(); // declare as private to prevent creating instances
    
    /*
     * Compute eigenvalues and corresponding eigenvectors.
     * \param m: input matrix.
     * \param lambda1: (OUTPUT) eigenvalue 1.
     * \param lambda2: (OUTPUT) eigenvalue 2.
     * \return: corresponding eigenvectors.
     */
    static std::pair<DBL2VECT, DBL2VECT> s_computeEigen(const Eigen::Matrix2d &m, double *lambda1, double *lambda2);
    
    /*
     * Compute C which is used on rescaling eigenvalues.
     * \param mu1: input parameter 1.
     * \param mu2: input parameter 2.
     * \return: C.
     */
    static double s_computeC(double mu1, double mu2);
    
    /*
     * Compute weight for weighted interpolation.
     * \param mu1: input distance.
     * \return: weight.
     */
    static double s_computeWeight(double r);
    
    /*
     * Compute anisotropic distance.
     * \param v0: input point 1.
     * \param mu2: input point 2.
     * \param T: anisotropic metric.
     * \return: distance betwee \param v0 and \param v1.
     */
    static double s_computeDistance(const double *v0, const double *v1,
                                    const Eigen::Matrix2d &T);
    
    /*
     * Compute chosen basis function based on s_basis.
     * \param r: input distance.
     * \return: value of radial basis function chosen.
     */
    static double s_phi(double r);
    
    /*
     * MQ radial basis function.
     * \param r: input distance.
     * \return: value of radial basis function.
     */
	static double s_basisMQ(double r);
    
    /*
     * IMQ radial basis function.
     * \param r: input distance.
     * \return: value of radial basis function.
     */
	static double s_basisIMQ(double r);
    
    /*
     * Gaussian radial basis function.
     * \param r: input distance.
     * \return: value of radial basis function.
     */
	static double s_basisGaussian(double r);
    
    /*
     * TPS radial basis function.
     * \param r: input distance.
     * \return: value of radial basis function.
     */
	static double s_basisTPS(double r);
    
    /*
     * Determine if a point \param x is in triangle labeled by \param triangleID.
     * \param x: input point.
     * \param triangleID: triangle ID.
     * \return: True if \param x is in triangle labeled by \param triangleID. False otherwise.
     */
    static bool s_isInTriangle(const double *x, int triangleID);
    
    /*
     * Compute square-root of mean error.
     * \return: RMSE.
     */
    static double s_computeRMSE();

private:
    static DBL7VECT *s_centers; // triangle centers
    static const PPMImage *s_imgOrig; // original image
    static PPMImage *s_imgSR; // SR image
    static const TriMesh *s_mesh; // Tri-mesh
	static BasisType s_basis; // type of basis function

    static std::vector<std::vector<int> > s_vf_1ring; // 1-ring VERTEX <-> FACE
    static std::vector<std::set<int> > s_vv_1ring; // 1-ring VERTEX <-> VERTEX
    static std::vector<std::set<int> > s_ff_1ring; // 1-ring FACE <-> FACE
    
    static std::vector<std::vector<int> > s_vf_2ring; // 2-ring VERTEX <-> FACE
    static std::vector<std::set<int> > s_vv_2ring; // 2-ring VERTEX <-> VERTEX
    static std::vector<std::set<int> > s_ff_2ring; // 2-ring FACE <-> FACE
    
    static std::vector<std::pair<double, double> > s_eigenvalues; // FACE <-> transformed eigenvalues
    static std::vector<Eigen::Matrix2d> s_T; // metric matrix
};

#endif /* defined(__arbf_image_super_resolution__ARBFInterpolator__) */
