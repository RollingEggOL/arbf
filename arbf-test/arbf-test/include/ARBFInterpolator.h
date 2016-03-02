//
//  ARBFInterpolator.h
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#ifndef __arbf_test_ARBFInterpolator__
#define __arbf_test_ARBFInterpolator__

#include <vector>
#include <set>
#include <tuple>
#include <Eigen/Dense>
#include "Common.h"
#include "PPMImage.h"
#include "TriMesh.h"

class ARBFInterpolator
{
public:
    static const Vertex* getCenters() {
        return s_centers;
    }
    
    static const Vertex* getEdgeCenters() {
        return s_edgeCenters;
    }
    
    static void setMesh(const TriMesh *mesh);
    
    static void setBasisType(const BasisType type=BasisType::MQ);
    
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
     * Calculate coordinates, intensity for triangle centers.
     */
    static void calculateCenters();
    
    /*
     * Calculate coordinates, intensity for triangle edge centers.
     */
    static void calculateEdgeCenters();
    
    /*
     * Find 1-ring VERTEX <-> VERTEX, VERTEX <-> FACE.
     */
    static void findNeighVertices1Ring();
    
    /*
     * Find 1-ring FACE <-> FACE.
     */
    static void findNeighFaces1Ring();
    
    /*
     * Find 2-ring VERTEX <-> VERTEX, VERTEX <-> FACE.
     */
    static void findNeighVertices2Ring();
    
    /*
     * Find 2-ring FACE <-> FACE.
     */
    static void findNeighFaces2Ring();
    
    /*
     * ARBF interpolation.
     * \return: tuple of <result_data, result_maxIntensity, result_dimX, result_dimY, result_dimZ>.
     */
    static std::tuple<double*, double, int, int, int> interpolate();
    
    /*
     * Compute Peak-Signal-Noise-Ratio between original image and reconstructed image.
     * \return: PSNR.
     */
    static double computePSNR(const PPMImage &img1, const PPMImage &img2);
    
    /*
     * Release resources.
     */
    static void clean();
    
private:
    ARBFInterpolator(); // declare as private to prevent creating instances
    
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
     * Compute Euclidean distance.
     * \param v0: input point 1.
     * \param mu2: input point 2.
     * \return: distance betwee \param v0 and \param v1.
     */
    static double s_computeDistance(const double *v0, const double *v1);
    
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
    static double s_computeRMSE(const PPMImage &img1, const PPMImage &img2);
    
    static void buildMatrix(int dim, int nv, int nf, Eigen::MatrixXd &m, Eigen::VectorXd &v);
    static std::vector<double> linspace(double a, double b, int n);

private:
    static Vertex *s_centers; // triangle centers
    static Vertex *s_edgeCenters; // triangle edge centers
    static const TriMesh *s_mesh; // Tri-mesh
	static BasisType s_basis; // type of basis function

    static std::vector<std::vector<int> > s_vf_1ring; // 1-ring VERTEX <-> FACE
    static std::vector<std::set<int> > s_vv_1ring; // 1-ring VERTEX <-> VERTEX
    static std::vector<std::set<int> > s_ff_1ring; // 1-ring FACE <-> FACE
    
    static std::vector<std::vector<int> > s_vf_2ring; // 2-ring VERTEX <-> FACE
    static std::vector<std::set<int> > s_vv_2ring; // 2-ring VERTEX <-> VERTEX
    static std::vector<std::set<int> > s_ff_2ring; // 2-ring FACE <-> FACE
};

#endif /* defined(__arbf_test_ARBFInterpolator__) */
