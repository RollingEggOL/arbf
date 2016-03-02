//
//  main.h
//  arbf-image-image
//
//  Created by Ke Liu on 1/28/14.
//  Copyright (c) 2014 Ke Liu. All rights reserved.
//

#ifndef ARBF_IMAGE_MAIN_H
#define ARBF_IMAGE_MAIN_H

#include <cassert>
#include <cstring>
#include <algorithm>
#include <utility>
#include <vector>
#include <set>
#include <Eigen/Dense>
#include "common.h"

//////////////////////////////////////////////////////
// other macros
#define BUF_SIZE (256) // for file reading

// global constants
//const double SHARPNESS_THRESHOLD = 0.125; // c_bar in paper
const double C0 = 0.2;
//const double GAMMA = 1.0;

// global variables
#ifdef _WIN32
const std::string ROOT = "C:\\Users\\keliu\\Documents\\Visual Studio 2013\\Projects\\arbf\\arbf-image-cuda\\";
#else
const std::string ROOT = "/Users/keliu/Documents/projects/arbf/arbf-image-cuda/";
#endif
unsigned char* orig_img = NULL; // original PPM image
//unsigned char* piecewise = NULL; // piecewise image
DBL7VECT *centers = NULL; // triangle centers
int X = 0, Y = 0; // recovering dimension
SurFaceMesh* samples = NULL; // input sampling mesh

#ifdef NEIGHBOR_1RING
std::vector<std::vector<int> > vertex_faces_1ring; // 1-ring VERTEX <-> FACE map
std::vector<std::set<int> > neigh_vertex_1ring; // 1-ring neighboring vertices for each vertex
std::vector<std::set<int> > face_faces_1ring; // 1-ring FACE <-> FACE map
#elif defined(NEIGHBOR_2RING)
std::vector<std::vector<int> > vertex_faces_2ring; // 2-ring VERTEX <-> FACE map
std::vector<std::set<int> > neigh_vertex_2ring; // 2-ring neighboring vertices for each vertex
std::vector<std::set<int> > face_faces_2ring; // 2-ring FACE <-> FACE map
#endif

std::vector<std::pair<double, double> > eigenvalues; // transformed eigenvalues for each vertex
std::vector<Eigen::Matrix2d> T; // metric matrix
//double *rho = NULL; // defines the size of support domain for each vertex
double *g_F = NULL; // interpolated intensities

// function declarations
void clean_up();
void find_neigh_vertex_1ring();
void find_neigh_face_1ring();
#ifdef NEIGHBOR_2RING
void find_neigh_vertex_2ring();
void find_neigh_face_2ring();
#endif
void transform_eigenvalues();
std::pair<DBL2VECT, DBL2VECT> compute_eigen(const Eigen::Matrix2d &m, double *lambda1, double *lambda2);
void compute_metrics();
void compute_rho();
void interpolate(); // ARBF interpolation
#ifdef WEIGHTED
void weight_interpolate(); // weighted ARBF interpolation
#endif
void find_enclosing_triangle();
void rescale_image();
void write_image();
void clean_mesh(SurFaceMesh **mesh);
void read_mesh_tri(const char *file_name, SurFaceMesh **surfmesh);
int init_device(); // find and initialize device
extern void interpolate_helper(int metricDimension, int maxNumNeighFace, const double *T, const int *face_face);

// for sharpness parameter c
inline
double compute_c(double mu1, double mu2)
{
//    double maxMu = std::max(mu1, mu2);
    double sumMu = mu1 + mu2;
    
//    if (maxMu > EPSILON || maxMu < -EPSILON)
    if (sumMu > EPSILON || sumMu < -EPSILON)
    {
//        return ((mu1 - mu2) / maxMu);
        return abs(mu1 - mu2) / sumMu;
    }
    else
    {
        return 1.0;
    }
}

// compute Rooted Mean Squared Error
inline
double compute_RMSE()
{
    double sum = 0.0;
    
    for (int j = 0; j < Y; ++j)
    {
        for (int i = 0; i < X; ++i)
        {
            sum += SQ((double) orig_img[j*X+i] - g_F[j*X+i]);
//            sum += SQ((double) orig_img[j*N1+i] - (double) piecewise[j*N1+i]);
        }
    }
    return sqrt(sum / (X * Y));
}

// compute Peak Signal-to-Noise Ratio
inline
double compute_PSNR()
{
    return 20.0 * log10(255.0 / compute_RMSE());
}

#ifdef WEIGHTED
inline
double compute_weight(double r)
{
    return exp(-SQ(r) / SQ(30.0));
}
#endif

#endif
