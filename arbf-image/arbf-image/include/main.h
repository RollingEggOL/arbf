//
//  main.h
//  arbf-image
//
//  Created by Ke Liu on 1/28/14.
//  Copyright (c) 2014 Ke Liu. All rights reserved.
//

#ifndef arbf_image_main_h
#define arbf_image_main_h

#include <cassert>
#include <cstring>
#include <algorithm>
#include <utility>
#include <vector>
#include <set>
#include <Eigen/Dense>
#include "../src/mesh_parser.c"

//////////////////////////////////////////////////////
// choose one of the two
#define NEIGHBOR_1RING // use 1-ring neighborhood
//#define NEIGHBOR_2RING // use 2-ring neighborhood

//////////////////////////////////////////////////////
//#define ISO_RBF // use isotropic RBF interpolation

//////////////////////////////////////////////////////
//#define WEIGHTED // take weighted average

//////////////////////////////////////////////////////
// other macros
#define SQ(x) ((x) * (x))

// typedefs
typedef struct
{
    int f;
    std::vector<int> x;
    std::vector<int> y;
} PixelTag;

// global constants
const double EPSILON = 1e-5;
//const double SHARPNESS_THRESHOLD = 0.125; // c_bar in paper
const double C0 = 0.2;
//const double GAMMA = 1.0;

// global variables
#ifdef _WIN32
const std::string ROOT = "C:\\Users\\keliu\\Documents\\Visual Studio 2013\\Projects\\arbf\\arbf-image\\";
#else
const std::string ROOT = "/Users/keliu/Documents/projects/arbf/arbf-image/";
#endif

int X = 0, Y = 0; // recovering dimension
SurFaceMesh* samples = NULL; // input sampling mesh
unsigned char* orig_img = NULL; // original PPM image
//unsigned char* piecewise = NULL; // piecewise image
DBL7VECT *centers = NULL; // triangle centers

std::vector<std::vector<int> > vertex_faces_1ring; // 1-ring VERTEX <-> FACE map
std::vector<std::set<int> > neigh_vertex_1ring; // 1-ring neighboring vertices for each vertex
std::vector<std::set<int> > face_faces_1ring; // 1-ring FACE <-> FACE map
#ifdef NEIGHBOR_2RING
std::vector<std::vector<int> > vertex_faces_2ring; // 2-ring VERTEX <-> FACE map
std::vector<std::set<int> > neigh_vertex_2ring; // 2-ring neighboring vertices for each vertex
std::vector<std::set<int> > face_faces_2ring; // 2-ring FACE <-> FACE map
#endif

std::vector<std::pair<double, double> > eigenvalues; // transformed eigenvalues for each vertex
std::vector<Eigen::Matrix2d> T; // metric matrix
//double *rho = NULL; // defines the size of support domain for each vertex
double *g_F = NULL; // product of interpolated intensities

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

// distance under metric T
inline
double compute_distance(const double* v0, const double *v1, const Eigen::Matrix2d &T)
{
    Eigen::RowVector2d dis(v0[0] - v1[0], v0[1] - v1[1]);
	Eigen::RowVector2d tmp = dis * T;
	
	return std::sqrt(tmp(0)*dis(0) + tmp(1)*dis(1));
}

// basis function MQ
inline
double phi(double r)
{
    return std::sqrt(SQ(r) + SQ(0.5));
}

// basis function IMQ
//inline
//double phi(double r)
//{
//    return 1.0 / (std::sqrt(SQ(r) + SQ(1.8)));
//}

// basis function Gaussian
//inline
//double phi(double r)
//{
//    double c = 0.01;
//    return std::exp(-SQ(c*r));
//}

// basis function TPS
//inline
//double phi(double r)
//{
//    if (abs(r) < EPSILON)
//        return 0;
//    else
//        return SQ(r) * std::log(r);
//}

inline
bool is_in_triangle(const int* x, int triangleID)
{
    DBL7VECT &a = samples->vertex[samples->face[triangleID].a];
    DBL7VECT &b = samples->vertex[samples->face[triangleID].b];
    DBL7VECT &c = samples->vertex[samples->face[triangleID].c];
    
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
    
    if (abs(area1 + area2 + area3 - area) < EPSILON)
        return true;
    else
        return false;
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
