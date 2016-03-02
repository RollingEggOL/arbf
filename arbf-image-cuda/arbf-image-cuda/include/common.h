//
//  common.h
//  arbf-image-image
//
//  Created by Ke Liu on 10/01/14.
//  Copyright (c) 2014 Ke Liu. All rights reserved.
//

#ifndef ARBF_IMAGE_COMMON_H
#define ARBF_IMAGE_COMMON_H

#include <vector>

//////////////////////////////////////////////////////
// choose one of the two
#define NEIGHBOR_1RING // use 1-ring neighborhood
//#define NEIGHBOR_2RING // use 2-ring neighborhood

//////////////////////////////////////////////////////
//#define ISO_RBF // use isotropic RBF interpolation

//////////////////////////////////////////////////////
//#define WEIGHTED // take weighted average

#define EPSILON (1e-5)
#define SQ(x) ((x) * (x))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// number of threads per thread block
#define BLOCK_SIZE (128)

// device heap size in bytes (for device side malloc and free)
#define DEVICE_HEAP (67108864) // 64 MB

// initial number of pixels per triangle. 
// Make sure this is large enough to hold maximum number of pixels per triangle
#define NUM_PIXEL_PER_FACE (128)

// number of neighboring triangles per triangle.
// Make sure this is large enough to hold (1 + maximum # of neighboring triangles per triangle)
#ifdef NEIGHBOR_1RING
#define NUM_NEIGHBOR_TRIANGLES (33)
#else
#define NUM_NEIGHBOR_TRIANGLES (65)
#endif

#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "ERROR: %s at line %d in file %s\n",				\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		exit(1);															\
		}																	\
}																			\

typedef struct
{
	int f;
	std::vector<int> x;
	std::vector<int> y;
} PixelTag;

typedef struct
{
	double x;
	double y;
} DBL2VECT;

typedef struct
{
	double x;
	double y;
	double intensity;
	double lambda1;
	double lambda2;
	DBL2VECT v1;
	DBL2VECT v2;
} DBL7VECT;

typedef struct
{
	int a;
	int b;
	int c;
	double intensity;
} INT3VECT;

typedef struct
{
	int nv; /* number of vertices */
	int nf; /* number of faces */
	DBL7VECT *vertex;
	INT3VECT *face;
} SurFaceMesh;

#endif