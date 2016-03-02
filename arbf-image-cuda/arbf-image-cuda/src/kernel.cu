//
//  kernel.cu
//  arbf-image-cuda
//
//  Created by Ke Liu on 10/01/14.
//  Copyright (c) 2014 Ke Liu. All rights reserved.
//

#include <vector>
#include <set>
#include <cuda_runtime.h>
#include "../include/common.h"

extern std::vector<PixelTag> g_pixel_tag;
extern SurFaceMesh* samples;
extern DBL7VECT* centers;
extern int X;
extern int Y;
extern double *g_F;

extern std::vector<std::vector<int> > vertex_faces_1ring;
extern std::vector<std::set<int> > neigh_vertex_1ring;
extern std::vector<std::set<int> > face_faces_1ring;
#ifdef NEIGHBOR_2RING
extern std::vector<std::vector<int> > vertex_faces_2ring;
extern std::vector<std::set<int> > neigh_vertex_2ring;
extern std::vector<std::set<int> > face_faces_2ring;
#endif

#define mat_elem(a, y, x, n) (a + ((y) * (n) + (x)))

class vector2d
{
public:
	__device__ vector2d(double d1, double d2)
	{
		_data.x = make_double2(d1, d2).x;
		_data.y = make_double2(d1, d2).y;
	}

	// getter
	__device__ const double2 data() const
	{
		return _data;
	}

	__device__ double x() const
	{
		return _data.x;
	}

	__device__ double y() const
	{
		return _data.y;
	}

	// scalar multiplication
	__device__ vector2d& mul(double s)
	{
		_data.x *= s;
		_data.y *= s;

		return *this;
	}

	// norm of vector of size 2
	__device__ double norm()
	{
		return sqrt(SQ(_data.x) + SQ(_data.y));
	}

	// dot product of vector of size 2
	__device__ double dot(const vector2d& v2)
	{
		return _data.x * v2.x() + _data.y * v2.y();
	}

private:
	double2 _data;
};

// swap row r1 and r2 of matrix A and corresponding elements of vector b
__device__ void swap_row(double *A, double *b, int r1, int r2, int n)
{
	double tmp, *p1, *p2;
	int i;
	
	if (r1 == r2)
		return;

	for (i = 0; i < n; i++)
	{
		p1 = mat_elem(A, r1, i, n);
		p2 = mat_elem(A, r2, i, n);
		tmp = *p1;
		*p1 = *p2;
		*p2 = tmp;
	}
	tmp = b[r1];
	b[r1] = b[r2];
	b[r2] = tmp;
}

// solve linear system A*x = b
__device__ void solve(double *x, double* A, double* b, int n)
{
#define A(y, x) (*mat_elem(A, y, x, n))
	int maxRow, row, col, dia, j;
	double max, tmp;

	// Gaussian elimination
	for (dia = 0; dia < n; dia++)
	{
		// find maxRow and max diagonal element
		maxRow = dia;
		max = A(dia, dia);
		for (row = dia + 1; row < n; row++)
		{
			if ((tmp = fabs(A(row, dia))) > max)
			{
				maxRow = row;
				max = tmp;
			}
		}
		swap_row(A, b, dia, maxRow, n);

		// eleminate
		for (row = dia + 1; row < n; row++)
		{
			tmp = A(row, dia) / A(dia, dia);
			for (col = dia + 1; col < n; col++)
				A(row, col) -= (tmp * A(dia, col));
			A(row, dia) = 0.0;
			b[row] -= (tmp * b[dia]);
		}
	}

	// back-substitution
	for (row = n - 1; row >= 0; row--)
	{
		tmp = b[row];
		for (j = n - 1; j > row; j--)
			tmp -= (x[j] * A(row, j));
		x[row] = tmp / A(row, row);
	}
#undef A
}

__device__ bool is_in_triangle(int _xx, int _xy, double _ax, double _ay, double _bx, double _by, double _cx, double _cy)
{
	vector2d ab(_bx - _ax, _by - _ay);
	vector2d bc(_cx - _bx, _cy - _by);
	vector2d ca(_ax - _cx, _ay - _cy);
	vector2d ax(_xx - _ax, _xy - _ay);
	vector2d bx(_xx - _bx, _xy - _by);
	vector2d cx(_xx - _cx, _xy - _cy);

	double area = 0.0, area1 = 0.0, area2 = 0.0, area3 = 0.0;
	double cosTheta = ab.dot(ca.mul(-1)) / (ab.norm()*ca.norm());
	area = 0.5*ab.norm()*ca.norm()*sqrt(1.0 - SQ(cosTheta));
	double cosTheta1 = ab.dot(ax) / (ab.norm()*ax.norm());
	area1 = 0.5*ab.norm()*ax.norm()*sqrt(1.0 - SQ(cosTheta1));
	double cosTheta2 = cx.dot(ca) / (cx.norm()*ca.norm());
	area2 = 0.5*cx.norm()*ca.norm()*sqrt(1.0 - SQ(cosTheta2));
	double cosTheta3 = bx.dot(bc) / (bx.norm()*bc.norm());
	area3 = 0.5*bx.norm()*bc.norm()*sqrt(1.0 - SQ(cosTheta3));

	if (isnan(area)) area = 0.0;
	if (isnan(area1)) area1 = 0.0;
	if (isnan(area2)) area2 = 0.0;
	if (isnan(area3)) area3 = 0.0;

	if (abs(area1 + area2 + area3 - area) < EPSILON)
		return true;
	else
		return false;
}

// distance under metric T
__device__ double compute_distance(const double *v0, const double *v1, const double *T)
{
	vector2d dis(v0[0] - v1[0], v0[1] - v1[1]);
	vector2d tmp(dis.x() * T[0] + dis.y() * T[2], dis.x() * T[1] + dis.y() * T[3]);
	
	return sqrt(dis.dot(tmp));
}

// basis function MQ
__device__ double phi(double r)
{
	return sqrt(SQ(r) + SQ(0.5));
}

// basis function IMQ
//__device__ double phi(double r)
//{
//    return 1.0 / (sqrt(SQ(r) + SQ(1.8)));
//}

// basis function Gaussian
//__device__ double phi(double r)
//{
//    double c = 0.01;
//    return exp(-SQ(c*r));
//}

// basis function TPS
//__device__ double phi(double r)
//{
//    if (abs(r) < EPSILON)
//        return 0;
//    else
//        return SQ(r) * log(r);
//}

__global__ void find_triangle_kernel(int* dx, int* dy, size_t stride, size_t blockStride,
	const DBL7VECT* vertex, const INT3VECT* face, int N)
{
	extern __shared__ int sdata[];

	int tx = blockIdx.x*blockDim.x + threadIdx.x;

	if (tx < N)
	{
		// face[tx].a = sdata[threadIdx.x*blockStride + 10]
		// face[tx].b = sdata[threadIdx.x*blockStride + 11]
		// face[tx].c = sdata[threadIdx.x*blockStride + 12]
		sdata[threadIdx.x*blockStride + 10] = face[tx].a;
		sdata[threadIdx.x*blockStride + 11] = face[tx].b;
		sdata[threadIdx.x*blockStride + 12] = face[tx].c;

		// v0x = sdata[threadIdx.x*blockStride]
		// v0y = sdata[threadIdx.x*blockStride+1]
		// v1x = sdata[threadIdx.x*blockStride+2]
		// v1y = sdata[threadIdx.x*blockStride+3]
		// v2x = sdata[threadIdx.x*blockStride+4]
		// v2y = sdata[threadIdx.x*blockStride+5]
		sdata[threadIdx.x*blockStride] = vertex[sdata[threadIdx.x*blockStride + 10]].x;
		sdata[threadIdx.x*blockStride + 1] = vertex[sdata[threadIdx.x*blockStride + 10]].y;
		sdata[threadIdx.x*blockStride + 2] = vertex[sdata[threadIdx.x*blockStride + 11]].x;
		sdata[threadIdx.x*blockStride + 3] = vertex[sdata[threadIdx.x*blockStride + 11]].y;
		sdata[threadIdx.x*blockStride + 4] = vertex[sdata[threadIdx.x*blockStride + 12]].x;
		sdata[threadIdx.x*blockStride + 5] = vertex[sdata[threadIdx.x*blockStride + 12]].y;

		// minX = sdata[threadIdx.x*blockStride+6]
		// maxX = sdata[threadIdx.x*blockStride+7]
		// minY = sdata[threadIdx.x*blockStride+8]
		// maxY = sdata[threadIdx.x*blockStride+9]
		sdata[threadIdx.x*blockStride + 6] = MIN(sdata[threadIdx.x*blockStride], sdata[threadIdx.x*blockStride + 2]);
		sdata[threadIdx.x*blockStride + 6] = MIN(sdata[threadIdx.x*blockStride + 6], sdata[threadIdx.x*blockStride + 4]);
		sdata[threadIdx.x*blockStride + 7] = MAX(sdata[threadIdx.x*blockStride], sdata[threadIdx.x*blockStride + 2]);
		sdata[threadIdx.x*blockStride + 7] = MAX(sdata[threadIdx.x*blockStride + 7], sdata[threadIdx.x*blockStride + 4]);
		sdata[threadIdx.x*blockStride + 8] = MIN(sdata[threadIdx.x*blockStride + 1], sdata[threadIdx.x*blockStride + 3]);
		sdata[threadIdx.x*blockStride + 8] = MIN(sdata[threadIdx.x*blockStride + 8], sdata[threadIdx.x*blockStride + 5]);
		sdata[threadIdx.x*blockStride + 9] = MAX(sdata[threadIdx.x*blockStride + 1], sdata[threadIdx.x*blockStride + 3]);
		sdata[threadIdx.x*blockStride + 9] = MAX(sdata[threadIdx.x*blockStride + 9], sdata[threadIdx.x*blockStride + 5]);

		int n = 0; // actual number of pixels in current triangle
		int k = 1; // indice of pixel coord

		for (int j = sdata[threadIdx.x*blockStride + 8]; j <= sdata[threadIdx.x*blockStride + 9]; ++j)
		{
			for (int i = sdata[threadIdx.x*blockStride + 6]; i <= sdata[threadIdx.x*blockStride + 7]; ++i)
			{
				int p1 = i;
				int p2 = j;
				if (is_in_triangle(p1, p2,
					sdata[threadIdx.x*blockStride], sdata[threadIdx.x*blockStride + 1], sdata[threadIdx.x*blockStride + 2], sdata[threadIdx.x*blockStride + 3],
					sdata[threadIdx.x*blockStride + 4], sdata[threadIdx.x*blockStride + 5]))
				{
					n++;
					dx[tx*stride] = n;
					dy[tx*stride] = n;
					dx[tx*stride + k] = i;
					dy[tx*stride + k] = j;
					k++;
				}
			}
		}
	}
}

__global__ void interpolate_kernel(double *F, int numFaces, int X, int Y, int numMetricElem, int maxNumNeighFace,
	const DBL7VECT* vertices, const INT3VECT* faces, const DBL7VECT* centers, const double *T, const int *face_face)
{
	extern __shared__ int s_face_face[]; // face_face on shared memory
	int tx = blockIdx.x * blockDim.x + threadIdx.x; // global thread ID
	int tid = threadIdx.x; // local thread ID
	int x[NUM_PIXEL_PER_FACE], y[NUM_PIXEL_PER_FACE]; // coords of pixels in current triangle

	if (tx < numFaces)
	{
		// find TRIANGLE <-> PIXEL mapping
		const INT3VECT *faceVertices = &(faces[tx]);
		const DBL7VECT *v0 = &(vertices[faceVertices->a]);
		const DBL7VECT *v1 = &(vertices[faceVertices->b]);
		const DBL7VECT *v2 = &(vertices[faceVertices->c]);
		double minX = MIN(MIN(v0->x, v1->x), v2->x);
		double maxX = MAX(MAX(v0->x, v1->x), v2->x);
		double minY = MIN(MIN(v0->y, v1->y), v2->y);
		double maxY = MAX(MAX(v0->y, v1->y), v2->y);
		int counter = 0; // stores how many pixels are actually in current face
		for (int j = (int) minY; j <= (int) ceil(maxY); ++j)
		{
			for (int i = (int) minX; i <= (int) ceil(maxX); ++i)
			{
				if (is_in_triangle(i, j, v0->x, v0->y, v1->x, v1->y, v2->x, v2->y))
				{
					x[counter] = i;
					y[counter] = j;
					counter ++;
				}
			}
		}

		// solve linear equations (ma * coeff = u)
		int startIdx_face = tid * maxNumNeighFace;
		s_face_face[startIdx_face] = tx * NUM_NEIGHBOR_TRIANGLES; // start index of current thread for face_face
		int numNeigh = face_face[s_face_face[startIdx_face]];
		double *ma = (double*) malloc(sizeof(double) * SQ(numNeigh)); // distance matrix, row-majored
		double *u = (double*) malloc(sizeof(double) * numNeigh); // right hand side vector
		double *coeff = (double*) malloc(sizeof(double) * numNeigh); // coefficient vector
		for (int i = 0; i < numNeigh; ++i)
		{
			s_face_face[startIdx_face + 1 + i] = face_face[s_face_face[startIdx_face] + 1 + i];
		}
		for (int i = 0; i < numNeigh; ++i)
		{
			const DBL7VECT *c1 = &(centers[s_face_face[startIdx_face + 1 + i]]);
			u[i] = c1->intensity;
			for (int j = 0; j < numNeigh; ++j)
			{
				int startIdx_T = s_face_face[startIdx_face + 1 + j] * numMetricElem;
				const DBL7VECT *c2 = &(centers[s_face_face[startIdx_face + 1 + j]]);
				double p1[2] = { c1->x, c1->y };
				double p2[2] = { c2->x, c2->y };
				double metric[4] = {T[startIdx_T], T[startIdx_T + 1], T[startIdx_T + 2], T[startIdx_T + 3]};
				double r = compute_distance(p1, p2, metric);
				ma[i * numNeigh + j] = phi(r);
			}
		}
		solve(coeff, ma, u, numNeigh); // solve coeff by Gaussian elimination
		// do interpolation
		for (int k1 = 0; k1 < counter; ++k1)
		{
			int m = x[k1];
			int n = y[k1];
			double p1[2] = { (double) m, (double) n };
			double intensity = 0.0;
			for (int k2 = 0; k2 < numNeigh; ++k2)
			{
				const DBL7VECT *center = &(centers[s_face_face[startIdx_face + 1 + k2]]);
				double p2[2] = { center->x, center->y };
				int startIdx_T = s_face_face[startIdx_face + 1 + k2] * numMetricElem;
				double metric[4] = { T[startIdx_T], T[startIdx_T + 1], T[startIdx_T + 2], T[startIdx_T + 3] };
				double r = compute_distance(p1, p2, metric);
				intensity += (coeff[k2] * phi(r));
			}
			F[n*X + m] = intensity; // sum interpolated intensity
		}
		free(ma);
		free(u);
		free(coeff);
	} // end if
}

__global__ void threshold_kernel(double *F, int N)
{
	int tx = blockIdx.x * blockDim.x + threadIdx.x;
	if (tx < N)
	{
		if (F[tx] < 0.0)
			F[tx] = 0.0;
		if (F[tx] > 255.0)
			F[tx] = 255.0;
	}
}

__global__ void testKernel(double *x)
{
	double a[] = { 
		1.00, 0.00, 0.00, 0.00, 0.00, 0.00,
		1.00, 0.63, 0.39, 0.25, 0.16, 0.10,
		1.00, 1.26, 1.58, 1.98, 2.49, 3.13,
		1.00, 1.88, 3.55, 6.70, 12.62, 23.80,
		1.00, 2.51, 6.32, 15.88, 39.90, 100.28,
		1.00, 3.14, 9.87, 31.01, 97.41, 306.02 };
	double b[] = { -0.01, 0.61, 0.91, 0.99, 0.60, 0.02 };
	solve(x, a, b, 6);
}

void interpolate_helper(int numMetricElem, int maxNumNeighFace, const double *T, const int *face_face)
{
	size_t bytesImage = sizeof(double) * X * Y; // size of image in bytes
	size_t bytesCenters = sizeof(DBL7VECT) * samples->nf; // size of triangle centers in bytes
	size_t bytesVertices = sizeof(DBL7VECT) * samples->nv; // size of vertices in bytes
	size_t bytesFaces = sizeof(INT3VECT) * samples->nf; // size of faces in bytes
	g_F = (double*) malloc(bytesImage);
	memset(g_F, 0, bytesImage);
	//CUDA_CHECK_RETURN(cudaHostAlloc(&g_F, bytesImage, cudaHostAllocDefault)); // allocate pinned-memory for g_F on host
	//CUDA_CHECK_RETURN(cudaMemset((void*)g_F, 0, bytesImage));
	double *d_F;
	DBL7VECT *d_Centers, *d_Vertices;
	INT3VECT *d_Faces;
	CUDA_CHECK_RETURN(cudaMalloc(&d_F, bytesImage)); // allocate memory on device for d_F
	CUDA_CHECK_RETURN(cudaMalloc(&d_Centers, bytesCenters)); // allocate memory on device for d_Centers
	CUDA_CHECK_RETURN(cudaMalloc(&d_Vertices, bytesVertices)); // allocate memory on device for d_Vertices
	CUDA_CHECK_RETURN(cudaMalloc(&d_Faces, bytesFaces)); // allocate memory on device for d_Faces
	CUDA_CHECK_RETURN(cudaMemcpy(d_Vertices, samples->vertex, bytesVertices, cudaMemcpyHostToDevice)); // copy samples->vertex
	CUDA_CHECK_RETURN(cudaMemcpy(d_Faces, samples->face, bytesFaces, cudaMemcpyHostToDevice)); // copy samples->face
	CUDA_CHECK_RETURN(cudaMemcpy(d_Centers, centers, bytesCenters, cudaMemcpyHostToDevice)); // copy centers

	// override default device heap size (8MB)
	CUDA_CHECK_RETURN(cudaDeviceSetLimit(cudaLimitMallocHeapSize, DEVICE_HEAP));

	// interpolate
	unsigned gridSize = samples->nf / BLOCK_SIZE + ((samples->nf % BLOCK_SIZE) == 0 ? 0 : 1);
	unsigned byteShardMem = sizeof(int) * BLOCK_SIZE * maxNumNeighFace; // size of shared memory in bytes
	interpolate_kernel <<<gridSize, BLOCK_SIZE, byteShardMem>>> (
		d_F, samples->nf, X, Y, numMetricElem, maxNumNeighFace, d_Vertices, d_Faces, d_Centers, T, face_face);
	/*
	double *h_x = (double*)malloc(sizeof(double) * 6);
	double *d_x;
	CUDA_CHECK_RETURN(cudaMalloc(&d_x, sizeof(double) * 6));
	testKernel <<<1, 1 >>> (d_x);
	CUDA_CHECK_RETURN(cudaDeviceSynchronize());
	CUDA_CHECK_RETURN(cudaMemcpy(h_x, d_x, sizeof(double) * 6, cudaMemcpyDeviceToHost));
	for (int i = 0; i < 6; ++i)
	{
		printf("%g\n", h_x[i]);
	}
	free(h_x);
	CUDA_CHECK_RETURN(cudaFree(d_x));*/

	// take average and set threshold
	int newBlockSize = 2 * BLOCK_SIZE;
	gridSize = X * Y / newBlockSize + ((X * Y % newBlockSize) == 0 ? 0 : 1);
	CUDA_CHECK_RETURN(cudaDeviceSynchronize());
	threshold_kernel << <gridSize, newBlockSize >> > (d_F, X * Y);
	CUDA_CHECK_RETURN(cudaDeviceSynchronize());

	// copy result back to host
	CUDA_CHECK_RETURN(cudaMemcpy(g_F, d_F, bytesImage, cudaMemcpyDeviceToHost));

	

	//for (int j = 0; j < Y; ++j)
	//{
	//	for (int i = 0; i < X; ++i)
	//	{
	//		if ((g_F[j*X + i]) < 0.0)
	//			g_F[j*X + i] = 0.0;
	//		if ((g_F[j*X + i]) > 255.0)
	//			g_F[j*X + i] = 255.0;
	//	}
	//}

	// clean up
	CUDA_CHECK_RETURN(cudaFree(d_F));
	CUDA_CHECK_RETURN(cudaFree(d_Centers));
	CUDA_CHECK_RETURN(cudaFree(d_Vertices));
	CUDA_CHECK_RETURN(cudaFree(d_Faces));
}
