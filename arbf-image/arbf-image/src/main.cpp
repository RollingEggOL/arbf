//
//  main.cpp
//  arbf-image
//
//  Created by Ke Liu on 1/23/14.
//  Copyright (c) 2014 Ke Liu. All rights reserved.
//

#include <iostream>
#include <iomanip> // for std::setw()
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include "../include/main.h"

//#define __TEST__

#ifdef __TEST__
#include <fstream>
#endif

using namespace std;
using namespace Eigen;

void clean_up()
{
	// clean mesh
	clean_mesh(&samples);
	
	// clean original image
	if (orig_img)
	{
		free(orig_img);
		orig_img = NULL;
	}
	
	// clean triangle centers
	if (centers)
	{
		free(centers);
		centers = NULL;
	}

	// clean g_F
	if (g_F)
    {
        free(g_F);
        g_F = NULL;
    }
}

void transform_eigenvalues()
{
	eigenvalues.resize(samples->nf);
	
#ifdef __TEST__
	ofstream ofs1((ROOT + string("DEBUG.transformed_eigenvalue.txt")).c_str());
#endif

	for (int f = 0; f < samples->nf; ++f)
	{
        double eig2 = 0.0;
        double c = compute_c(centers[f].lambda1, centers[f].lambda2);
		if (c != 0.0)
        {
            eig2 = 1.0 - exp(-3.8 * SQ(C0/c));
        }
        else
        {
            eig2 = 1.0;
        }
        
#ifdef ISO_RBF
        eig2 = 1.0;
#endif
        
		eigenvalues[f] = make_pair(1.0, eig2);
		
#ifdef __TEST__
		ofs1 << "f" << f << ": "
        << eigenvalues[f].first << "\t"
        << eigenvalues[f].second << "\n";
#endif
	}
	
#ifdef __TEST__
	ofs1.close();
#endif
}

pair<DBL2VECT, DBL2VECT> compute_eigen(const Eigen::Matrix2d &m, double *lambda1, double *lambda2)
{
	assert(abs(m(0, 1) - m(1, 0)) < EPSILON); // m should be symmetric
	
	// compute eigenvalues
	// matrix looks like
	//      | a    c|
	// A =  |       |
	//      | c    b|
	
	double a = m(0, 0), b = m(1, 1), c = m(0, 1);
	*lambda1 = 0.5 * (a + b) + sqrt(SQ(c) + 0.25 * SQ(a - b));
	*lambda2 = 0.5 * (a + b) - sqrt(SQ(c) + 0.25 * SQ(a - b));
	
	// compute eigenvectors
	DBL2VECT v1, v2;
	if (abs(c) > EPSILON) // c != 0
	{
//        double u1 = c;
//        double w1 = *lambda1 - a;
		
		double u1 = *lambda1 - b;
		double w1 = c;
		
		double p = sqrt(SQ(u1) + SQ(w1));
		v1.x = u1 / p;
		v1.y = w1 / p;
        
        v2.x = v1.y;
        v2.y = -v1.x;
	}
	else
	{
		v1.x = 1.0;
		v1.y = 0.0;
		v2.x = 0.0;
		v2.y = 1.0;
	}
	
	return make_pair(v1, v2);
}

void compute_metrics()
{
	T.resize(samples->nf);
	
	for (int f = 0; f < samples->nf; ++f)
	{
		DBL7VECT *center = &centers[f];
		//             |lambda1 0      | |v1T|
		// T = [v1 v2] |               | |   |
		//             |0       lambda2| |v2T|
        
		Matrix2d t1, t2;
		t1(0, 0) = center->v1.x;
        t1(1, 0) = center->v1.y;
		t1(0, 1) = center->v2.x;
		t1(1, 1) = center->v2.y;
		
		t2(0, 0) = eigenvalues[f].first;
		t2(0, 1) = 0.0;
		t2(1, 0) = 0.0;
		t2(1, 1) = eigenvalues[f].second;
		
		T[f] = t1 * t2;
		T[f] = T[f] * t1.transpose();
	}
}

#ifdef NEIGHBOR_2RING
void find_neigh_face_2ring()
{
    face_faces_2ring.resize(samples->nf);
    
    for (int f = 0; f < samples->nf; ++f)
    {
        int a = samples->face[f].a;
        int b = samples->face[f].b;
        int c = samples->face[f].c;
        
        for (int j = 0; j < vertex_faces_2ring[a].size(); ++j)
        {
            face_faces_2ring[f].insert(vertex_faces_2ring[a][j]);
        }
        
        for (int j = 0; j < vertex_faces_2ring[b].size(); ++j)
        {
            face_faces_2ring[f].insert(vertex_faces_2ring[b][j]);
        }
        
        for (int j = 0; j < vertex_faces_2ring[c].size(); ++j)
        {
            face_faces_2ring[f].insert(vertex_faces_2ring[c][j]);
        }
    }
    
#ifdef __TEST__
    ofstream ofs(ROOT + "DEBUG.face_faces_2ring.txt");
    
	if (! ofs.is_open())
	{
		cerr << "ERROR: open DEBUG.face_faces_2ring.txt failed ..." << endl;
		exit(EXIT_FAILURE);
	}
    
    for (int f = 0; f < face_faces_2ring.size(); ++f)
    {
        ofs << "[f <-> f] [" << face_faces_2ring[f].size() << "] f" << f << ": ";
        for (set<int>::const_iterator it = face_faces_2ring[f].begin();
             it != face_faces_2ring[f].end(); ++it)
        {
            ofs << *it << " ";
        }
        ofs << endl;
    }
    
    ofs.close();
#endif
}

void find_neigh_vertex_2ring()
{
	neigh_vertex_2ring.resize(samples->nv);
	vertex_faces_2ring.resize(samples->nv);
	
	// build 2-ring VERTEX <-> FACE mapping
	vertex_faces_2ring.assign(vertex_faces_1ring.begin(), vertex_faces_1ring.end());
	for (int i = 0; i < samples->nf; ++i)
	{
		// 3 vertices of current face
		int a = samples->face[i].a;
		int b = samples->face[i].b;
		int c = samples->face[i].c;
		
		// neighboring vertices of 3 vertices
		const std::set<int> &neib1 = neigh_vertex_1ring[a];
		const std::set<int> &neib2 = neigh_vertex_1ring[b];
		const std::set<int> &neib3 = neigh_vertex_1ring[c];
		
		// go over each neighboring vertices set, check if current face is in
		// the neighboring faces list
		std::set<int>::const_iterator it;
		for (it = neib1.begin(); it != neib1.end(); ++it)
		{
			std::vector<int> &neib_faces = vertex_faces_1ring[*it];
			if (find(neib_faces.begin(), neib_faces.end(), i) == neib_faces.end())
			{
				// current face is not in the list, add it
				vertex_faces_2ring[*it].push_back(i);
			}
		}
		
		for (it = neib2.begin(); it != neib2.end(); ++it)
		{
			std::vector<int> &neib_faces = vertex_faces_1ring[*it];
			if (find(neib_faces.begin(), neib_faces.end(), i) == neib_faces.end())
			{
				// current face is not in the list, add it
				vertex_faces_2ring[*it].push_back(i);
			}
		}
		
		for (it = neib3.begin(); it != neib3.end(); ++it)
		{
			std::vector<int> &neib_faces = vertex_faces_1ring[*it];
			if (find(neib_faces.begin(), neib_faces.end(), i) == neib_faces.end())
			{
				// current face is not in the list, add it
				vertex_faces_2ring[*it].push_back(i);
			}
		}
	}
	
	// find 2-ring neighboring vertices
	neigh_vertex_2ring.assign(neigh_vertex_1ring.begin(), neigh_vertex_1ring.end());
	for (int i = 0; i < samples->nv; ++i)
	{
		for (int j = 0; j < vertex_faces_2ring[i].size(); ++j)
		{
			neigh_vertex_2ring[i].insert(samples->face[vertex_faces_2ring[i][j]].a);
			neigh_vertex_2ring[i].insert(samples->face[vertex_faces_2ring[i][j]].b);
			neigh_vertex_2ring[i].insert(samples->face[vertex_faces_2ring[i][j]].c);
		}
	}
	
#ifdef __TEST__
	ofstream ofs("DEBUG.vertex_faces_2ring.txt");
	
	if (! ofs.is_open())
	{
		cerr << "ERROR: open DEBUG.vertex_faces_2ring.txt failed ..." << endl;
		exit(EXIT_FAILURE);
	}
	
	for (int i = 0; i < vertex_faces_2ring.size(); ++i)
	{
		ofs << "[v <-> f] [" << vertex_faces_2ring[i].size() << "] v" << i << ": ";
		for (int j = 0; j < vertex_faces_2ring[i].size(); ++j)
		{
			ofs << vertex_faces_2ring[i][j] << " ";
		}
		ofs << endl;
	}
	
	ofs.close();
#endif
	
#ifdef __TEST__
	ofstream ofs2("DEBUG.neigh_vertex_2ring.txt");
	
	if (! ofs2.is_open())
	{
		cerr << "ERROR: open DEBUG.neigh_vertex_2ring.txt failed ..." << endl;
		exit(EXIT_FAILURE);
	}
	
	for (int i = 0; i < neigh_vertex_2ring.size(); ++i)
	{
		ofs2 << "[v <-> v] [" << neigh_vertex_2ring[i].size() << "] v" << i << ": ";
		for (set<int>::const_iterator it = neigh_vertex_2ring[i].begin();
			 it != neigh_vertex_2ring[i].end(); ++it)
		{
			ofs2 << *it << " ";
		}
		ofs2 << endl;
	}
	
	ofs2.close();
#endif
}
#endif

void find_neigh_face_1ring()
{
    face_faces_1ring.resize(samples->nf);
    
    for (int f = 0; f < samples->nf; ++f)
    {
        int a = samples->face[f].a;
        int b = samples->face[f].b;
        int c = samples->face[f].c;
        
        for (size_t j = 0; j < vertex_faces_1ring[a].size(); ++j)
        {
            face_faces_1ring[f].insert(vertex_faces_1ring[a][j]);
        }
        
        for (size_t j = 0; j < vertex_faces_1ring[b].size(); ++j)
        {
            face_faces_1ring[f].insert(vertex_faces_1ring[b][j]);
        }
        
        for (size_t j = 0; j < vertex_faces_1ring[c].size(); ++j)
        {
            face_faces_1ring[f].insert(vertex_faces_1ring[c][j]);
        }
    }
    
#ifdef __TEST__
    ofstream ofs((ROOT + string("DEBUG.face_faces_1ring.txt")).c_str());
    
	if (! ofs.is_open())
	{
		cerr << "ERROR: open DEBUG.face_faces_1ring.txt failed ..." << endl;
		exit(EXIT_FAILURE);
	}
    
    for (int f = 0; f < face_faces_1ring.size(); ++f)
    {
        ofs << "[f <-> f] [" << face_faces_1ring[f].size() << "] f" << f << ": ";
        for (set<int>::const_iterator it = face_faces_1ring[f].begin();
             it != face_faces_1ring[f].end(); ++it)
        {
            ofs << *it << " ";
        }
        ofs << endl;
    }

    ofs.close();
#endif
}

void find_neigh_vertex_1ring()
{
	neigh_vertex_1ring.resize(samples->nv);
	vertex_faces_1ring.resize(samples->nv);
	
	// build 1-ring VERTEX <-> FACE mapping
	for (int i = 0; i < samples->nf; ++i)
	{
		vertex_faces_1ring[samples->face[i].a].push_back(i);
		vertex_faces_1ring[samples->face[i].b].push_back(i);
		vertex_faces_1ring[samples->face[i].c].push_back(i);
	}
	
	// find 1-ring neighboring vertices
	for (int i = 0; i < samples->nv; ++i)
	{
		for (size_t j = 0; j < vertex_faces_1ring[i].size(); ++j)
		{
			neigh_vertex_1ring[i].insert(samples->face[vertex_faces_1ring[i][j]].a);
			neigh_vertex_1ring[i].insert(samples->face[vertex_faces_1ring[i][j]].b);
			neigh_vertex_1ring[i].insert(samples->face[vertex_faces_1ring[i][j]].c);
		}
	}

#ifdef __TEST__
	ofstream ofs((ROOT + string("DEBUG.vertex_faces_1ring.txt")).c_str());
	
	if (! ofs.is_open())
	{
		cerr << "ERROR: open debug_vertex_faces_1ring.txt failed ..." << endl;
		exit(EXIT_FAILURE);
	}
	
	for (int i = 0; i < vertex_faces_1ring.size(); ++i)
	{
		ofs << "[v <-> f] [" << vertex_faces_1ring[i].size() << "] v" << i << ": ";
		for (int j = 0; j < vertex_faces_1ring[i].size(); ++j)
		{
			ofs << vertex_faces_1ring[i][j] << " ";
		}
		ofs << endl;
	}
	
	ofs.close();
#endif
	
#ifdef __TEST__
	ofstream ofs2((ROOT + string("DEBUG.neigh_vertex_1ring.txt")).c_str());
	
	if (! ofs2.is_open())
	{
		cerr << "ERROR: open debug_neigh_vertex_1ring.txt failed ..." << endl;
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < neigh_vertex_1ring.size(); ++i)
	{
		ofs2 << "[v <-> v] [" << neigh_vertex_1ring[i].size() << "] v" << i << ": ";
		for (set<int>::const_iterator it = neigh_vertex_1ring[i].begin();
			 it != neigh_vertex_1ring[i].end(); ++it)
		{
			ofs2 << *it << " ";
		}
		ofs2 << endl;
	}
	
	ofs2.close();
#endif
}

//void compute_rho()
//{
//	rho = (double*) malloc(sizeof(double) * samples->nf);
//	if (! rho)
//	{
//		cerr << "ERROR: allocate rho failed." << endl;
//	}
//	
//	double *dist = NULL;
//	for (int f = 0; f < samples->nf; ++f)
//	{
//		int counter = 0;
//#ifdef NEIGHBOR_1RING
//		set<int> &nb = face_faces_1ring[f]; // neighbors of centers[f]
//#endif
//#ifdef NEIGHBOR_2RING
//		set<int> &nb = face_faces_2ring[f]; // neighbors of centers[f]
//#endif
//		size_t sz = nb.size();
//		dist = (double*) malloc(sizeof(double) * sz);
//		for (set<int>::const_iterator it = nb.begin(); it != nb.end(); ++it)
//		{
//			DBL7VECT &v0 = centers[f];
//			DBL7VECT &v1 = centers[*it];
//			dist[counter++] = compute_distance(v0, v1, T[f]);
//		}
//		rho[f] = *max_element(dist, dist + sz);
//		free(dist);
//	}
//	
//#ifdef __TEST__
//	ofstream ofs(ROOT + "DEBUG.rho.txt");
//	for (int f = 0; f < samples->nf; ++f)
//	{
//		ofs << "c" << f << ": " << rho[f] << "\n";
//	}
//	ofs << endl;
//	ofs.close();
//#endif
//}

#ifdef WEIGHTED
void weight_interpolate()
{
    g_F = allocate_2D<double>(X, Y);
	for (int j = 0; j < Y; ++j)
	{
		memset(g_F[j], 0.0, sizeof(double) * X);
	}
    
    for (int i = 0; i < samples->nf; ++i)
	{
#if defined(NEIGHBOR_1RING)
		set<int> &neighborhood = face_faces_1ring[i]; // neighboring vertices
#elif defined(NEIGHBOR_2RING)
		set<int> &neighborhood = face_faces_2ring[i]; // neighboring vertices
#endif
		size_t num_neig = neighborhood.size();
		MatrixXd ma(num_neig, num_neig); // distance matrix
		VectorXd u(num_neig); // intensities
		VectorXd coeff(num_neig); // solution
		
		// assembly distance matrix and vector
		size_t counter1 = 0, counter2 = 0;
		for (set<int>::const_iterator it = neighborhood.begin();
			 it != neighborhood.end(); ++it)
		{
			u(counter1) = centers[*it].intensity;
			counter2 = 0;
			for (set<int>::const_iterator it2 = neighborhood.begin();
				 it2 != neighborhood.end(); ++it2)
			{
				double r = compute_distance(centers[*it], centers[*it2], T[*it2]);
				ma(counter1, counter2++) = phi(r);
			}
			counter1 ++;
		}
		
		// solve interpolation coefficents for each local domain
		coeff = ma.fullPivLu().solve(u);
        
        for (int m = 0; m < N1; ++m)
		{
			for (int n = 0; n < N2; ++n)
			{
				INT2VECT p;
				p.x = m; p.y = n;
                
                if (i == get<0>(g_pixel_tag[m][n]))
				{
                    // if pixel[m][n] is in current triangle i
                    
                    DBL7VECT pp[3] = {
                        samples->vertex[samples->face[i].a],
                        samples->vertex[samples->face[i].b],
                        samples->vertex[samples->face[i].c]
                    };
					
					// interpolate intensities for 3 vertices of trianlge i
					double intensities[3] = {0.0, 0.0, 0.0};
                    
                    counter1 = 0;
					for (set<int>::const_iterator it = neighborhood.begin();
						 it != neighborhood.end(); ++it)
					{
						double r = compute_distance(pp[0], centers[*it], T[*it]);
						
						// compute interpolated intensities
						intensities[0] += (coeff(counter1++) * phi(r));
					}
                    
                    counter1 = 0;
                    for (set<int>::const_iterator it = neighborhood.begin();
						 it != neighborhood.end(); ++it)
					{
						double r = compute_distance(pp[1], centers[*it], T[*it]);
						
						// compute interpolated intensities
						intensities[1] += (coeff(counter1++) * phi(r));
					}
                    
                    counter1 = 0;
                    for (set<int>::const_iterator it = neighborhood.begin();
						 it != neighborhood.end(); ++it)
					{
						double r = compute_distance(pp[2], centers[*it], T[*it]);
						
						// compute interpolated intensities
						intensities[2] += (coeff(counter1++) * phi(r));
					}
                    
                    // transform eigenvalues for 3 vertices of triangle i
                    double eig2[3] = {0.0, 0.0, 0.0};
                    Matrix2d t1[3], t2[3], T[3];
                    double weights[3];
                    
                    double c = compute_c(pp[0].lambda1, pp[0].lambda2);
                    if (c != 0.0)
                    {
                        eig2[0] = 1.0 - expf(-3.8 * SQ(C0/c));
                    }
                    else
                    {
                        eig2[0] = 1.0;
                    }
#ifdef ISO_RBF
                    eig2[0] = 1.0;
#endif
                    t1[0](0, 0) = pp[0].v1.x;
                    t1[0](1, 0) = pp[0].v1.y;
                    t1[0](0, 1) = pp[0].v2.x;
                    t1[0](1, 1) = pp[0].v2.y;
                    t2[0](0, 0) = 1.0;
                    t2[0](0, 1) = 0.0;
                    t2[0](1, 0) = 0.0;
                    t2[0](1, 1) = eig2[0];
                    T[0] = t1[0] * t2[0];
                    T[0] = T[0] * t1[0].transpose();
                    
                    c = compute_c(pp[1].lambda1, pp[1].lambda2);
                    if (c != 0.0)
                    {
                        eig2[1] = 1.0 - expf(-3.8 * SQ(C0/c));
                    }
                    else
                    {
                        eig2[1] = 1.0;
                    }
#ifdef ISO_RBF
                    eig2[1] = 1.0;
#endif
                    t1[1](0, 0) = pp[1].v1.x;
                    t1[1](1, 0) = pp[1].v1.y;
                    t1[1](0, 1) = pp[1].v2.x;
                    t1[1](1, 1) = pp[1].v2.y;
                    t2[1](0, 0) = 1.0;
                    t2[1](0, 1) = 0.0;
                    t2[1](1, 0) = 0.0;
                    t2[1](1, 1) = eig2[1];
                    T[1] = t1[1] * t2[1];
                    T[1] = T[1] * t1[1].transpose();
                    
                    c = compute_c(pp[2].lambda1, pp[2].lambda2);
                    if (c != 0.0)
                    {
                        eig2[2] = 1.0 - expf(-3.8 * SQ(C0/c));
                    }
                    else
                    {
                        eig2[2] = 1.0;
                    }
#ifdef ISO_RBF
                    eig2[2] = 1.0;
#endif
                    t1[2](0, 0) = pp[2].v1.x;
                    t1[2](1, 0) = pp[2].v1.y;
                    t1[2](0, 1) = pp[2].v2.x;
                    t1[2](1, 1) = pp[2].v2.y;
                    t2[2](0, 0) = 1.0;
                    t2[2](0, 1) = 0.0;
                    t2[2](1, 0) = 0.0;
                    t2[2](1, 1) = eig2[2];
                    T[2] = t1[2] * t2[2];
                    T[2] = T[2] * t1[2].transpose();
                    
                    // compute weights
                    DBL7VECT p; p.x = m; p.y = n;
                    
                    double r = compute_distance(p, pp[0], T[0]);
                    weights[0] = compute_weight(r);
                    r = compute_distance(p, pp[1], T[1]);
                    weights[1] = compute_weight(r);
                    r = compute_distance(p, pp[2], T[2]);
                    weights[2] = compute_weight(r);
					
                    double sum_weights = weights[0] + weights[1] + weights[2];
                    
					// sum interpolated intensity
                    g_F[m][n] = (weights[0] * intensities[0] + weights[1] * intensities[1] + weights[2] * intensities[2]) / sum_weights;
				} // end if
			} // end for n
		} // end for m
    }
}
#endif

void interpolate()
{
	g_F = (double*) malloc(sizeof(double) * X * Y);
    memset((void*) g_F, 0, sizeof(double) * X * Y);
	
	for (int f = 0; f < samples->nf; ++f)
	{
        // find TRIANGLE <-> PIXEL mapping
        PixelTag tag;
        double minX, maxX, minY, maxY;
        DBL7VECT &v0 = samples->vertex[samples->face[f].a];
        DBL7VECT &v1 = samples->vertex[samples->face[f].b];
        DBL7VECT &v2 = samples->vertex[samples->face[f].c];
        minX = min(v0.x, v1.x);
        minX = min(minX, v2.x);
        maxX = max(v0.x, v1.x);
        maxX = max(maxX, v2.x);
        minY = min(v0.y, v1.y);
        minY = min(minY, v2.y);
        maxY = max(v0.y, v1.y);
        maxY = max(maxY, v2.y);
        for (int j = (int) minY; j <= (int) ceil(maxY); ++j)
        {
            for (int i = (int) minX; i <= (int) ceil(maxX); ++i)
            {
                int p[2] = {i, j};
                if (is_in_triangle(p, f))
                {
                    tag.f = f;
                    tag.x.push_back(i);
                    tag.y.push_back(j);
                }
            }
        }
        
        // Solve linear equations
#if defined(NEIGHBOR_1RING)
		set<int> &neighborhood = face_faces_1ring[f]; // neighboring vertices
#elif defined(NEIGHBOR_2RING)
		set<int> &neighborhood = face_faces_2ring[f]; // neighboring vertices
#endif
		size_t num_neig = neighborhood.size();
		MatrixXd ma(num_neig, num_neig); // distance matrix
		VectorXd u(num_neig); // intensities
		VectorXd coeff(num_neig); // solution
		// assembly distance matrix and vector
		size_t counter1 = 0, counter2 = 0;
		for (set<int>::const_iterator it = neighborhood.begin();
             it != neighborhood.end(); ++it)
		{
			u(counter1) = centers[*it].intensity;
			counter2 = 0;
			for (set<int>::const_iterator it2 = neighborhood.begin();
				 it2 != neighborhood.end(); ++it2)
			{
                double p1[2] = {centers[*it].x, centers[*it].y};
                double p2[2] = {centers[*it2].x, centers[*it2].y};
				double r = compute_distance(p1, p2, T[*it2]);
				ma(counter1, counter2++) = phi(r);
			}
			counter1 ++;
		}
		coeff = ma.fullPivLu().solve(u); // solve interpolation coefficents
		
		// interpolation
        assert(tag.x.size() == tag.y.size());
        for (size_t k = 0; k < tag.x.size(); ++k)
        {
            int m = tag.x[k];
            int n = tag.y[k];
            double p1[2] = {(double) m, (double) n};
            
            // compute interpolated intensity
            counter1 = 0;
            double intensity = 0.0;
            for (set<int>::const_iterator it = neighborhood.begin();
                 it != neighborhood.end(); ++it)
            {
                double p2[2] = {centers[*it].x, centers[*it].y};
                double r = compute_distance(p1, p2, T[*it]);
                
                // compute interpolated intensities
                intensity += (coeff(counter1++) * phi(r));
            }
            
            g_F[n*X+m] = intensity; // sum interpolated intensity
        } // end for k
    } // end for f in triangles

	// take average and set threshold
    for (int j = 0; j < Y; ++j)
    {
        for (int i = 0; i < X; ++i)
        {
            if ((g_F[j*X+i]) < 0.0)
                g_F[j*X+i] = 0.0;
            if ((g_F[j*X+i]) > 255.0)
                g_F[j*X+i] = 255.0;
        }
    }
}

void read_image(unsigned char** img, const char* file_name, int *rows, int *cols)
{
	const int SZ = 1024;
	int c;
	int i, j;
	int xdim = 0, ydim = 0;
	int maxraw;
	unsigned char *line;
	char buf[SZ];
	unsigned char *temp = NULL;
	
	printf("Reading PPM image \"%s\"\n", file_name);
	FILE* fp = NULL;
	
#ifdef _WIN32
	if (fopen_s(&fp, file_name, "rb"))
#else
	if ((fp = fopen(file_name, "rb")) == NULL)
#endif
	{
		fprintf(stderr, "ERROR: read image ...\n");
		exit(EXIT_FAILURE);
	}
	
	while ((c = fgetc(fp)) == '#')
		fgets(buf, SZ, fp);
	ungetc(c, fp);
#ifdef _WIN32
	if (fscanf_s(fp, "P%d\n", &c) != 1) {
#else
	if (fscanf(fp, "P%d\n", &c) != 1) {
#endif
		fprintf (stderr, "ERROR: read image ...\n");
		exit(EXIT_FAILURE);
	}
	if (c != 6 && c != 2) {
		fprintf(stderr, "ERROR: read image ...\n");
		exit(EXIT_FAILURE);
	}
	
	if (c == 6)
	{
		while ((c = fgetc(fp)) == '#')
			fgets(buf, SZ, fp);
		ungetc(c, fp);
#ifdef _WIN32
		if (fscanf_s(fp, "%d%d%d", &xdim, &ydim, &maxraw) != 3) {
#else
		if (fscanf(fp, "%d%d%d", &xdim, &ydim, &maxraw) != 3) {
#endif
			fprintf(stderr, "ERROR: failed to read width/height/max\n");
			exit(EXIT_FAILURE);
		}
		printf("\tWidth = %d, Height = %d, Maximum = %d ", xdim, ydim, maxraw);
		
		temp = (unsigned char*) malloc(sizeof(unsigned char) * xdim * ydim);
		getc(fp);
		line = (unsigned char*) malloc(sizeof(unsigned char) * xdim * 3);
		for (j = 0; j < ydim; j++)
		{
			fread(line, 1, xdim*3, fp);
			for (i = 0; i < xdim*3; i += 3)
			{
				temp[j * xdim + i / 3] = line[i];
			}
		}
		free(line);
	}
	
	*cols = xdim;
	*rows = ydim;
	*img = temp;
	
	fclose(fp);
	printf("\tdone!\n");
}

#ifdef __TEST__
void write(const char* file_name)
{
	FILE* fp = NULL;
	
	if ((fp = fopen(file_name, "wb")) == NULL)
	{
		fprintf(stderr, "ERROR: write image ...\n");
		exit(EXIT_FAILURE);
	}
	
	/* write to a PPM image */
	fprintf(fp, "P6\n%d %d\n%d\n", X, Y, 255);
	
	for (int j = Y-1; j >= 0; j--)
	{
		for (int i = 0; i < X; i++)
		{
			bool isCenter = false;
			for (int f = 0; f < samples->nf; ++f)
			{
				if (roundf(centers[f].x) == i && roundf(centers[f].y) == j)
				{
					fputc((int) centers[f].intensity, fp);
					fputc((int) centers[f].intensity, fp);
					fputc((int) centers[f].intensity, fp);
					isCenter = true;
					break;
				}
			}
			
			if (! isCenter)
			{
				fputc(0, fp);
				fputc(0, fp);
				fputc(0, fp);
			}
		}
	}
	
	fclose(fp);
}
#endif

void write_image(const char* file_name)
{  
	FILE* fp = NULL;
	
#ifdef _WIN32
	if (fopen_s(&fp, file_name, "wb"))
#else
	if ((fp = fopen(file_name, "wb")) == NULL)
#endif
	{
		fprintf(stderr, "ERROR: write image ...\n");
		exit(EXIT_FAILURE);
	}
	
	printf("Writing data ... ");
	
	/* write to a PPM image */
	fprintf(fp, "P6\n%d %d\n%d\n", X, Y, 255);
	
    for (int j = Y-1; j >= 0; j--)
	{
		for (int i = 0; i < X; i++)
		{
//            if (abs(orig_img[j*N1+i] - g_F[i][j]) >= 60)
//            {
//                fputc(255, fp);
//                fputc(0, fp);
//                fputc(0, fp);
//            }
//            else if (abs(orig_img[j*N1+i] - g_F[i][j]) >= 40)
//            {
//                fputc(0, fp);
//                fputc(255, fp);
//                fputc(0, fp);
//            }
//            else if (abs(orig_img[j*N1+i] - g_F[i][j]) >= 20)
//            {
//                fputc(0, fp);
//                fputc(0, fp);
//                fputc(255, fp);
//            }
//            else
//            {
                fputc((int) g_F[j*X+i], fp);
                fputc((int) g_F[j*X+i], fp);
                fputc((int) g_F[j*X+i], fp);
//            }
		}
	}
	
	fclose(fp);
	printf("done!\n");
}

void rescale_image()
{
	double maxI = -999999.9;
	double minI = 999999.9;
	
	for (int i = 0; i < X; ++i)
	{
		for (int j = 0; j < Y; ++j)
		{
			if (g_F[j*X+i] > maxI)
			{
				maxI = g_F[j*X+i];
			}
			
			if (g_F[j*X+i] < minI)
			{
				minI = g_F[j*X+i];
			}
		}
	}
    
    cout << "minI = " << minI << ", maxI = " << maxI << endl;
	
	for (int i = 0; i < X; ++i)
	{
		for (int j = 0; j < Y; ++j)
		{
			g_F[j*X+i] = (g_F[j*X+i] - minI) / (maxI - minI) * 255.0;
		}
	}
}

void initialize_centers()
{
	centers = (DBL7VECT*) malloc(sizeof(DBL7VECT) * samples->nf);
	if (! centers)
	{
		cerr << "ERROR: allocate memory for triangle centers failed. " << endl;
	}
	
	for (int f = 0; f < samples->nf; ++f)
	{
		int a = samples->face[f].a;
		int b = samples->face[f].b;
		int c = samples->face[f].c;
		
		// compute coordinates
		centers[f].x = (samples->vertex[a].x + samples->vertex[b].x + samples->vertex[c].x) / 3.0;
		centers[f].y = (samples->vertex[a].y + samples->vertex[b].y + samples->vertex[c].y) / 3.0;
		
		// compute intensity
		centers[f].intensity = samples->face[f].intensity;
		
		// compute eigenvalues and eigenvectors
		// use double-precision to avoid precision loss
		Vector2d a1(samples->vertex[a].v1.x, samples->vertex[a].v1.y);
		Vector2d a2(samples->vertex[a].v2.x, samples->vertex[a].v2.y);
		Vector2d b1(samples->vertex[b].v1.x, samples->vertex[b].v1.y);
		Vector2d b2(samples->vertex[b].v2.x, samples->vertex[b].v2.y);
		Vector2d c1(samples->vertex[c].v1.x, samples->vertex[c].v1.y);
		Vector2d c2(samples->vertex[c].v2.x, samples->vertex[c].v2.y);
		Matrix2d ma = samples->vertex[a].lambda1 * a1 * a1.transpose() + samples->vertex[a].lambda2 * a2 * a2.transpose();
		Matrix2d mb = samples->vertex[b].lambda1 * b1 * b1.transpose() + samples->vertex[b].lambda2 * b2 * b2.transpose();
		Matrix2d mc = samples->vertex[c].lambda1 * c1 * c1.transpose() + samples->vertex[c].lambda2 * c2 * c2.transpose();
		Matrix2d mcenter = (ma + mb + mc);
		
		pair<DBL2VECT, DBL2VECT> eigenvectors = compute_eigen(mcenter, &centers[f].lambda1, &centers[f].lambda2);
		centers[f].v1 = eigenvectors.first;
		centers[f].v2 = eigenvectors.second;
	}
	
	// testing
	
#ifdef __TEST__
	string filename = ROOT + "DEBUG.centers.ppm";
	write(filename.c_str());
	ofstream ofs((ROOT + string("DEBUG.centers.txt")).c_str());
	ofs << "X\tY\tINTENSITY\tLAMBDA1\tLAMBDA2\tV1X\tV1Y\tV2X\tV2Y\n";
	for (int f = 0; f < samples->nf; ++f)
	{
		ofs << centers[f].x << "\t" << centers[f].y << "\t" << centers[f].intensity << "\t"
			<< centers[f].lambda1 << "\t" << centers[f].lambda2 << "\t"
			<< centers[f].v1.x << "\t" << centers[f].v1.y << "\t"
			<< centers[f].v2.x << "\t" << centers[f].v2.y << "\n";
	}
	ofs.close();
#endif
}

int main(int argc, const char * argv[])
{
    if (argc != 4)
        cerr << "WRONG: incorrect command format. Correct is :" << endl
        << "\t arbf-image <filename_prefix> <width> <height>\n";
    
    clock_t start, span, total;
    
    // read mesh
    start = total = clock();
    string file = argv[1];
#ifdef _WIN32
	string filename = ROOT + "data\\" + file + ".tri";
#else
	string filename = ROOT + "data/" + file + ".tri";
#endif
    read_mesh_tri(filename.c_str(), &samples);
	cout << "nv = " << samples->nv << ", nf = " << samples->nf << endl;
    span = clock() - start;
    total += span;
    printf("Reading mesh ... %f ms\n",
           (double) span * 1e3 / CLOCKS_PER_SEC);
    
    // read original image
    start = clock();
#ifdef _WIN32
    filename = ROOT + "data\\" + file + ".ppm";
#else
    filename = ROOT + "data/" + file + ".ppm";
#endif
    int x = 0, y = 0;
    read_image(&orig_img, filename.c_str(), &y, &x);
    cout << "rows = " << y << ", cols = " << x << endl;
    span = clock() - start;
    total += span;
    printf("Reading original image ... %f ms\n",
           (double) span * 1e3 / CLOCKS_PER_SEC);
    
    X = atoi(argv[2]);
    Y = atoi(argv[3]);
	
	// initialize triangle centers for each face
    start = clock();
	initialize_centers();
    span = clock() - start;
    total += span;
    printf("Init centers ... %f ms\n",
           (double) span * 1e3 / CLOCKS_PER_SEC);
	
	// find neighboring vertices for each vertex and
    // neighboring faces for each face (center)
    start = clock();
    find_neigh_vertex_1ring();
    find_neigh_face_1ring();
    span = clock() - start;
    total += span;
    printf("Finding 1-Ring neighbors ... %f ms\n",
           (double) span * 1e3 / CLOCKS_PER_SEC);
#ifdef NEIGHBOR_2RING
    start = clock();
	find_neigh_vertex_2ring();
	find_neigh_face_2ring();
    span = clock() - start;
    total += span;
    printf("Finding 2-Ring neighbors ... %f ms\n",
           (double) span * 1e3 / CLOCKS_PER_SEC);
#endif

    // transform eigenvalues
    start = clock();
	transform_eigenvalues();
    span = clock() - start;
    total += span;
    printf("Transforming eigenvalues ... %f ms\n",
           (double) span * 1e3 / CLOCKS_PER_SEC);
	
	// compute metrics
    start = clock();
	compute_metrics();
    span = clock() - start;
    total += span;
    printf("Computing tensor T ... %f ms\n",
           (double) span * 1e3 / CLOCKS_PER_SEC);
	
	// interpolate
#ifdef WEIGHTED
    start = clock();
    weight_interpolate();
    span = clock() - start;
    total += span;
    printf("Weighted interpolating ... %f ms\n",
           (double) span * 1e3 / CLOCKS_PER_SEC);
#else
    start = clock();
	interpolate();
    span = clock() - start;
    total += span;
    printf("Interpolating ... %f ms\n",
           (double) span * 1e3 / CLOCKS_PER_SEC);
#endif
    
    // compute PSNR
    start = clock();
    cout << "PSNR = " << compute_PSNR() << endl;
    span = clock() - start;
    total += span;
    printf("Computing PSNR ... %f ms\n",
           (double) span * 1e3 / CLOCKS_PER_SEC);
    
	// output
    start = clock();
	string outputName = ROOT + "out.ppm";
	write_image(outputName.c_str());
    span = clock() - start;
    total += span;
    printf("Printing result ... %f ms\n",
           (double) span * 1e3 / CLOCKS_PER_SEC);
	
	// clean up
    start = clock();
	clean_up();
    span = clock() - start;
    total += span;
    printf("Cleaning ... %f ms\n",
           (double) span * 1e3 / CLOCKS_PER_SEC);
    
    cout << "Running time = " << (double) total / CLOCKS_PER_SEC
        << " seconds." << endl;
	
	return 0;
}
