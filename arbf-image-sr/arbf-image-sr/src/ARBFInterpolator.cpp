//
//  ARBFInterpolator.cpp
//  arbf-image-super_resolution
//
//  Created by Ke Liu on 2/27/15.
//  Copyright (c) 2015 Ke Liu. All rights reserved.
//

#include <cstdio>
#include <cfenv> // for std::fesetround()
#include <cassert>
//#include <glog/logging.h>
//#include <boost/format.hpp>
#include "../include/Config.h"
#include "../include/ARBFInterpolator.h"

using namespace std;
//using namespace boost;
using namespace Eigen;

DBL7VECT* ARBFInterpolator::s_centers = nullptr;
const PPMImage* ARBFInterpolator::s_imgOrig = nullptr;
PPMImage* ARBFInterpolator::s_imgSR = nullptr;
const TriMesh* ARBFInterpolator::s_mesh = nullptr;
BasisType ARBFInterpolator::s_basis = BasisType::MQ;

vector<pair<double, double> > ARBFInterpolator::s_eigenvalues;
vector<vector<int> > ARBFInterpolator::s_vf_1ring;
vector<set<int> > ARBFInterpolator::s_vv_1ring;
vector<set<int> > ARBFInterpolator::s_ff_1ring;
vector<Matrix2d> ARBFInterpolator::s_T;
#ifdef NEIGHBOR_2RING
vector<vector<int> > ARBFInterpolator::s_vf_2ring;
vector<set<int> > ARBFInterpolator::s_vv_2ring;
vector<set<int> > ARBFInterpolator::s_ff_2ring;
#endif

ARBFInterpolator::ARBFInterpolator() {}

inline const DBL7VECT* ARBFInterpolator::getCenters()
{
    return s_centers;
}

void ARBFInterpolator::setOriginalImage(const PPMImage *img)
{
    s_imgOrig = img;
}

void ARBFInterpolator::setSRImage(PPMImage *img)
{
    s_imgSR = img;
}

void ARBFInterpolator::setMesh(const TriMesh *mesh)
{
    s_mesh = mesh;
}

inline vector<vector<int> >& ARBFInterpolator::getNeigborVertexFace1Ring()
{
    return s_vf_1ring;
}

void ARBFInterpolator::setNeighborVertexFace1Ring(const std::vector<std::vector<int> > &vf)
{
    s_vf_1ring = vf;
}

inline vector<set<int> >& ARBFInterpolator::getNeighborVertexVertex1Ring()
{
    return s_vv_1ring;
}

void ARBFInterpolator::setNeighborVertexVertex1Ring(const std::vector<std::set<int> > &vv)
{
    s_vv_1ring = vv;
}

inline vector<set<int> >& ARBFInterpolator::getNeighborFaceFace1Ring()
{
    return s_ff_1ring;
}

inline void ARBFInterpolator::setNeighborFaceFace1Ring(std::vector<std::set<int> > &ff)
{
    s_ff_1ring = ff;
}

void ARBFInterpolator::setBasisType(const BasisType type)
{
	s_basis = type;
}

#ifdef NEIGHBOR_2RING
inline vector<vector<int> >& ARBFInterpolator::getNeigborVertexFace2Ring()
{
    return s_vf_2ring;
}

void ARBFInterpolator::setNeighborVertexFace2Ring(const std::vector<std::vector<int> > &vf)
{
    s_vf_2ring = vf;
}

inline vector<set<int> >& ARBFInterpolator::getNeighborVertexVertex2Ring()
{
    return s_vv_2ring;
}

void ARBFInterpolator::setNeighborVertexVertex2Ring(const std::vector<std::set<int> > &vv)
{
    s_vv_2ring = vv;
}

inline vector<set<int> >& ARBFInterpolator::getNeighborFaceFace2Ring()
{
    return s_ff_2ring;
}

void ARBFInterpolator::setNeighborFaceFace2Ring(std::vector<std::set<int> > &ff)
{
    s_ff_2ring = ff;
}

int ARBFInterpolator::findNeighVertices2Ring()
{
    s_vv_2ring.resize(s_mesh->numVertices());
    s_vf_2ring.resize(s_mesh->numVertices());
    
    // build 2-ring VERTEX <-> FACE mapping
    s_vf_2ring.assign(s_vf_1ring.begin(), s_vf_1ring.end());
    for (int i = 0; i < s_mesh->numFaces(); ++i) {
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
    for (int i = 0; i < s_mesh->numVertices(); ++i) {
        for (int j = 0; j < s_vf_2ring[i].size(); ++j) {
            s_vv_2ring[i].insert(s_mesh->getFaces()[s_vf_2ring[i][j]].a);
            s_vv_2ring[i].insert(s_mesh->getFaces()[s_vf_2ring[i][j]].b);
            s_vv_2ring[i].insert(s_mesh->getFaces()[s_vf_2ring[i][j]].c);
        }
    }
    
#ifdef __TEST__
    string filename = "DEBUG.vertex_faces_2ring.txt";
    FILE *fp = nullptr;

#ifdef _WIN32
    if (fopen_s(&fp, filename.c_str(), "w")) {
#else
    if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
#endif
        fprintf(stderr, "ERROR: open DEBUG.vertex_faces_2ring.txt failed ...\n");
        
        return 0;
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
    
    filename = "DEBUG.neigh_vertex_2ring.txt";
    
#ifdef _WIN32
        if (fopen_s(&fp, filename.c_str(), "w")) {
#else
        if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
#endif
        fprintf(stderr, "ERROR: open DEBUG.neigh_vertex_2ring.txt failed ...");
        
        return 0;
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
#endif
    
    return 1;
}

int ARBFInterpolator::findNeighFaces2Ring()
{
    s_ff_2ring.resize(s_mesh->numFaces());
    
    for (int f = 0; f < s_mesh->numFaces(); ++f) {
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
    
#ifdef __TEST__
    string filename = "DEBUG.face_faces_2ring.txt";
    FILE *fp = nullptr;
    
#ifdef _WIN32
    if (fopen_s(&fp, filename.c_str(), "w")) {
#else
    if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
#endif
        fprintf(stderr, "ERROR: open DEBUG.face_faces_2ring.txt failed ...");
        
        return 0;
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
#endif
    
    return 1;
}
#endif

int ARBFInterpolator::calculateCenters()
{
    s_centers = new DBL7VECT[s_mesh->numFaces()];
    if (s_centers == nullptr) {
        fprintf(stderr, "ERROR: allocate memory for triangle centers failed.\n");
        return 0;
    }
    
    for (int f = 0; f < s_mesh->numFaces(); ++f)
    {
        int a = s_mesh->getFaces()[f].a;
        int b = s_mesh->getFaces()[f].b;
        int c = s_mesh->getFaces()[f].c;
        
        // compute coordinates
        s_centers[f].x = (s_mesh->getVertices()[a].x
                          + s_mesh->getVertices()[b].x
                          + s_mesh->getVertices()[c].x) / 3.0;
        s_centers[f].y = (s_mesh->getVertices()[a].y
                          + s_mesh->getVertices()[b].y
                          + s_mesh->getVertices()[c].y) / 3.0;
        
        // compute intensity
        s_centers[f].intensity = s_mesh->getFaces()[f].intensity;
        
        // compute eigenvalues and eigenvectors
        Vector2d a1(s_mesh->getVertices()[a].v1.x, s_mesh->getVertices()[a].v1.y);
        Vector2d a2(s_mesh->getVertices()[a].v2.x, s_mesh->getVertices()[a].v2.y);
        Vector2d b1(s_mesh->getVertices()[b].v1.x, s_mesh->getVertices()[b].v1.y);
        Vector2d b2(s_mesh->getVertices()[b].v2.x, s_mesh->getVertices()[b].v2.y);
        Vector2d c1(s_mesh->getVertices()[c].v1.x, s_mesh->getVertices()[c].v1.y);
        Vector2d c2(s_mesh->getVertices()[c].v2.x, s_mesh->getVertices()[c].v2.y);
        Matrix2d ma = s_mesh->getVertices()[a].lambda1 * a1 * a1.transpose()
            + s_mesh->getVertices()[a].lambda2 * a2 * a2.transpose();
        Matrix2d mb = s_mesh->getVertices()[b].lambda1 * b1 * b1.transpose()
            + s_mesh->getVertices()[b].lambda2 * b2 * b2.transpose();
        Matrix2d mc = s_mesh->getVertices()[c].lambda1 * c1 * c1.transpose()
            + s_mesh->getVertices()[c].lambda2 * c2 * c2.transpose();
        Matrix2d mcenter = (ma + mb + mc);
        
        pair<DBL2VECT, DBL2VECT> eigenvectors = s_computeEigen(mcenter, &s_centers[f].lambda1, &s_centers[f].lambda2);
        s_centers[f].v1 = eigenvectors.first;
        s_centers[f].v2 = eigenvectors.second;
        s_centers[f].gradient.x = 0;
        s_centers[f].gradient.y = 0;
    }
    
//#ifdef __TEST__
    string filename = "DEBUG.centers.ppm";
    PPMImage imgCenter(s_imgOrig->X(), s_imgOrig->Y());
    imgCenter.setPath(filename.c_str());
    double *data = imgCenter.getImageData();
    data = new double[s_imgOrig->X() * s_imgOrig->Y()];
    for (int j = s_imgOrig->Y() - 1; j >= 0; j--) {
        for (int i = 0; i <= s_imgOrig->X(); i++) {
            for (int f = 0; f < s_mesh->numFaces(); f++) {
                if ((round(s_centers[f].x) == i) && (round(s_centers[f].y) == j))
                    data[j*s_imgOrig->X()+i] = s_centers[f].intensity;
            }
        }
    }
    imgCenter.write();
    filename = "DEBUG.centers.txt";
    FILE *fp = nullptr;
#ifdef _WIN32
    if (fopen_s(&fp, filename.c_str(), "w")) {
#else
    if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
#endif
        fprintf(stderr, "ERROR: write centers ...\n");
        
        return 0;
    }

    fprintf(fp, "X\tY\tINTENSITY\tLAMBDA1\tLAMBDA2\tV1X\tV1Y\tV2X\tV2Y\tGRADX\tGRADY\n");
    for (int f = 0; f < s_mesh->numFaces(); ++f) {
        fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                s_centers[f].x, s_centers[f].y, s_centers[f].intensity,
                s_centers[f].lambda1, s_centers[f].lambda2,
                s_centers[f].v1.x, s_centers[f].v1.y,
                s_centers[f].v2.x, s_centers[f].v2.y,
                s_centers[f].gradient.x, s_centers[f].gradient.y);
    }
    fclose(fp);
//#endif
    
    return 1;
}

void ARBFInterpolator::computeMetrics()
{
    s_T.resize(s_mesh->numFaces());
    
    for (int f = 0; f < s_mesh->numFaces(); ++f) {
        DBL7VECT *center = &s_centers[f];
        //             |lambda1 0      | |v1T|
        // T = [v1 v2] |               | |   |
        //             |0       lambda2| |v2T|
        
        Matrix2d t1, t2;
        t1(0, 0) = center->v1.x;
        t1(1, 0) = center->v1.y;
        t1(0, 1) = center->v2.x;
        t1(1, 1) = center->v2.y;
        
        t2(0, 0) = s_eigenvalues[f].first;
        t2(0, 1) = 0.0;
        t2(1, 0) = 0.0;
        t2(1, 1) = s_eigenvalues[f].second;
        
        s_T[f] = t1 * t2;
        s_T[f] = s_T[f] * t1.transpose();
    }
}

int ARBFInterpolator::interpolate()
{
	int srX = s_imgSR->X();
	int srY = s_imgSR->Y();
	double stepX = 1.0 / s_imgSR->getScaleX();
	double stepY = 1.0 / s_imgSR->getScaleY();
	fesetround(FE_TONEAREST); // set round direction to to-nearest

    for (int f = 0; f < s_mesh->numFaces(); ++f) {
        // find TRIANGLE <-> PIXEL mapping
        vector<double> x;
        vector<double> y;
        double minX, maxX, minY, maxY;
        const DBL7VECT &v0 = s_mesh->getVertices()[s_mesh->getFaces()[f].a];
        const DBL7VECT &v1 = s_mesh->getVertices()[s_mesh->getFaces()[f].b];
        const DBL7VECT &v2 = s_mesh->getVertices()[s_mesh->getFaces()[f].c];
        minX = min(v0.x, v1.x);
        minX = min(minX, v2.x);
        maxX = max(v0.x, v1.x);
        maxX = max(maxX, v2.x);
        minY = min(v0.y, v1.y);
        minY = min(minY, v2.y);
        maxY = max(v0.y, v1.y);
        maxY = max(maxY, v2.y);
		for (double j = minY; j <= ceil(maxY + stepY); j += stepY) {
			for (double i = minX; i <= ceil(maxX + stepX); i += stepX) {
                double p[2] = {i, j};
                if (s_isInTriangle(p, f)) {
                    x.push_back(i);
                    y.push_back(j);
                }
            }
        }
        
        set<int> &neighborhood = (Config::isUsing2Ring) ? s_ff_2ring[f] : s_ff_1ring[f]; // neighboring vertices
        size_t num_neig = neighborhood.size();
//        MatrixXd ma(num_neig+6, num_neig); // distance matrix
//        VectorXd u(num_neig+6); // intensities
        
        MatrixXd ma(num_neig, num_neig); // distance matrix
        VectorXd u(num_neig); // intensities
        
        VectorXd coeff(num_neig); // solution
        
        // assembly distance matrix and vector
        size_t counter1 = 0, counter2 = 0;
		set<int>::const_iterator it;
		set<int>::const_iterator it2;
        for (it = neighborhood.begin(); it != neighborhood.end(); ++it) {
            u(counter1) = s_centers[*it].intensity;
            counter2 = 0;
            for (it2 = neighborhood.begin(); it2 != neighborhood.end(); ++it2) {
                double p1[2] = {s_centers[*it].x, s_centers[*it].y};
                double p2[2] = {s_centers[*it2].x, s_centers[*it2].y};
//                double r = s_computeDistance(p1, p2, s_T[*it2]);
                double r  = sqrt(SQ(p1[0]-p2[0]) + SQ(p1[1]-p2[1]));
                
                ma(counter1, counter2++) = s_phi(r);
            }
            counter1 ++;
        }
        
//        u(num_neig) = v0.gradient.x;
//        u(num_neig+1) = v0.gradient.y;
//        u(num_neig+2) = v1.gradient.x;
//        u(num_neig+3) = v1.gradient.y;
//        u(num_neig+4) = v2.gradient.x;
//        u(num_neig+5) = v2.gradient.y;

//        counter2 = 0;
//        for (it2 = neighborhood.begin(); it2 != neighborhood.end(); ++it2) {
//            double p2[2] = {s_centers[*it2].x, s_centers[*it2].y};
//            
//            ma(num_neig, counter2++) = (v0.x-p2[0]) / pow(sqrt(SQ(v0.x-p2[0])+SQ(v0.y-p2[1])+SQ(0.5)),3);
//        }
//        
//        counter2 = 0;
//        for (it2 = neighborhood.begin(); it2 != neighborhood.end(); ++it2) {
//            double p2[2] = {s_centers[*it2].x, s_centers[*it2].y};
//            
//            ma(num_neig+1, counter2++) = (v0.y-p2[0]) / pow(sqrt(SQ(v0.x-p2[0])+SQ(v0.y-p2[1])+SQ(0.5)),3);
//        }
//        
//        counter2 = 0;
//        for (it2 = neighborhood.begin(); it2 != neighborhood.end(); ++it2) {
//            double p2[2] = {s_centers[*it2].x, s_centers[*it2].y};
//            
//            ma(num_neig+2, counter2++) = (v1.x-p2[0]) / pow(sqrt(SQ(v1.x-p2[0])+SQ(v1.y-p2[1])+SQ(0.5)),3);
//        }
//        
//        counter2 = 0;
//        u(num_neig+1) = v0.gradient.y;
//        for (it2 = neighborhood.begin(); it2 != neighborhood.end(); ++it2) {
//            double p2[2] = {s_centers[*it2].x, s_centers[*it2].y};
//            
//            ma(num_neig+3, counter2++) = (v1.y-p2[0]) / pow(sqrt(SQ(v1.x-p2[0])+SQ(v1.y-p2[1])+SQ(0.5)),3);
//        }
//        
//        counter2 = 0;
//        u(num_neig+1) = v0.gradient.y;
//        for (it2 = neighborhood.begin(); it2 != neighborhood.end(); ++it2) {
//            double p2[2] = {s_centers[*it2].x, s_centers[*it2].y};
//            
//            ma(num_neig+1, counter2++) = (v2.x-p2[0]) / pow(sqrt(SQ(v2.x-p2[0])+SQ(v2.y-p2[1])+SQ(0.5)),3);
//        }
//        
//        counter2 = 0;
//        u(num_neig+1) = v0.gradient.y;
//        for (it2 = neighborhood.begin(); it2 != neighborhood.end(); ++it2) {
//            double p2[2] = {s_centers[*it2].x, s_centers[*it2].y};
//            
//            ma(num_neig+1, counter2++) = (v2.y-p2[0]) / pow(sqrt(SQ(v2.x-p2[0])+SQ(v2.y-p2[1])+SQ(0.5)),3);
//        }
        
        
        // solve interpolation coefficents for each local domain
        coeff = ma.fullPivLu().solve(u);
//        coeff = ma.jacobiSvd(ComputeThinU | ComputeThinV).solve(u);
        
        // interpolation
        assert(x.size() == y.size());
        for (size_t k = 0; k < x.size(); ++k) {
			int m = (int)nearbyint(x[k] * s_imgSR->getScaleX());
			int n = (int)nearbyint(y[k] * s_imgSR->getScaleY());
			double p1[2] = {x[k], y[k]};
            
            // compute interpolated intensity
            counter1 = 0;
            double intensity = 0.0;
            for (it = neighborhood.begin(); it != neighborhood.end(); ++it) {
                double p2[2] = {s_centers[*it].x, s_centers[*it].y};
//                double r = s_computeDistance(p1, p2, s_T[*it]);
                double r  = sqrt(SQ(p1[0]-p2[0]) + SQ(p1[1]-p2[1]));
                
                // compute interpolated intensities
                intensity += (coeff(counter1++) * s_phi(r));
            }
            
            // sum interpolated intensity
            s_imgSR->getImageData()[n*srX+m] = intensity;
        } // end for k
    } // end for f in triangles
    
    // take average and set threshold
    for (int j = 0; j < srY; ++j) {
        for (int i = 0; i < srX; ++i) {
            if ((s_imgSR->getImageData()[j*srX+i]) < 0.0)
                s_imgSR->getImageData()[j*srX+i] = 0.0;
            if ((s_imgSR->getImageData()[j*srX+i]) > s_imgOrig->getMaxIntensity())
                s_imgSR->getImageData()[j*srX+i] = s_imgOrig->getMaxIntensity();
        }
    }
    
    return 1;
}

int ARBFInterpolator::weightedInterpolate()
{
    int srX = s_imgSR->X();
    int srY = s_imgSR->Y();
    double stepX = 1.0 / s_imgSR->getScaleX();
    double stepY = 1.0 / s_imgSR->getScaleY();
    fesetround(FE_TONEAREST); // set round direction to to-nearest
    
    for (int f = 0; f < s_mesh->numFaces(); ++f) {
        // find TRIANGLE <-> PIXEL mapping
        vector<double> x;
        vector<double> y;
        double minX, maxX, minY, maxY;
        const DBL7VECT &v0 = s_mesh->getVertices()[s_mesh->getFaces()[f].a];
        const DBL7VECT &v1 = s_mesh->getVertices()[s_mesh->getFaces()[f].b];
        const DBL7VECT &v2 = s_mesh->getVertices()[s_mesh->getFaces()[f].c];
        minX = min(v0.x, v1.x);
        minX = min(minX, v2.x);
        maxX = max(v0.x, v1.x);
        maxX = max(maxX, v2.x);
        minY = min(v0.y, v1.y);
        minY = min(minY, v2.y);
        maxY = max(v0.y, v1.y);
        maxY = max(maxY, v2.y);
        for (double j = minY; j <= ceil(maxY + stepY); j += stepY) {
            for (double i = minX; i <= ceil(maxX + stepX); i += stepX) {
                double p[2] = {i, j};
                if (s_isInTriangle(p, i)) {
                    x.push_back(i);
                    y.push_back(j);
                }
            }
        }
        
        set<int> &neighborhood = (Config::isUsing2Ring) ? s_ff_2ring[f] : s_ff_1ring[f]; // neighboring vertices
        size_t num_neig = neighborhood.size();
        MatrixXd ma(num_neig, num_neig); // distance matrix
        VectorXd u(num_neig); // intensities
        VectorXd coeff(num_neig); // solution
        
        // assembly distance matrix and vector
        size_t counter1 = 0, counter2 = 0;
        set<int>::const_iterator it;
        set<int>::const_iterator it2;
        for (it = neighborhood.begin(); it != neighborhood.end(); ++it) {
            u(counter1) = s_centers[*it].intensity;
            counter2 = 0;
            for (it2 = neighborhood.begin(); it2 != neighborhood.end(); ++it2) {
                double p1[2] = {s_centers[*it].x, s_centers[*it].y};
                double p2[2] = {s_centers[*it2].x, s_centers[*it2].y};
                double r = s_computeDistance(p1, p2, s_T[*it2]);
                ma(counter1, counter2++) = s_phi(r);
            }
            counter1 ++;
        }
        
        // solve interpolation coefficents for each local domain
        coeff = ma.fullPivLu().solve(u);
        
        assert(x.size() == y.size());
        for (size_t k = 0; k < x.size(); ++k) {
            int m = (int)nearbyint(x[k] * s_imgSR->getScaleX());
            int n = (int)nearbyint(y[k] * s_imgSR->getScaleY());
            double p[2] = {x[k], y[k]};
            const DBL7VECT &v0 = s_mesh->getVertices()[s_mesh->getFaces()[f].a];
            const DBL7VECT &v1 = s_mesh->getVertices()[s_mesh->getFaces()[f].b];
            const DBL7VECT &v2 = s_mesh->getVertices()[s_mesh->getFaces()[f].c];
            double p0[2] = {v0.x, v0.y};
            double p1[2] = {v1.x, v1.y};
            double p2[2] = {v2.x, v2.y};
            
            // interpolate intensities for 3 vertices of trianlge f
            double intensities[3] = {0.0, 0.0, 0.0};
            counter1 = 0;
            for (it = neighborhood.begin(); it != neighborhood.end(); ++it) {
                double pc[2] = {s_centers[*it].x, s_centers[*it].y};
                double r = s_computeDistance(p0, pc, s_T[*it]);
                // compute interpolated intensities
                intensities[0] += (coeff(counter1++) * s_phi(r));
            }
            
            counter1 = 0;
            for (it = neighborhood.begin(); it != neighborhood.end(); ++it) {
                double pc[2] = {s_centers[*it].x, s_centers[*it].y};
                double r = s_computeDistance(p1, pc, s_T[*it]);
                // compute interpolated intensities
                intensities[1] += (coeff(counter1++) * s_phi(r));
            }
            
            counter1 = 0;
            for (it = neighborhood.begin(); it != neighborhood.end(); ++it) {
                double pc[2] = {s_centers[*it].x, s_centers[*it].y};
                double r = s_computeDistance(p2, pc, s_T[*it]);
                
                // compute interpolated intensities
                intensities[2] += (coeff(counter1++) * s_phi(r));
            }
            
            // transform eigenvalues for 3 vertices of triangle i
            double eig2[3] = {0.0, 0.0, 0.0};
            Matrix2d t1[3], t2[3], T[3];
            double weights[3];
            
            if (Config::isUsingISORBF) {
                eig2[0] = 1.0;
            } else {
                double c = s_computeC(v0.lambda1, v0.lambda2);
                if (c != 0.0) {
                    eig2[0] = 1.0 - exp(-3.8 * SQ(Config::C0 / c));
                } else {
                    eig2[0] = 1.0;
                }
            }
            
            t1[0](0, 0) = v0.v1.x;
            t1[0](1, 0) = v0.v1.y;
            t1[0](0, 1) = v0.v2.x;
            t1[0](1, 1) = v0.v2.y;
            t2[0](0, 0) = 1.0;
            t2[0](0, 1) = 0.0;
            t2[0](1, 0) = 0.0;
            t2[0](1, 1) = eig2[0];
            T[0] = t1[0] * t2[0];
            T[0] = T[0] * t1[0].transpose();
            
            if (Config::isUsingISORBF) {
                eig2[1] = 1.0;
            } else {
                double c = s_computeC(v1.lambda1, v1.lambda2);
                if (c != 0.0) {
                    eig2[1] = 1.0 - exp(-3.8 * SQ(Config::C0 / c));
                } else {
                    eig2[1] = 1.0;
                }
            }
            
            t1[1](0, 0) = v1.v1.x;
            t1[1](1, 0) = v1.v1.y;
            t1[1](0, 1) = v1.v2.x;
            t1[1](1, 1) = v1.v2.y;
            t2[1](0, 0) = 1.0;
            t2[1](0, 1) = 0.0;
            t2[1](1, 0) = 0.0;
            t2[1](1, 1) = eig2[1];
            T[1] = t1[1] * t2[1];
            T[1] = T[1] * t1[1].transpose();
            
            if (Config::isUsingISORBF) {
                eig2[2] = 1.0;
            } else {
                double c = s_computeC(v2.lambda1, v2.lambda2);
                if (c != 0.0) {
                    eig2[2] = 1.0 - expf(-3.8 * SQ(Config::C0 / c));
                } else {
                    eig2[2] = 1.0;
                }
            }

            t1[2](0, 0) = v2.v1.x;
            t1[2](1, 0) = v2.v1.y;
            t1[2](0, 1) = v2.v2.x;
            t1[2](1, 1) = v2.v2.y;
            t2[2](0, 0) = 1.0;
            t2[2](0, 1) = 0.0;
            t2[2](1, 0) = 0.0;
            t2[2](1, 1) = eig2[2];
            T[2] = t1[2] * t2[2];
            T[2] = T[2] * t1[2].transpose();
            
            // compute weights
            double r = s_computeDistance(p, p0, s_T[0]);
            weights[0] = s_computeWeight(r);
            r = s_computeDistance(p, p1, s_T[1]);
            weights[1] = s_computeWeight(r);
            r = s_computeDistance(p, p2, s_T[2]);
            weights[2] = s_computeWeight(r);
            
            double sum_weights = weights[0] + weights[1] + weights[2];
            
            // sum interpolated intensity
            s_imgSR->getImageData()[n*srX+m] = (weights[0] * intensities[0] + weights[1] * intensities[1] + weights[2] * intensities[2]) / sum_weights;
        } // end for k
    } // end for f
    
    // take average and set threshold
    for (int j = 0; j < srY; ++j) {
        for (int i = 0; i < srX; ++i) {
            if ((s_imgSR->getImageData()[j*srX+i]) < 0.0)
                s_imgSR->getImageData()[j*srX+i] = 0.0;
            if ((s_imgSR->getImageData()[j*srX+i]) > s_imgOrig->getMaxIntensity())
                s_imgSR->getImageData()[j*srX+i] = s_imgOrig->getMaxIntensity();
        }
    }
    
    return 1;
}
    
double ARBFInterpolator::computePSNR()
{
    return 20.0 * log10(255.0 / s_computeRMSE());
}

int ARBFInterpolator::findNeighVertices1Ring()
{
    s_vv_1ring.resize(s_mesh->numVertices());
    s_vf_1ring.resize(s_mesh->numVertices());
    
    // build 1-ring VERTEX <-> FACE mapping
    for (int f = 0; f < s_mesh->numFaces(); ++f) {
        s_vf_1ring[s_mesh->getFaces()[f].a].push_back(f);
        s_vf_1ring[s_mesh->getFaces()[f].b].push_back(f);
        s_vf_1ring[s_mesh->getFaces()[f].c].push_back(f);
    }
    
    // find 1-ring VERTEX <-> VERTEX mapping
    for (int v = 0; v < s_mesh->numVertices(); ++v) {
        for (size_t j = 0; j < s_vf_1ring[v].size(); ++j) {
            s_vv_1ring[v].insert(s_mesh->getFaces()[s_vf_1ring[v][j]].a);
            s_vv_1ring[v].insert(s_mesh->getFaces()[s_vf_1ring[v][j]].b);
            s_vv_1ring[v].insert(s_mesh->getFaces()[s_vf_1ring[v][j]].c);
        }
    }
    
#ifdef __TEST__
    string filename = string("DEBUG.vertex_faces_1ring.txt");
    FILE *fp = nullptr;
#ifdef _WIN32
    if (fopen_s(&fp, filename.c_str(), "w")) {
#else
    if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
#endif
        fprintf(stderr, "ERROR: open debug_vertex_faces_1ring.txt failed ...\n");
        
        return 0;
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
        
    filename = string("DEBUG.neigh_vertex_1ring.txt");
#ifdef _WIN32
    if (fopen_s(&fp, filename.c_str(), "w")) {
#else
    if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
#endif
        fprintf(stderr, "ERROR: open debug_neigh_vertex_1ring.txt failed ...\n");
        
        return 0;
    }
    
    for (int i = 0; i < s_vv_1ring.size(); ++i)
    {
        fprintf(fp, "[v <-> v] [%ld] v%d: ", s_vv_1ring[i].size(), i);
        for (set<int>::const_iterator it = s_vv_1ring[i].begin();
             it != s_vv_1ring[i].end(); ++it) {
            fprintf(fp, "%d ", *it);
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);
#endif
    
    return 1;
}

int ARBFInterpolator::findNeighFaces1Ring()
{
    s_ff_1ring.resize(s_mesh->numFaces());
    
    for (int f = 0; f < s_mesh->numFaces(); ++f) {
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
    
#ifdef __TEST__
    string filename = string("DEBUG.face_faces_1ring.txt");
    FILE *fp = nullptr;
#ifdef _WIN32
    if (fopen_s(&fp, filename.c_str(), "w")) {
#else
    if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
#endif
        fprintf(stderr, "ERROR: open DEBUG.face_faces_1ring.txt failed ...\n");

        return 0;
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
#endif
    
    return 1;
}

int ARBFInterpolator::transformEigenvalues()
{
    s_eigenvalues.resize(s_mesh->numFaces());
    
#ifdef __TEST__
    string filename = "DEBUG.transformed_eigenvalue.txt";
    FILE *fp = nullptr;
#ifdef _WIN32
    if (fopen_s(&fp, filename.c_str(), "w")) {
#else
    if ((fp = fopen(filename.c_str(), "w")) == nullptr) {
#endif
        fprintf(stderr, "ERROR: open DEBUG.transformed_eigenvalue.txt failed ...\n");
    
        return 0;
    }
#endif
    
    for (int f = 0; f < s_mesh->numFaces(); ++f) {
        double eig2 = 0.0;
        double c = s_computeC(s_centers[f].lambda1, s_centers[f].lambda2);
        if (c != 0.0) {
            eig2 = 1.0 - exp(-3.8 * SQ(Config::C0/c));
        } else {
            eig2 = 1.0;
        }
        
#ifdef ISO_RBF
        eig2 = 1.0;
#endif
        
        s_eigenvalues[f] = make_pair(1.0, eig2);
        
#ifdef __TEST__
        fprintf(fp, "f%d: %lf\t%lf\n", f, s_eigenvalues[f].first,
                s_eigenvalues[f].second);
#endif
    }
    
#ifdef __TEST__
    fclose(fp);
#endif
        
    return 1;
}

void ARBFInterpolator::clean()
{
    if (s_centers != nullptr)
        delete [] s_centers;
    s_centers = nullptr;
}
    
pair<DBL2VECT, DBL2VECT>
ARBFInterpolator::s_computeEigen(const Matrix2d &m, double *lambda1, double *lambda2)
{
    assert(abs(m(0, 1) - m(1, 0)) < Config::EPSILON); // m should be symmetric
    
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
    if (abs(c) > Config::EPSILON) { // c != 0
        //        double u1 = c;
        //        double w1 = *lambda1 - a;
        
        double u1 = *lambda1 - b;
        double w1 = c;
        
        double p = sqrt(SQ(u1) + SQ(w1));
        v1.x = u1 / p;
        v1.y = w1 / p;
        
        v2.x = v1.y;
        v2.y = -v1.x;
    } else {
        v1.x = 1.0;
        v1.y = 0.0;
        v2.x = 0.0;
        v2.y = 1.0;
    }
    
    return make_pair(v1, v2);
}

inline double ARBFInterpolator::s_computeC(double mu1, double mu2)
{
    double sumMu = mu1 + mu2;
    
    if (sumMu > Config::EPSILON || sumMu < -Config::EPSILON) {
        return abs(mu1 - mu2) / sumMu;
    } else {
        return 1.0;
    }
}

inline double ARBFInterpolator::s_computeWeight(double r)
{
    return exp(-SQ(r) / SQ(30.0));
}

double ARBFInterpolator::s_computeDistance(const double *v0, const double *v1,
                            const Eigen::Matrix2d &T)
{
    Eigen::RowVector2d dis(v0[0] - v1[0], v0[1] - v1[1]);
    Eigen::RowVector2d tmp = dis * T;
    
    return std::sqrt(tmp(0)*dis(0) + tmp(1)*dis(1));
}

inline double ARBFInterpolator::s_phi(double r)
{
	switch (s_basis)
	{
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

// basis function MQ
inline double ARBFInterpolator::s_basisMQ(double r)
{
    return std::sqrt(SQ(r) + SQ(0.5));
}

// basis function IMQ
inline double ARBFInterpolator::s_basisIMQ(double r)
{
    return 1.0 / (sqrt(SQ(r) + SQ(1.8)));
}

// basis function Gaussian
inline double ARBFInterpolator::s_basisGaussian(double r)
{
    double c = 0.01;
    return exp(-SQ(c*r));
}

// basis function TPS
inline double ARBFInterpolator::s_basisTPS(double r)
{
    if (abs(r) < Config::EPSILON)
        return 0;
    else
        return SQ(r) * log(r);
}

bool ARBFInterpolator::s_isInTriangle(const double* x, int triangleID)
{
    const DBL7VECT &a = s_mesh->getVertices()[s_mesh->getFaces()[triangleID].a];
    const DBL7VECT &b = s_mesh->getVertices()[s_mesh->getFaces()[triangleID].b];
    const DBL7VECT &c = s_mesh->getVertices()[s_mesh->getFaces()[triangleID].c];
    
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
    
inline double ARBFInterpolator::s_computeRMSE()
{
    double sum = 0.0;
    //assert(s_imgOrig->X() == s_imgSR->X() &&
    //       s_imgOrig->Y() == s_imgSR->Y());
    
    for (int j = 0; j < s_imgOrig->Y(); ++j) {
        for (int i = 0; i < s_imgOrig->X(); ++i)
        {
            sum += SQ((double) s_imgOrig->getImageData()[j*s_imgOrig->X()+i] - s_imgSR->getImageData()[j*s_imgSR->X()+i]);
            //            sum += SQ((double) orig_img[j*N1+i] - (double) piecewise[j*N1+i]);
        }
    }
    
    return sqrt(sum / (s_imgOrig->X() * s_imgOrig->Y()));
}
