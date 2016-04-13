//
//  main.cpp
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <string>
#include <tuple>
//#include <iostream>
//#include <boost/format.hpp>
#include "../include/Config.h"
#include "../include/Common.h"
#include "../include/PPMImage.h"
#include "../include/RawivImage.h"
#include "../include/TriMesh.h"
#include "../include/ARBFInterpolator.h"

using namespace std;
//using namespace boost;
using namespace Eigen;

//void initLogger(const string& loggerName)
//{
//    Config::loggerName = "logs/";
//    
//    string info_log = Config::loggerName + Config::loggerINFO;
//    string warning_log = Config::loggerName + Config::loggerWARNING;
//    string error_log = Config::loggerName + Config::loggerERROR;
//    string fatal_log = Config::loggerName + Config::loggerFATAL;
//    
//    google::InitGoogleLogging(loggerName.c_str()); // glog init
//    google::SetLogDestination(google::INFO, info_log.c_str());
//    google::SetLogDestination(google::WARNING, warning_log.c_str());
//    google::SetLogDestination(google::ERROR, error_log.c_str());
//    google::SetLogDestination(google::FATAL, fatal_log.c_str());
//}

int main(int argc, const char * argv[])
{
	if (argc != 5) {
		fprintf(stderr, "ERROR: incorrect command format. Usage: arbf-test <proj_dir> <mesh_filename> <neighborhood_size> <model_dimension>\n");
		return EXIT_FAILURE;
	}
    
    if (argv[1] == nullptr || strcmp(argv[1], "") == 0) {
        fprintf(stderr, "ERROR: incorrect project directory.\n");
        return EXIT_FAILURE;
    }
    
    if (argv[2] == nullptr || strcmp(argv[2], "") == 0) {
        fprintf(stderr, "ERROR: incorrect mesh filename.\n");
        return EXIT_FAILURE;
    }
    
    if (atoi(argv[3]) != 0 && atoi(argv[3]) != 1 && atoi(argv[3]) != 2) {
        fprintf(stderr, "ERROR: incorrect neighborhood size.\n\t0: no neighborhood, 1: 1-ring, 2: 2-ring.\n");
        return EXIT_FAILURE;
    }
    
    if (atoi(argv[4]) != 2 && atoi(argv[4]) != 3) {
        fprintf(stderr, "ERROR: incorrect model dimension.\n\t2: use 2D model, 3: use 3D model.\n");
        return EXIT_FAILURE;
    }
    
    const string PROJ_DIR = argv[1];
    
    if (atoi(argv[3]) == 2) {
        Config::isUsing2Ring = true;
    } else {
        Config::isUsing2Ring = false;
    }
    
    if (atoi(argv[4]) == 3) {
        Config::dim = 3;
    }
    
    clock_t start, span, total;
    
    // read mesh
    start = total = clock();
    string file_path = PROJ_DIR + argv[2];
    TriMesh* mesh = new TriMesh(file_path.c_str());
    printf("nv = %d, nf = %d\n", mesh->getNumVertices(), mesh->getNumFaces());
    span = clock() - start;
    total += span;
    printf("Reading mesh ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
    
    ARBFInterpolator::setMesh(mesh);
	ARBFInterpolator::setBasisType(BasisType::MQ);
    
    // initialize triangle centers for each face
    start = clock();
    ARBFInterpolator::calculateCenters();
    span = clock() - start;
    total += span;
    printf("Computing centers ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
    
    // initialize triangle edge centers, for 2D model only
    if (Config::dim == 2) {
        start = clock();
        ARBFInterpolator::calculateEdgeCenters();
        span = clock() - start;
        total += span;
        printf("Computing edge centers ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
    }
    
	// find neighboring vertices for each vertex and
    // neighboring faces for each face (center)
    start = clock();
    ARBFInterpolator::findNeighVertices1Ring();
    ARBFInterpolator::findNeighFaces1Ring();
    span = clock() - start;
    total += span;
    printf("Finding 1-Ring neighbors ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
    
    if (Config::isUsing2Ring) {
        start = clock();
        ARBFInterpolator::findNeighVertices2Ring();
        ARBFInterpolator::findNeighFaces2Ring();
        span = clock() - start;
        total += span;
        printf("Finding 2-Ring neighbors ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
    }
    
    // compute metrics
    start = clock();
    ARBFInterpolator::computeMetric();
    span = clock() - start;
    total += span;
    printf("Computing metrics ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
    
    // interpolate
    start = clock();
    tuple<double*, double, int, int, int> result = ARBFInterpolator::interpolate();
    
    if (Config::dim == 2) {
        PPMImage* output = new PPMImage(get<2>(result), get<3>(result));
        output->setMaxIntensity(get<1>(result));
        string output_path = PROJ_DIR + "output.ppm";
        output->setPath(output_path.c_str());
        output->setImageData(get<0>(result));
        output->write();
    } else {
        int dim[] = { get<2>(result), get<3>(result), get<4>(result) };
        float min[] = { (float) mesh->getMinX(), (float) mesh->getMinY(), (float) mesh->getMinZ() };
        float max[] = { (float) mesh->getMaxX(), (float) mesh->getMaxY(), (float) mesh->getMaxZ() };
        RawivImage* output = new RawivImage(dim, min, max);
        string output_path = PROJ_DIR + "output.rawiv";
        output->setPath(output_path.c_str());
        
        float *data = new float[dim[0]*dim[1]*dim[2]];
        
        for (int k = 0; k < dim[2]; k++) {
            for (int j = 0; j < dim[1]; j++) {
                for (int i = 0; i < dim[0]; i++) {
                    data[k*dim[1]*dim[0] + j*dim[0] + i] = (float) get<0>(result)[k*dim[1]*dim[0] + j*dim[0] + i];
                }
            }
        }
        
        output->setImageData(data);
        printf("dim = [%d, %d, %d]\n", dim[0], dim[1], dim[2]);
        printf("min = [%f, %f, %f]\n", min[0], min[1], min[2]);
        printf("max = [%f, %f, %f]\n", max[0], max[1], max[2]);
        
        output->write();
        
        delete [] data;
    }
    
    span = clock() - start;
    total += span;
    printf("Interpolating ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
    
    // clean up
    start = clock();
    ARBFInterpolator::clean();
    span = clock() - start;
    total += span;
    printf("Cleaning ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
    
    printf("Running time = %lf seconds.\n", (double) total / CLOCKS_PER_SEC);
	
	return EXIT_SUCCESS;
}
