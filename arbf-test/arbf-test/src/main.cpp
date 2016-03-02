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
#include "../include/PLYImage.h"
#include "../include/TriMesh.h"
#include "../include/ARBFInterpolator.h"

//#define __TEST__

using namespace std;
//using namespace boost;
using namespace Eigen;

//void initLogger(const string& loggerName)
//{
//#ifdef _WIN32
//    Config::loggerName = "logs\\";
//#else
//    Config::loggerName = "logs/";
//#endif
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
        PLYImage* output = new PLYImage(get<2>(result), get<3>(result), get<4>(result));
        output->setMaxIntensity(get<1>(result));
        string output_path = PROJ_DIR + "output3D.ply";
        output->setPath(output_path.c_str());
        output->setImageData(get<0>(result));
        output->write();
        
        const int NUM_PPMs = 10;
        int step = get<4>(result) / NUM_PPMs;
        double *outData = new double[get<2>(result) * get<3>(result)];
        
        for (int k = 0; k < get<4>(result); k += step) {
            for (int j = 0; j < get<3>(result); j++) {
                for (int i = 0; i < get<2>(result); i++) {
                    outData[j*get<2>(result)+i] = get<0>(result)[k*get<3>(result)*get<2>(result) + j*get<2>(result) + i];
                }
            }
            
            PPMImage *out = new PPMImage(get<2>(result), get<3>(result));
            out->setMaxIntensity(get<1>(result));
            string outPath = PROJ_DIR + "output3D_z" + to_string(k) + ".ppm";
            out->setPath(outPath.c_str());
            out->setImageData(outData);
            out->write();
        }
        
        delete [] outData;
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
