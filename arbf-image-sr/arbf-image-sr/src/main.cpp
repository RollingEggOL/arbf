//
//  main.cpp
//  arbf-image
//
//  Created by Ke Liu on 1/23/14.
//  Copyright (c) 2014 Ke Liu. All rights reserved.
//

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <string>
//#include <iostream>
//#include <glog/logging.h>
//#include <boost/format.hpp>
#include "../include/Config.h"
#include "../include/Common.h"
#include "../include/PPMImage.h"
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
		fprintf(stderr, "ERROR: Incorrect command format. Correct is:\n"
			"\tarbf-image-super_resolution <filename_prefix> <width_scaling> <height_scaling> <neighborhood>\n");

		return EXIT_FAILURE;
	}
    
    if (argv[1] == nullptr ||
        strcmp(argv[1], "") == 0) {
        fprintf(stderr, "ERROR: Incorrect filename_prefix.\n");
        
        return EXIT_FAILURE;
    }

	if ((stod(argv[2]) < Config::EPSILON) || (stod(argv[3]) < Config::EPSILON)) {
		fprintf(stderr, "ERROR: Scaling factors must be positive.\n");

		return EXIT_FAILURE;
	}
    
    if (atoi(argv[4]) != 1 && atoi(argv[4]) != 2) {
        fprintf(stderr, "ERROR: Incorrect neighborhood number.\n"
                "\t1 for 1-ring, 2 for 2-ring.\n");
        
        return EXIT_FAILURE;
    }
    
//    initLogger(argv[0]);
    clock_t start, span, total;
    
    // read mesh
    start = total = clock();
    string file = argv[1];
    
#ifdef _WIN32
    string filename = "data\\" + file + ".tri";
#else
    string filename = "data/" + file + ".tri";
#endif

    TriMesh* mesh = new TriMesh(filename.c_str());
    printf("nv = %d, nf = %d\n", mesh->numVertices(), mesh->numFaces());
    
    span = clock() - start;
    total += span;
    printf("Reading mesh ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
    
    // read original image
    start = clock();
#ifdef _WIN32
    filename = "data\\" + file + ".ppm";
#else
    filename = "data/" + file + ".ppm";
#endif
    PPMImage* imgOriginal = new PPMImage(filename.c_str());
    printf("Original dimension: %d X %d\n", imgOriginal->X(), imgOriginal->Y());
    span = clock() - start;
    total += span;
    printf("Reading original image ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
	PPMImage* imgSR = new PPMImage(stod(argv[2]) * imgOriginal->X(),
                                   stod(argv[3]) * imgOriginal->Y(),
                                   stod(argv[2]), stod(argv[3]));
    imgSR->setPath((string("out.ppm")).c_str());
    imgSR->setMaxIntensity(imgOriginal->getMaxIntensity());
    printf("Recovering dimension: %d X %d\n", imgSR->X(), imgSR->Y());
    
    ARBFInterpolator::setMesh(mesh);
    ARBFInterpolator::setOriginalImage(imgOriginal);
    ARBFInterpolator::setSRImage(imgSR);
	ARBFInterpolator::setBasisType(BasisType::MQ);
    if (atoi(argv[4]) == 2) {
        Config::isUsing2Ring = true;
    } else {
        Config::isUsing2Ring = false;
    }
    
    // initialize triangle centers for each face
    start = clock();
    ARBFInterpolator::calculateCenters();
    span = clock() - start;
    total += span;
    printf("Computing centers ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
	
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

    // transform eigenvalues
    start = clock();
	ARBFInterpolator::transformEigenvalues();
    span = clock() - start;
    total += span;
    printf("Transforming eigenvalues ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
	
	// compute metrics
    start = clock();
    ARBFInterpolator::computeMetrics();
    span = clock() - start;
    total += span;
    printf("Computing tensor T ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
	
	// interpolate
    start = clock();
    if (Config::isUsingWeightedInterpolation) {
        ARBFInterpolator::weightedInterpolate();
        span = clock() - start;
        total += span;
        printf("Weighted interpolating ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
    } else {
        ARBFInterpolator::interpolate();
        span = clock() - start;
        total += span;
        printf("Interpolating ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
    }
	
    // compute PSNR
    start = clock();
    printf("PSNR = %lf\n", ARBFInterpolator::computePSNR());
    span = clock() - start;
    total += span;
    printf("Computing PSNR ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
    
	// outputm
    start = clock();
    imgSR->write();
    span = clock() - start;
    total += span;
    printf("Printing result ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
	
	// clean up
    start = clock();
    ARBFInterpolator::clean();
    span = clock() - start;
    total += span;
    printf("Cleaning ... %lf ms\n", (double) span * 1e3 / CLOCKS_PER_SEC);
    
    printf("Running time = %lf seconds.\n", (double) total / CLOCKS_PER_SEC);
	
	return EXIT_SUCCESS;
}
