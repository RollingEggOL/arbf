//
//  PLYImage.cpp
//  arbf-test
//
//  Created by Ke Liu on 2/29/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#include <cstdio>
#include <cstring>
#include <cassert>
#include "../include/PLYImage.h"

using namespace std;

PLYImage::PLYImage() : m_data(nullptr), m_x(0), m_y(0), m_z(0), m_maxRaw(0), m_path(string()) {}

PLYImage::PLYImage(int x, int y, int z) {
    assert((x >= 0) && (y >= 0) && (z >= 0));
    m_data = new double[x * y * z];
    memset(m_data, 0, sizeof(double) * x * y * z);
    m_x = x;
    m_y = y;
    m_z = z;
    m_maxRaw = 0;
    m_path = string();
}

PLYImage::PLYImage(const PLYImage &other) {
    m_data = new double[other.m_x * other.m_y * other.m_z];
	memcpy(m_data, other.m_data, sizeof(double) * other.m_x * other.m_y * other.m_z);
    m_x = other.m_x;
    m_y = other.m_y;
    m_z = other.m_z;
    m_maxRaw = other.m_maxRaw;
    m_path = other.m_path;
}

PLYImage::~PLYImage() {
    delete [] m_data;
    m_data = nullptr;
}

PLYImage& PLYImage::operator=(const PLYImage &other) {
    if (this == &other)
        return *this;
    
    this->~PLYImage();
	m_data = new double[other.m_x * other.m_y * other.m_z];
	memcpy(m_data, other.m_data, sizeof(double) * other.m_x * other.m_y * other.m_z);
    m_x = other.m_x;
    m_y = other.m_y;
    m_z = other.m_z;
    m_maxRaw = other.m_maxRaw;
    m_path = other.m_path;
    return *this;
}

void PLYImage::write() const {
    FILE* fp = nullptr;
    
    if ((fp = fopen(m_path.c_str(), "w")) == nullptr) {
        fprintf(stderr, "ERROR: write image %s \n", m_path.c_str());
        exit(EXIT_FAILURE);
    }
        
    fprintf(fp, "ply\nformat ascii 1.0\n");
    fprintf(fp, "element vertex %d\n", m_x * m_y * m_z);
    fprintf(fp, "property float x\nproperty float y\nproperty float z\n");
    fprintf(fp, "property uchar red\nproperty uchar green\nproperty uchar blue\n");
    fprintf(fp, "end_header\n");

    for (int k = 0; k < m_z; k++) {
        for (int j = 0; j < m_y; j++) {
            for (int i = 0; i < m_x; i++) {
                int intensity = (int) m_data[k*m_y*m_x + j*m_x + i];
                fprintf(fp, "%d %d %d %d %d %d\n", i, j, k, intensity, intensity, intensity);
            }
        }
    }

    fclose(fp);
}