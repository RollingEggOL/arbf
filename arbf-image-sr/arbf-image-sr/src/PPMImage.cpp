//
//  PPMImage.cpp
//  arbf-image-super_resolution
//
//  Created by Ke Liu on 2/26/15.
//  Copyright (c) 2015 Ke Liu. All rights reserved.
//

#include <cstdio>
#include <cstring>
#include <cassert>
//#include <glog/logging.h>
//#include <boost/format.hpp>
#include "../include/PPMImage.h"

using namespace std;
//using namespace boost;

const int PPMImage::BUF_SIZE = 1024;

PPMImage::PPMImage() : m_data(nullptr), m_x(0), m_y(0), m_maxRaw(0), m_scaleX(1.0), m_scaleY(1.0), m_path(string()) {}

PPMImage::PPMImage(int x, int y, double scaleX, double scaleY)
{
    assert((x >= 0) && (y >= 0));
    m_data = new double[x * y];
    memset(m_data, 0, sizeof(double) * x * y);
    m_x = x;
    m_y = y;
    m_maxRaw = 0;
	m_scaleX = scaleX;
	m_scaleY = scaleY;
    m_path = string();
}

PPMImage::PPMImage(const char *path) : m_data(nullptr), m_x(0), m_y(0), m_maxRaw(0), m_scaleX(1.0), m_scaleY(1.0)
{
    assert(path != nullptr);
    m_path = string(path);
    read();
}

PPMImage::PPMImage(const PPMImage &other)
{
    m_data = new double[other.m_x * other.m_y];
	memcpy(m_data, other.m_data, sizeof(double) * other.m_x * other.m_y);
    m_x = other.m_x;
    m_y = other.m_y;
    m_maxRaw = other.m_maxRaw;
	m_scaleX = other.m_scaleX;
	m_scaleY = other.m_scaleY;
    m_path = other.m_path;
}

PPMImage::~PPMImage()
{
    delete [] m_data;
    m_data = nullptr;
}

PPMImage& PPMImage::operator=(const PPMImage &other)
{
    if (this == &other)
        return *this;
    
    if (m_data != nullptr) {
        delete [] m_data;
    }
	m_data = new double[other.m_x * other.m_y];
	memcpy(m_data, other.m_data, sizeof(double) * other.m_x * other.m_y);
    m_x = other.m_x;
    m_y = other.m_y;
	m_scaleX = other.m_scaleX;
	m_scaleY = other.m_scaleY;
    m_maxRaw = other.m_maxRaw;
    
    return *this;
}

int PPMImage::X() const
{
    return m_x;
}

void PPMImage::setX(int x)
{
    m_x = x;
}

int PPMImage::Y() const
{
    return m_y;
}

void PPMImage::setY(int y)
{
    m_y = y;
}

int PPMImage::getMaxIntensity() const
{
    return m_maxRaw;
}

void PPMImage::setMaxIntensity(int max)
{
    m_maxRaw = max;
}

double PPMImage::getScaleX() const
{
	return m_scaleX;
}

void PPMImage::setScaleX(double scaleX)
{
	m_scaleX = scaleX;
}

double PPMImage::getScaleY() const
{
	return m_scaleY;
}

void PPMImage::setScaleY(double scaleY)
{
	m_scaleY = scaleY;
}

const char* PPMImage::getPath() const
{
    return m_path.c_str();
}

void PPMImage::setPath(const char *path)
{
    m_path = string(path);
}

double* PPMImage::getImageData()
{
    return m_data;
}

const double* PPMImage::getImageData() const
{
    return m_data;
}

int PPMImage::read()
{
    int c;
    int i, j;
    unsigned char *line;
    char buf[BUF_SIZE];
    
    FILE* fp = nullptr;
    
#ifdef _WIN32
    if (fopen_s(&fp, m_path.c_str(), "rb")) {
#else
    if ((fp = fopen(m_path.c_str(), "rb")) == nullptr) {
#endif
        fprintf(stderr, "ERROR: read image %s \n", m_path.c_str());
        
        return 0;
    }
    
    while ((c = fgetc(fp)) == '#') {
        fgets(buf, BUF_SIZE, fp);
    }
    ungetc(c, fp);
#ifdef _WIN32
    if (fscanf_s(fp, "P%d\n", &c) != 1) {
#else
    if (fscanf(fp, "P%d\n", &c) != 1) {
#endif
        fprintf(stderr, "ERROR: read image ... \n");
        
        return 0;
    }
    if (c != 6 && c != 2) {
        fprintf(stderr, "ERROR: read image ...\n");
        
        return 0;
    }
    
    if (c == 6) {
        while ((c = fgetc(fp)) == '#') {
            fgets(buf, BUF_SIZE, fp);
        }
        ungetc(c, fp);
#ifdef _WIN32
        if (fscanf_s(fp, "%d%d%d", &m_x, &m_y, &m_maxRaw) != 3) {
#else
        if (fscanf(fp, "%d%d%d", &m_x, &m_y, &m_maxRaw) != 3) {
#endif
            fprintf(stderr, "ERROR: failed to read width/height/max\n");
            
            return 0;
        }
        
		m_data = new double[m_x * m_y];
        getc(fp);
        line = new unsigned char[m_x * 3];
        for (j = 0; j < m_y; j++) {
            fread(line, 1, m_x*3, fp);
            for (i = 0; i < m_x*3; i += 3) {
                m_data[j * m_x + i / 3] = line[i];
            }
        }
        delete [] line;
    }
        
    fclose(fp);
    
    return 1;
}

int PPMImage::write() const
{
    FILE* fp = nullptr;
    
#ifdef _WIN32
    if (fopen_s(&fp, m_path.c_str(), "wb")) {
#else
    if ((fp = fopen(m_path.c_str(), "wb")) == nullptr) {
#endif
        fprintf(stderr, "ERROR: write image %s \n", m_path.c_str());
        return 0;
    }
        
    fprintf(fp, "P6\n%d %d\n%d\n", m_x, m_y, m_maxRaw);

    for (int j = m_y-1; j >= 0; j--) {
        for (int i = 0; i < m_x; i++) {
            fputc((int) m_data[j*m_x+i], fp);
            fputc((int) m_data[j*m_x+i], fp);
            fputc((int) m_data[j*m_x+i], fp);
        }
    }

    fclose(fp);
        
    return 1;
}