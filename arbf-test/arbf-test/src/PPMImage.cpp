//
//  PPMImage.cpp
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#include <cstdio>
#include <cstring>
#include <cassert>
#include "../include/PPMImage.h"

using namespace std;

PPMImage::PPMImage() : m_data(nullptr), m_x(0), m_y(0), m_maxRaw(0.0), m_path(string()) {}

PPMImage::PPMImage(int x, int y) {
    assert((x >= 0) && (y >= 0));
    m_data = new double[x * y];
    memset(m_data, 0, sizeof(double) * x * y);
    m_x = x;
    m_y = y;
    m_maxRaw = 0;
    m_path = string();
}

PPMImage::PPMImage(const char *path) : m_data(nullptr), m_x(0), m_y(0), m_maxRaw(0.0) {
    assert(path != nullptr);
    m_path = string(path);
    read();
}

PPMImage::PPMImage(const PPMImage &other) {
    m_data = new double[other.m_x * other.m_y];
	memcpy(m_data, other.m_data, sizeof(double) * other.m_x * other.m_y);
    m_x = other.m_x;
    m_y = other.m_y;
    m_maxRaw = other.m_maxRaw;
    m_path = other.m_path;
}

PPMImage::~PPMImage() {
    delete [] m_data;
    m_data = nullptr;
}

PPMImage& PPMImage::operator=(const PPMImage &other) {
    if (this == &other)
        return *this;
    
    this->~PPMImage();
	m_data = new double[other.m_x * other.m_y];
	memcpy(m_data, other.m_data, sizeof(double) * other.m_x * other.m_y);
    m_x = other.m_x;
    m_y = other.m_y;
    m_maxRaw = other.m_maxRaw;
    m_path = other.m_path;
    return *this;
}

void PPMImage::read() {
    int c;
    int i, j;
    unsigned char *line;
    const int BUF_SIZE = 1024;
    char buf[BUF_SIZE];
    FILE* fp = nullptr;
    
    if ((fp = fopen(m_path.c_str(), "rb")) == nullptr) {
        fprintf(stderr, "ERROR: read image %s \n", m_path.c_str());
        exit(EXIT_FAILURE);
    }
    
    while ((c = fgetc(fp)) == '#') {
        fgets(buf, BUF_SIZE, fp);
    }
    
    ungetc(c, fp);
    
    if (fscanf(fp, "P%d\n", &c) != 1) {
        fprintf(stderr, "ERROR: read image ... \n");
        exit(EXIT_FAILURE);
    }
    
    if (c != 6 && c != 2) {
        fprintf(stderr, "ERROR: read image ...\n");
        exit(EXIT_FAILURE);
    }
    
    // PPM format
    if (c == 6) {
        while ((c = fgetc(fp)) == '#') {
            fgets(buf, BUF_SIZE, fp);
        }
        
        ungetc(c, fp);
        
        if (fscanf(fp, "%d %d\n%lf\n", &m_x, &m_y, &m_maxRaw) != 3) {
            fprintf(stderr, "ERROR: failed to read width/height/max\n");
            exit(EXIT_FAILURE);
        }
        
		m_data = new double[m_x * m_y];
        getc(fp);
        line = new unsigned char[m_x * 3];
        for (j = 0; j < m_y; j++) {
            fread(line, 1, m_x*3, fp);
            for (i = 0; i < m_x*3; i += 3) {
                m_data[j*m_x+i/3] = line[i];
            }
        }
        delete [] line;
    }
        
    fclose(fp);
}

void PPMImage::write() const {
    double isoValues[] = {64.0, 32.0, 16.0};
    
    FILE* fp = nullptr;
    
    if ((fp = fopen(m_path.c_str(), "wb")) == nullptr) {
        fprintf(stderr, "ERROR: write image %s \n", m_path.c_str());
        exit(EXIT_FAILURE);
    }
        
    fprintf(fp, "P6\n%d %d\n%lf\n", m_x, m_y, m_maxRaw);

    for (int j = 0; j < m_y; j++) {
        for (int i = 0; i < m_x; i++) {
            
            if (m_data[j*m_x+i] <= isoValues[0] && m_data[j*m_x+i] > isoValues[1]) {
                fputc(255, fp);
                fputc(0, fp);
                fputc(0, fp);
            } else if (m_data[j*m_x+i] <= isoValues[1] && m_data[j*m_x+i] > isoValues[2]) {
                fputc(0, fp);
                fputc(255, fp);
                fputc(0, fp);
            } else if (m_data[j*m_x+i] <= isoValues[2]) {
                fputc(0, fp);
                fputc(0, fp);
                fputc(255, fp);
            } else {
                fputc((int) m_data[j*m_x+i], fp);
                fputc((int) m_data[j*m_x+i], fp);
                fputc((int) m_data[j*m_x+i], fp);
            }
        }
    }

    fclose(fp);
}