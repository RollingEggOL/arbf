//
//  RawivImage.cpp
//  arbf-test
//
//  Created by Ke Liu on 3/2/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#include <cstdio>
#include <cstring>
#include <cassert>
#include "../include/RawivImage.h"

using namespace std;

RawivImage::RawivImage() : m_data(nullptr), m_path(string()) {
    m_dim[0] = 0;
    m_dim[1] = 0;
    m_dim[2] = 0;
}

RawivImage::RawivImage(const unsigned dim[3], const float min[3], const float max[3]) {
    assert((dim[0] >= 0) && (dim[1] >= 0) && (dim[2] >= 0));
    m_data = new float[dim[0] * dim[1] * dim[2]];
    memset(m_data, 0, sizeof(float) * dim[0] * dim[1] * dim[2]);
    memcpy(m_dim, dim, sizeof(int) * 3);
    memcpy(m_min, min, sizeof(float) * 3);
    memcpy(m_max, max, sizeof(float) * 3);
//    m_dim[0] = dim[0];
//    m_dim[1] = dim[1];
//    m_dim[2] = dim[2];
//    
//    m_min[0] = min[0];
//    m_min[1] = min[1];
//    m_min[2] = min[2];
//    
//    m_max[0] = max[0];
//    m_max[1] = max[1];
//    m_max[2] = max[2];
    m_path = string();
}

RawivImage::RawivImage(const RawivImage &other) {
    m_data = new float[other.m_dim[0] * other.m_dim[1] * other.m_dim[2]];
	memcpy(m_data, other.m_data, sizeof(float) * other.m_dim[0] * other.m_dim[1] * other.m_dim[2]);
    memcpy(m_dim, other.m_dim, sizeof(int) * 3);
    memcpy(m_min, other.m_min, sizeof(float) * 3);
    memcpy(m_max, other.m_max, sizeof(float) * 3);
//    m_dim[0] = other.m_dim[0];
//    m_dim[1] = other.m_dim[1];
//    m_dim[2] = other.m_dim[2];
    
//    m_min[0] = other.m_min[0];
//    m_min[1] = other.m_min[1];
//    m_min[2] = other.m_min[2];
    
//    m_max[0] = other.m_max[0];
//    m_max[1] = other.m_max[1];
//    m_max[2] = other.m_max[2];
    m_path = other.m_path;
}

RawivImage::~RawivImage() {
    delete [] m_data;
    m_data = nullptr;
}

RawivImage& RawivImage::operator=(const RawivImage &other) {
    if (this == &other)
        return *this;
    
    this->~RawivImage();
    m_data = new float[other.m_dim[0] * other.m_dim[1] * other.m_dim[2]];
    memcpy(m_data, other.m_data, sizeof(float) * other.m_dim[0] * other.m_dim[1] * other.m_dim[2]);
    memcpy(m_dim, other.m_dim, sizeof(int) * 3);
    memcpy(m_min, other.m_min, sizeof(float) * 3);
    memcpy(m_max, other.m_max, sizeof(float) * 3);
//    m_dim[0] = other.m_dim[0];
//    m_dim[1] = other.m_dim[1];
//    m_dim[2] = other.m_dim[2];
    
//    m_min[0] = other.m_min[0];
//    m_min[1] = other.m_min[1];
//    m_min[2] = other.m_min[2];
    
//    m_max[0] = other.m_max[0];
//    m_max[1] = other.m_max[1];
//    m_max[2] = other.m_max[2];
    m_path = other.m_path;
    return *this;
}

void RawivImage::write() {
    FILE* fp = nullptr;
    
    if ((fp = fopen(m_path.c_str(), "wb")) == nullptr) {
        fprintf(stderr, "ERROR: write image %s \n", m_path.c_str());
        exit(EXIT_FAILURE);
    }
    
    float min[] = { m_min[0], m_min[1], m_min[2] };
    float max[] = { m_max[0], m_max[1], m_max[2] };
    int dim[] = { m_dim[0], m_dim[1], m_dim[2] };
    float origin[] = { m_min[0], m_min[1], m_min[2] };
    float span[] = {
        span[0] = (m_max[0] - m_min[0]) / (float) (m_dim[0]-1),
        span[1] = (m_max[1] - m_min[1]) / (float) (m_dim[1]-1),
        span[2] = (m_max[2] - m_min[2]) / (float) (m_dim[2]-1)
    };
    int numVertices = m_dim[0] * m_dim[1] * m_dim[2];
    int numCells = (m_dim[0]-1) * (m_dim[1]-1) * (m_dim[2]-1);
    float *data = new float[m_dim[0] * m_dim[1] * m_dim[2]];
    memcpy(data, m_data, sizeof(float) * m_dim[0] * m_dim[1] * m_dim[2]);
    
    m_swap_buffer((char *) min, 3, sizeof(float));
    m_swap_buffer((char *) max, 3, sizeof(float));
    m_swap_buffer((char *) &numVertices, 1, sizeof(unsigned int));
    m_swap_buffer((char *) &numCells, 1, sizeof(unsigned int));
    m_swap_buffer((char *) dim, 3, sizeof(int));
    m_swap_buffer((char *) origin, 3, sizeof(float));
    m_swap_buffer((char *) span, 3, sizeof(float));
    m_swap_buffer((char *) data, m_dim[0] * m_dim[1] * m_dim[2], sizeof(float));
    
    fwrite(min, sizeof(float), 3, fp);
    fwrite(max, sizeof(float), 3, fp);
    fwrite(&numVertices, sizeof(unsigned int), 1, fp);
    fwrite(&numCells, sizeof(unsigned int), 1, fp);
    fwrite(dim, sizeof(int), 3, fp);
    fwrite(origin, sizeof(float), 3, fp);
    fwrite(span, sizeof(float), 3, fp);

    for (int k = 0; k < m_dim[2]; k++) {
        for (int j = 0; j < m_dim[1]; j++) {
            for (int i = 0; i < m_dim[0]; i++) {
                fwrite(&data[k*m_dim[1]*m_dim[0] + j*m_dim[0] + i], sizeof(float), 1, fp);
            }
        }
    }

    delete [] data;
    fclose(fp);
}

void RawivImage::m_swap_buffer(char *buffer, int count, int typesize) {
    char sbuf[4];
    int i;
    int temp = 1;
    unsigned char* chartempf = (unsigned char*) &temp;
    
    if (chartempf[0] > '\0') {
        // swapping isn't necessary on single byte data
        if (typesize == 1)
            return;
        
        
        for (i = 0; i < count; i++){
            memcpy(sbuf, buffer+(i*typesize), typesize);
            
            switch (typesize) {
                case 2: {
                    buffer[i*typesize] = sbuf[1];
                    buffer[i*typesize+1] = sbuf[0];
                    break;
                }
                case 4: {
                    buffer[i*typesize] = sbuf[3];
                    buffer[i*typesize+1] = sbuf[2];
                    buffer[i*typesize+2] = sbuf[1];
                    buffer[i*typesize+3] = sbuf[0];
                    break;
                }
                default:
                    break;
            }
        }
    }
}
