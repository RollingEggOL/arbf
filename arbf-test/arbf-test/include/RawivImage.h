//
//  RawivImage.h
//  arbf-test
//
//  Created by Ke Liu on 2/29/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#ifndef __arbf_test_RawivImage_h__
#define __arbf_test_RawivImage_h__

#include <string>

class RawivImage
{
public:
    /*
     * Default constructor
     */
    RawivImage();
    
    /*
     * Create image and allocate memory of \param x * \param y.
     * \param dim: dimension in X, Y, Z
     * \param min: minimum coordinates
     * \param max: maximum coordinates
     */
    RawivImage(const unsigned dim[3], const float min[3], const float max[3]);
    
    /*
     * Copy constructor
     */
    RawivImage(const RawivImage &other);
    
    /*
     * Destructor
     */
    virtual ~RawivImage();
    
    /*
     * Assignment operator
     */
    RawivImage& operator=(const RawivImage &other);
    
    const int* getDimension() const {
        return m_dim;
    }
    
    void setDimension(const int dim[3]) {
        m_dim[0] = dim[0];
        m_dim[1] = dim[1];
        m_dim[2] = dim[2];
    }
    
    const char* getPath() const {
        return m_path.c_str();
    }
    
    void setPath(const char* path) {
        m_path = std::string(path);
    }
    
    float* getImageData() const {
        return m_data;
    }
    
    void setImageData(float *data) {
        m_data = data;
    }
    
    void write();

private:
    void m_swap_buffer(char *buffer, int count, int typesize);
    
private:
	float *m_data;
    int m_dim[3]; // dimension
    float m_min[3]; // minimum coordinates
    float m_max[3]; // maximum coordinates
    std::string m_path; // image path
};

#endif /* defined(__arbf_test_RawivImage_h__) */
