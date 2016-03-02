//
//  PLYImage.h
//  arbf-test
//
//  Created by Ke Liu on 2/29/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#ifndef __arbf_test_PLYImage_h__
#define __arbf_test_PLYImage_h__

#include <string>

class PLYImage
{
public:
    /*
     * Default constructor
     */
    PLYImage();
    
    /*
     * Create image and allocate memory of \param x * \param y.
     * \param x: # of columns
     * \param y: # of rows
     * \param z: # of levels
     */
    PLYImage(int x, int y, int z);
    
    /*
     * Copy constructor
     */
    PLYImage(const PLYImage &other);
    
    /*
     * Destructor
     */
    virtual ~PLYImage();
    
    /*
     * Assignment operator
     */
    PLYImage& operator=(const PLYImage &other);
    
    int getX() const {
        return m_x;
    }
    
    void setX(int x) {
        m_x = x;
    }
    
    int getY() const {
        return m_y;
    }
    
    void setY(int y) {
        m_y = y;
    }
    
    int getZ() const {
        return m_z;
    }
    
    void setZ(int z) {
        m_z = z;
    }
    
    int getMaxIntensity() const {
        return m_maxRaw;
    }
    
    void setMaxIntensity(int max) {
        m_maxRaw = max;
    }
    
    const char* getPath() const {
        return m_path.c_str();
    }
    
    void setPath(const char* path) {
        m_path = std::string(path);
    }
    
    double* getImageData() const {
        return m_data;
    }
    
    void setImageData(double *data) {
        m_data = data;
    }
    
    void write() const;
    
private:
	double *m_data;
    int m_x; // # of cols
    int m_y; // # of rows
    int m_z; // # of levels
    int m_maxRaw; // maximum intensity
    std::string m_path; // image path
};

#endif /* defined(__arbf_test_PLYImage_h__) */
