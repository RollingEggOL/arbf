//
//  PPMImage.h
//  arbf-test
//
//  Created by Ke Liu on 2/19/16.
//  Copyright (c) 2016 Ke Liu. All rights reserved.
//

#ifndef __arbf_test_PPMImage_h__
#define __arbf_test_PPMImage_h__

#include <string>

class PPMImage
{
public:
    /*
     * Default constructor
     */
    PPMImage();
    
    /*
     * Create image and allocate memory of \param x * \param y.
     * \param x: # of columns
     * \param y: # of rows
     */
    PPMImage(int x, int y);
    
    /*
     * Create image by reading image file.
     * \param path: image file path
     */
    PPMImage(const char *path);
    
    /*
     * Copy constructor
     */
    PPMImage(const PPMImage &other);
    
    /*
     * Destructor
     */
    virtual ~PPMImage();
    
    /*
     * Assignment operator
     */
    PPMImage& operator=(const PPMImage &other);
    
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
    
    int getMaxIntensity() const {
        return m_maxRaw;
    }
    
    void setMaxIntensity(double max) {
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
    void read();
    
private:
	double *m_data;
    int m_x; // # of cols
    int m_y; // # of rows
    double m_maxRaw; // maximum intensity
    std::string m_path; // image path
};

#endif /* defined(__arbf_test_PPMImage_h__) */
