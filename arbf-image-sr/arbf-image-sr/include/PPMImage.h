//
//  PPMImage.h
//  arbf-image-super_resolution
//
//  Created by Ke Liu on 2/26/15.
//  Copyright (c) 2015 Ke Liu. All rights reserved.
//

#ifndef __arbf_image_super_resolution__PPMImage__
#define __arbf_image_super_resolution__PPMImage__

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
     * \param x: # of columns of image
     * \param y: # of rows of image
     * \param scaleX: scaling factor for columns
     * \param scaleY: scaling factor for rows
     */
    PPMImage(int x, int y, double scaleX = 1.0, double scaleY = 1.0);
    
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
    
    int X() const;
    void setX(int x);
    
    int Y() const;
    void setY(int y);
    
    int getMaxIntensity() const;
    void setMaxIntensity(int max);
    
	double getScaleX() const;
	void setScaleX(double scaleX);

	double getScaleY() const;
	void setScaleY(double scaleY);

    const char* getPath() const;
    void setPath(const char* path);
    
    double* getImageData();
    const double* getImageData() const;
    
    int read();
    int write() const;
    
private:
	double *m_data;
    int m_x; // # of cols
    int m_y; // # of rows
    int m_maxRaw; // maximum intensity
	double m_scaleX; // scaling factor of cols
	double m_scaleY; // scaling factor of rows
    std::string m_path; // image path
    
private:
    static const int BUF_SIZE; // size of buffer when reading image
};

#endif /* defined(__arbf_image_super_resolution__PPMImage__) */
