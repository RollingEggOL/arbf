//
//  Common.h
//  arbf-image-super_resolution
//
//  Created by Ke Liu on 2/27/15.
//  Copyright (c) 2015 Ke Liu. All rights reserved.
//

#ifndef arbf_image_super_resolution_Common_h
#define arbf_image_super_resolution_Common_h

#define SQ(x) ((x) * (x))

typedef struct _DBL2VECT
{
    double x;
    double y;
} DBL2VECT;

typedef struct _DBL7VECT
{
    double x;
    double y;
    double intensity;
    double lambda1;
    double lambda2;
    DBL2VECT v1;
    DBL2VECT v2;
    DBL2VECT gradient;
} DBL7VECT;

typedef struct _INT3VECT
{
    int a;
    int b;
    int c;
    double intensity;
} INT3VECT;

enum BasisType { MQ, IMQ, Gaussian, TPS }; // type of basis functions

#endif
