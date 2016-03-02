//
//  TriMesh.h
//  arbf-image-super_resolution
//
//  Created by Ke Liu on 2/26/15.
//  Copyright (c) 2015 Ke Liu. All rights reserved.
//

#ifndef __arbf_image_super_resolution__TriMesh__
#define __arbf_image_super_resolution__TriMesh__

#include "Common.h"

class TriMesh
{
public:
    /*
     * Default constructor
     */
    TriMesh();
    
    /*
     * Create mesh by reading mesh file.
     * \param path: mesh file path
     */
    TriMesh(const char *path);
    
    /*
     * Copy constructor
     */
    TriMesh(const TriMesh &other);
    
    /*
     * Destructor
     */
    virtual ~TriMesh();
    
    /*
     * Assignment operator
     */
    TriMesh& operator=(const TriMesh &other);
    
    int numVertices() const;
    void setNumVertices(int nv);
    
    int numFaces() const;
    void setNumFaces(int nf);
    
    DBL7VECT* getVertices();
    const DBL7VECT* getVertices() const;
    
    INT3VECT* getFaces();
    const INT3VECT* getFaces() const;
    
    int read(const char *path);
    
private:
    DBL7VECT *m_vertex; // vertices
    INT3VECT *m_face; // faces
    int m_nv; // # of vertices
    int m_nf; // # of faces
    
private:
    static const int BUF_SIZE; // size of buffer when reading TriMesh
};

#endif /* defined(__arbf_image_super_resolution__TriMesh__) */
