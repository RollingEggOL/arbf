/*
 * mesh_parser.c
 *
 *  Created on: Jan 24, 2014
 *      Author: keliu
 */

/******************************************************************
 To use this reader, you need to make sure your input mesh
 has the keyword "OFF" in the beginning, followed by the mesh.
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <memory.h>

#define BUF_SIZE (256)

typedef struct
{
	double x;
	double y;
} DBL2VECT;

typedef struct
{
	double x;
	double y;
	double intensity;
	double lambda1;
	double lambda2;
	DBL2VECT v1;
	DBL2VECT v2;
} DBL7VECT;

typedef struct
{
	int a;
	int b;
	int c;
	double intensity;
} INT3VECT;

typedef struct
{
	int nv; /* number of vertices */
	int nf; /* number of faces */
	DBL7VECT *vertex;
	INT3VECT *face;
} SurFaceMesh;

void read_mesh_tri(const char *file_name, SurFaceMesh **surfmesh)
{
	int n, m, dim;
	int a, b, c, d;
	double x, y, lambda1, lambda2, vx1, vy1, z;
	char line[BUF_SIZE];
	FILE *fin;
    
#ifdef _WIN32
	if (fopen_s(&fin, file_name, "r"))
#else
	if ((fin = fopen(file_name, "r")) == NULL)
#endif
	{
		fprintf(stderr, "ERROR: reading mesh file ... \n");
		exit(EXIT_FAILURE);
	};
    
	/* TRI format */
	while (fgets(line, BUF_SIZE, fin) != NULL)
	{
		if (line[0] == 'T' && line[1] == 'R' && line[2] == 'I')
		{
			break;
		}
        
		if (line[0] == '#')
		{
			break;
		}
	}

#ifdef _WIN32
	fscanf_s(fin, "%d\n", &dim);
#else
    fscanf(fin, "%d\n", &dim);
#endif
    if(dim != 2) // support 2D TRI
    {
        fprintf(stderr, "ERROR: corrupted TRI format ... \n");
    }
    
#ifdef _WIN32
	fscanf_s(fin, "%d %d\n", &m, &n); // # of vertices, # of faces
#else
	fscanf(fin, "%d %d\n", &m, &n); // # of vertices, # of faces
#endif
	(*surfmesh) = (SurFaceMesh*) malloc(sizeof(SurFaceMesh));
	(*surfmesh)->nv = m;
	(*surfmesh)->nf = n;
	(*surfmesh)->vertex = (DBL7VECT*) malloc(sizeof(DBL7VECT) * (*surfmesh)->nv);
	(*surfmesh)->face = (INT3VECT*) malloc(sizeof(INT3VECT) * (*surfmesh)->nf);
    
	/* read vertices */
	for (n = 0; n < (*surfmesh)->nv; n++)
	{
#ifdef _WIN32
		fscanf_s(fin, "%lf %lf %lf %lf %lf %lf\n", &x, &y, &lambda1, &lambda2, &vx1, &vy1);
#else
		fscanf(fin, "%lf %lf %lf %lf %lf %lf\n", &x, &y, &lambda1, &lambda2, &vx1, &vy1);
#endif
		(*surfmesh)->vertex[n].x = x;
		(*surfmesh)->vertex[n].y = y;
        (*surfmesh)->vertex[n].intensity = 0;
		(*surfmesh)->vertex[n].lambda1 = lambda1;
		(*surfmesh)->vertex[n].lambda2 = lambda2;
		(*surfmesh)->vertex[n].v1.x = vx1;
		(*surfmesh)->vertex[n].v1.y = vy1;
		(*surfmesh)->vertex[n].v2.x = -vy1;
		(*surfmesh)->vertex[n].v2.y = vx1;
	}
    
	/* read faces */
	for (n = 0; n < (*surfmesh)->nf; n++)
	{
#ifdef _WIN32
		fscanf_s(fin, "%d %d %d %d %lf\n", &a, &b, &c, &d, &z);
#else
		fscanf(fin,"%d %d %d %d %lf\n",&a, &b, &c, &d, &z);
#endif
		(*surfmesh)->face[n].a = b;
		(*surfmesh)->face[n].b = c;
		(*surfmesh)->face[n].c = d;
        (*surfmesh)->face[n].intensity = z;
	}
    
	fclose(fin);
}

void clean_mesh(SurFaceMesh **mesh)
{
	if ((*mesh))
	{
		if ((*mesh)->vertex)
		{
			free((*mesh)->vertex);
		}
        
		if ((*mesh)->face)
		{
			free((*mesh)->face);
		}

		free(*mesh);
	}
	*mesh = NULL;
}
