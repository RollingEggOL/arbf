//
//  main.cpp
//  generate-mesh
//
//  Created by Ke Liu on 3/10/14.
//  Copyright (c) 2014 Ke Liu. All rights reserved.
//

#include <cstdlib>
#include <cstdio>
#include <string>

using namespace std;

#define XDIM (256)
#define YDIM (256)

const string ROOT = "/Users/keliu/Documents/projects/arbf/generate-mesh/";
char data[XDIM][YDIM];
int x1, y1, x2, y2;

// vertices
const int NUM_VERTICES = 14;
const int NUM_TRIANGLES = 22;
int v0[2], v1[2], v2[2], v3[2], v4[2], v5[2], v6[2], v7[2], v8[2], v9[2],
    v10[2], v11[2], v12[2], v13[2];

void set_data()
{
    x1 = XDIM / 4;
    y1 = YDIM / 4;
    x2 = 3 * x1;
    y2 = 3 * y1;

    for (int j = 0; j < YDIM; ++j)
    {
        for (int i = 0; i < XDIM; ++i)
        {
            if (i >= x1 && i <= x2 &&
                j >= y1 && j <= y2)
            {
                data[i][j] = (char) 0;
            }
            else
            {
                data[i][j] = (char) 255;
            }
        }
    }
}

void write_PPM(const char* file_name)
{
    FILE* fp;
    
    if ((fp = fopen(file_name, "wb")) == NULL)
    {
        fprintf(stderr, "ERROR: write PPM ...\n");
        exit(EXIT_FAILURE);
    }
    
    printf("Writing PPM ... \n");
    
    /* write to a PPM image */
    fprintf(fp, "P6\n%d %d\n%d\n", XDIM, YDIM, 255);
    
    for (int j = 0; j < YDIM; j++)
    {
        for (int i = 0; i < XDIM; i++)
        {
            fputc(data[i][j], fp);
            fputc(data[i][j], fp);
            fputc(data[i][j], fp);
        }
    }
    
    fclose(fp);
}

void write_OFF(const char* file_name)
{
    v0[0] = 0; v0[1] = 0;
    v1[0] = XDIM - 1; v1[1] = 0;
    v2[0] = XDIM - 1; v2[1] = YDIM - 1;
    v3[0] = 0; v3[1] = YDIM - 1;
    v4[0] = x1; v4[1] = y1;
    v5[0] = x2; v5[1] = y1;
    v6[0] = x2; v6[1] = y2;
    v7[0] = x1; v7[1] = y2;
    v8[0] = (int) (v4[0] + v5[0] + v7[0]) / 3;
    v8[1] = (int) (v4[1] + v5[1] + v7[1]) / 3;
    v9[0] = (int) (v5[0] + v6[0] + v7[0]) / 3;
    v9[1] = (int) (v5[1] + v6[1] + v7[1]) / 3;
    v10[0] = (int) (v4[0] + v8[0] + v7[0]) / 3;
    v10[1] = (int) (v4[1] + v8[1] + v7[1]) / 3;
    v11[0] = (int) (v4[0] + v8[0] + v5[0]) / 3;
    v11[1] = (int) (v4[1] + v8[1] + v5[1]) / 3;
    v12[0] = (int) (v5[0] + v6[0] + v9[0]) / 3;
    v12[1] = (int) (v5[1] + v6[1] + v9[1]) / 3;
    v13[0] = (int) (v6[0] + v7[0] + v9[0]) / 3;
    v13[1] = (int) (v6[1] + v7[1] + v9[1]) / 3;
    
    FILE* fp;
    if ((fp = fopen(file_name, "w")) == NULL)
    {
        fprintf(stderr, "ERROR: write OFF ...\n");
        exit(EXIT_FAILURE);
    }
    
    printf("Writing OFF ... \n");
    
    /* write to a OFF image */
    fprintf(fp, "OFF\n%d %d %d\n\n", NUM_VERTICES, NUM_TRIANGLES, 0);
    
    fprintf(fp, "%d %d\n", v0[0], v0[1]);
    fprintf(fp, "%d %d\n", v1[0], v1[1]);
    fprintf(fp, "%d %d\n", v2[0], v2[1]);
    fprintf(fp, "%d %d\n", v3[0], v3[1]);
    fprintf(fp, "%d %d\n", v4[0], v4[1]);
    fprintf(fp, "%d %d\n", v5[0], v5[1]);
    fprintf(fp, "%d %d\n", v6[0], v6[1]);
    fprintf(fp, "%d %d\n", v7[0], v7[1]);
    fprintf(fp, "%d %d\n", v8[0], v8[1]);
    fprintf(fp, "%d %d\n", v9[0], v9[1]);
    fprintf(fp, "%d %d\n", v10[0], v10[1]);
    fprintf(fp, "%d %d\n", v11[0], v11[1]);
    fprintf(fp, "%d %d\n", v12[0], v12[1]);
    fprintf(fp, "%d %d\n", v13[0], v13[1]);
    fprintf(fp, "\n");
    
    fprintf(fp, "3 %d %d %d\n", 0, 4, 1);
    fprintf(fp, "3 %d %d %d\n", 1, 4, 5);
    fprintf(fp, "3 %d %d %d\n", 2, 1, 5);
    fprintf(fp, "3 %d %d %d\n", 5, 6, 2);
    fprintf(fp, "3 %d %d %d\n", 2, 6, 3);
    fprintf(fp, "3 %d %d %d\n", 6, 7, 3);
    fprintf(fp, "3 %d %d %d\n", 0, 3, 7);
    fprintf(fp, "3 %d %d %d\n", 7, 4, 0);
    
    fprintf(fp, "3 %d %d %d\n", 7, 5, 8);
    fprintf(fp, "3 %d %d %d\n", 9, 5, 7);
    
    fprintf(fp, "3 %d %d %d\n", 10, 4, 7);
    fprintf(fp, "3 %d %d %d\n", 7, 8, 10);
    fprintf(fp, "3 %d %d %d\n", 10, 8, 4);
    fprintf(fp, "3 %d %d %d\n", 4, 8, 11);
    fprintf(fp, "3 %d %d %d\n", 11, 5, 4);
    fprintf(fp, "3 %d %d %d\n", 11, 8, 5);
    fprintf(fp, "3 %d %d %d\n", 12, 6, 5);
    fprintf(fp, "3 %d %d %d\n", 5, 9, 12);
    fprintf(fp, "3 %d %d %d\n", 12, 9, 6);
    fprintf(fp, "3 %d %d %d\n", 13, 7, 6);
    fprintf(fp, "3 %d %d %d\n", 6, 9, 13);
    fprintf(fp, "3 %d %d %d\n", 13, 9, 7);
    
    fclose(fp);
}

int main(int argc, const char * argv[])
{
    set_data();
    write_PPM(string(ROOT + "simple.ppm").c_str());
    write_OFF(string(ROOT + "simple.off").c_str());
    
    return 0;
}

