#ifndef OTHERS_H
#define OTHERS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>

#include "vec4.h"
#include "mat4.h"
#include "tiny_obj_loader.h"
///-------------------------
///***** Image Type *****///
///-------------------------
typedef struct pixel_struct {
    unsigned char r, g, b;
}pixel_t;

typedef struct img_struct {
    pixel_t *data;
    int w, h;
} img_t;

typedef struct boundingbox{
    float maxx,maxy,minx,miny;
    std::vector<vec4> positions;
    std::vector<vec4> normals;
    float diffuse[3];
}bbox;

typedef struct allmat
{
    mat4 F,ViewO,ViewT; //F is projection matrix, ViewO and ViewT are view matrix
} mat;

mat readcam(const char *fname); //read camera.txt

typedef struct myface
{
    std::vector<unsigned int> indices;
}Face;

typedef struct mytdshapes
{
    std::vector<vec4> positions;
    std::vector<vec4> normals;
}tdshapes;


///-------------------------------------
///*** Matrix function and type *****///
///-------------------------------------

std::vector<Face> getface(std::vector<tinyobj::shape_t> shapes);

std::vector<tdshapes> rasterize(std::vector<Face> faces, std::vector<tinyobj::shape_t> shapes, mat cammat, int h, int w);
///---------------------------------------
///***** image processing Function*****///
///---------------------------------------

std::vector<bbox> bounding(std::vector<tdshapes> tdcoord, int w, int h, std::vector<tinyobj::shape_t> shapes, std::vector<tinyobj::material_t> materials); //compute 2-D bounding box and get material

img_t *scanline(img_t *img, std::vector<bbox> box, int h, int w, char *argv);

img_t *new_img(int w, int h);

void destroy_img(img_t **img);

void  write_ppm(const img_t *img, const char *fname);

#endif // OTHERS_H
