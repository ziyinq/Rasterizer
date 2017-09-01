#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <others.h>
#include <mat4.h>

#include "tiny_obj_loader.h"
#include "others.h"

using namespace std;

int main(int argc, char *argv[])
{
    if (strcmp(argv[1],"rasterize")==0)
    {
        int w=atoi(argv[4]);
        int h=atoi(argv[5]);
        if (w*h>520000)
        {
            cout << "Warning! Increasing w or h may crash the program!" << endl;   //zbuffer is an 1D array and has a limit
        }
        std::vector<tinyobj::shape_t> shapes;
        std::vector<tinyobj::material_t> materials;
        const char *filename= argv[2];
        const char *mtl_basepath = NULL;

        string a=tinyobj::LoadObj(shapes, materials, filename , mtl_basepath);

        mat cammat=readcam(argv[3]); //call function to return F matrix, View matrices
        if (a.empty()) //check if the obj file name is correct
        {
            std::vector<Face> faces=getface(shapes);
            std::vector<tdshapes> tdcoord=rasterize(faces,shapes,cammat,h,w);
            std::vector<bbox> box=bounding(tdcoord,w,h,shapes,materials);
            img_t *img=new_img(w, h);
            img=scanline(img, box, h, w, argv[7]);
            write_ppm(img, argv[6]);
            destroy_img(&img);
        }
        else
        {
            cout << a << endl;
            assert(a.empty()!=0);
        }
    }
    else
    {
        cout << "Wrong Input! Should be rasterize!" <<endl;
    }
    return 0;
}
