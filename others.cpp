#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include <math.h>
#include "others.h"
#include "vec4.h"
#include "mat4.h"
using namespace std;


///-------------------------------------
///*** Functions and type *****///
///-------------------------------------
mat readcam(const char *fname)
{
    float l,r,t,b,n,f;
    float eye_x,eye_y,eye_z;
    float center_x,center_y,center_z;
    float up_x,up_y,up_z;

    mat cammat;

    assert(fname != NULL);

    FILE *fc = fopen(fname, "rb");
    assert(fc != NULL);

    fscanf(fc,"%f %f %f %f %f %f",&l,&r,&b,&t,&n,&f);
    fscanf(fc,"%f %f %f", &eye_x,&eye_y,&eye_z);
    fscanf(fc,"%f %f %f", &center_x,&center_y,&center_z);
    fscanf(fc,"%f %f %f", &up_x,&up_y,&up_z);

    cammat.F[0][0]=2*n/(r-l);
    cammat.F[1][1]=2*n/(t-b);
    cammat.F[2][0]=(r+l)/(r-l);
    cammat.F[2][1]=(t+b)/(t-b);
    cammat.F[2][2]=f/(f-n);
    cammat.F[2][3]=1;
    cammat.F[3][2]=-f*n/(f-n);
    cammat.F[3][3]=0;

    vec4 fvec(center_x-eye_x,center_y-eye_y,center_z-eye_z,float(0));
    fvec.norm();
    vec4 uvec(up_x,up_y,up_z,0);
    vec4 rvec=cross(fvec,uvec);
    vec4 lvec(0,0,0,1);


    cammat.ViewO=cammat.ViewO.transpose();
    cammat.ViewO(0)=rvec;
    cammat.ViewO(1)=uvec;
    cammat.ViewO(2)=fvec;
    cammat.ViewO(3)=lvec;
    cammat.ViewO=cammat.ViewO.transpose();

    vec4 trans(-eye_x,-eye_y,-eye_z,1);
    cammat.ViewT(3)=trans;
    return cammat;
}

std::vector<Face> getface(std::vector<tinyobj::shape_t> shapes)
{
    std::vector<Face> faces;
    for (unsigned int i=0;i<shapes[0].mesh.indices.size() ;i=i+3)
    {
        Face f;
        f.indices.push_back(shapes[0].mesh.indices[i]);
        f.indices.push_back(shapes[0].mesh.indices[i+1]);
        f.indices.push_back(shapes[0].mesh.indices[i+2]);
        faces.push_back(f);
    }
    return faces;
}

std::vector<tdshapes> rasterize(std::vector<Face> faces, std::vector<tinyobj::shape_t> shapes, mat cammat, int h, int w)
{
    std::vector<tdshapes> tdcoord;
    std::vector<tdshapes> origtdcoord;
    for (unsigned int i=0;i<faces.size();i++)
    {
        tdshapes myshape;
        for (unsigned j=0;j<3;j++)
        {
            float x=shapes[0].mesh.positions[faces[i].indices[j]*3];
            float y=shapes[0].mesh.positions[faces[i].indices[j]*3+1];
            float z=shapes[0].mesh.positions[faces[i].indices[j]*3+2];
            float nx=shapes[0].mesh.normals[faces[i].indices[j]*3];
            float ny=shapes[0].mesh.normals[faces[i].indices[j]*3+1];
            float nz=shapes[0].mesh.normals[faces[i].indices[j]*3+2];

            vec4 vec(x,y,z,1);
            vec4 nvec(nx,ny,nz,0);
            vec=cammat.F*cammat.ViewO*cammat.ViewT*vec;
            vec=vec/vec[3];
            vec[0]=(vec[0]+1)*w/2;                 //convert to pixel coordinate
            vec[1]=(1-vec[1])*h/2;
            nvec=cammat.ViewO*nvec;

            myshape.positions.push_back(vec);
            myshape.normals.push_back(nvec);
        }
        tdcoord.push_back(myshape);
    }
    return tdcoord;
}

///-------------------------------
///***** image processing *****///
///-------------------------------
///

std::vector<bbox> bounding(std::vector<tdshapes> tdcoord, int w, int h, std::vector<tinyobj::shape_t> shapes, std::vector<tinyobj::material_t> materials)
{
    std::vector<bbox> Boundbox;
    for (unsigned int i=0;i<tdcoord.size();i++)
    {
        bbox boundbox;
        boundbox.maxx=max(tdcoord[i].positions[0][0],tdcoord[i].positions[1][0]);
        boundbox.maxx=max(boundbox.maxx,tdcoord[i].positions[2][0]);
        boundbox.minx=min(tdcoord[i].positions[0][0],tdcoord[i].positions[1][0]);
        boundbox.minx=min(boundbox.minx,tdcoord[i].positions[2][0]);
        boundbox.maxy=max(tdcoord[i].positions[0][1],tdcoord[i].positions[1][1]);
        boundbox.maxy=max(boundbox.maxy,tdcoord[i].positions[2][1]);
        boundbox.miny=min(tdcoord[i].positions[0][1],tdcoord[i].positions[1][1]);
        boundbox.miny=min(boundbox.miny,tdcoord[i].positions[2][1]);
        if (boundbox.maxx<0)        //if out of pixel
        {
            continue;
        }
        else if (boundbox.minx>w)   //if out of pixel
        {
            continue;
        }
        else if (boundbox.maxy<0)   //if out of pixel
        {
            continue;
        }
        else if (boundbox.miny>h)   //if out of pixel
        {
            continue;
        }
        else                        //clamp
        {
            boundbox.minx=max((float)0,boundbox.minx);
            boundbox.maxx=min((float)w,boundbox.maxx);
            boundbox.miny=max((float)0,boundbox.miny);
            boundbox.maxy=min((float)w,boundbox.maxy);
        }
        boundbox.positions.push_back(tdcoord[i].positions[0]);
        boundbox.positions.push_back(tdcoord[i].positions[1]);
        boundbox.positions.push_back(tdcoord[i].positions[2]);
        boundbox.diffuse[0]=materials[shapes[0].mesh.material_ids[i]].diffuse[0];
        boundbox.diffuse[1]=materials[shapes[0].mesh.material_ids[i]].diffuse[1];
        boundbox.diffuse[2]=materials[shapes[0].mesh.material_ids[i]].diffuse[2];
        boundbox.normals.push_back(tdcoord[i].normals[0]);
        boundbox.normals.push_back(tdcoord[i].normals[1]);
        boundbox.normals.push_back(tdcoord[i].normals[2]);
        Boundbox.push_back(boundbox);
    }
    return Boundbox;
}

img_t *scanline(img_t *img, std::vector<bbox> box, int h, int w, char *argv)
{
    float zbuffer[w*h];
    std::fill_n(zbuffer, w*h, 2); //fill in zbuffer with 2
    for (unsigned int i=0;i<box.size();i++)
    {
        float pz;
        float tol=1;
        float area=(box[i].positions[2][0]-box[i].positions[0][0])*(box[i].positions[1][1]-box[i].positions[0][1])-(box[i].positions[2][1]-box[i].positions[0][1])*(box[i].positions[1][0]-box[i].positions[0][0]);
        for (int j=floor(box[i].miny);j<ceil(box[i].maxy);j++)
        {
            for (int k=floor(box[i].minx);k<ceil(box[i].maxx);k++)
            {
                float px=k+0.5;
                float py=j+0.5;
                // edge function to check if a point is in triangle
                float Ea=(px-box[i].positions[0][0])*(box[i].positions[1][1]-box[i].positions[0][1])-(py-box[i].positions[0][1])*(box[i].positions[1][0]-box[i].positions[0][0]);
                float Eb=(px-box[i].positions[1][0])*(box[i].positions[2][1]-box[i].positions[1][1])-(py-box[i].positions[1][1])*(box[i].positions[2][0]-box[i].positions[1][0]);
                float Ec=(px-box[i].positions[2][0])*(box[i].positions[0][1]-box[i].positions[2][1])-(py-box[i].positions[2][1])*(box[i].positions[0][0]-box[i].positions[2][0]);
                // coefficients for barycentric
                float lambda2=Ea/area;
                float lambda0=Eb/area;
                float lambda1=Ec/area;
                if (Ea>=-tol && Eb>=-tol && Ec>=-tol) // if a point is in triangle
                {
                    pz=1/(lambda0/box[i].positions[0][2]+lambda1/box[i].positions[1][2]+lambda2/box[i].positions[2][2]); //get z value
                    if (pz<zbuffer[j*img->w+k]&&pz<1&&pz>0)
                    {
                        zbuffer[j*img->w+k]=pz;  //update zbuffer
                        if (argv==NULL)  //no input, diffuse color
                        {
                            img->data[j*img->w+k].r=round(box[i].diffuse[0]*255);
                            img->data[j*img->w+k].g=round(box[i].diffuse[1]*255);
                            img->data[j*img->w+k].b=round(box[i].diffuse[2]*255);
                        }
                        else if (strcmp(argv,"--white")==0)
                        {
                            img->data[j*img->w+k].r=255;
                            img->data[j*img->w+k].g=255;
                            img->data[j*img->w+k].b=255;
                        }
                        else if (strcmp(argv,"--norm_flat")==0)  //use first indices normal to get color
                        {
                            img->data[j*img->w+k].r=round((box[i].normals[0][0]+1)*255/2);
                            img->data[j*img->w+k].g=round((box[i].normals[0][1]+1)*255/2);
                            img->data[j*img->w+k].b=round((box[i].normals[0][2]+1)*255/2);
                        }
                        else if (strcmp(argv,"--norm_bary")==0) //barycentric interpolation without perspective correct
                        {
                            img->data[j*img->w+k].r=round(lambda2*(box[i].normals[2][0]+1)*255/2+lambda0*(box[i].normals[0][0]+1)*255/2+lambda1*(box[i].normals[1][0]+1)*255/2);
                            img->data[j*img->w+k].g=round(lambda2*(box[i].normals[2][1]+1)*255/2+lambda0*(box[i].normals[0][1]+1)*255/2+lambda1*(box[i].normals[1][1]+1)*255/2);
                            img->data[j*img->w+k].b=round(lambda2*(box[i].normals[2][2]+1)*255/2+lambda0*(box[i].normals[0][2]+1)*255/2+lambda1*(box[i].normals[1][2]+1)*255/2);
                        }
                        else if (strcmp(argv,"--norm_bary_z")==0) //barycentric interpolation with perspective correct
                        {
                            img->data[j*img->w+k].r=round((pz*(lambda2*box[i].normals[2][0]/box[i].positions[2][2]+lambda0*box[i].normals[0][0]/box[i].positions[0][2]+lambda1*box[i].normals[1][0]/box[i].positions[1][2])+1)*255/2);
                            img->data[j*img->w+k].g=round((pz*(lambda2*box[i].normals[2][1]/box[i].positions[2][2]+lambda0*box[i].normals[0][1]/box[i].positions[0][2]+lambda1*box[i].normals[1][1]/box[i].positions[1][2])+1)*255/2);
                            img->data[j*img->w+k].b=round((pz*(lambda2*box[i].normals[2][2]/box[i].positions[2][2]+lambda0*box[i].normals[0][2]/box[i].positions[0][2]+lambda1*box[i].normals[1][2]/box[i].positions[1][2])+1)*255/2);
                        }
                        else if (strcmp(argv,"--norm_gouraud")==0)  //gouraud shading wihout perspective correct
                        {
                            // slope for three lines
                            float m0=(box[i].positions[0][1]-box[i].positions[1][1])/(box[i].positions[0][0]-box[i].positions[1][0]);
                            float m1=(box[i].positions[1][1]-box[i].positions[2][1])/(box[i].positions[1][0]-box[i].positions[2][0]);
                            float m2=(box[i].positions[2][1]-box[i].positions[0][1])/(box[i].positions[2][0]-box[i].positions[0][0]);
                            // three intersection points
                            float x0=box[i].positions[0][0]+((j+0.5)-box[i].positions[0][1])/m0;
                            float x1=box[i].positions[1][0]+((j+0.5)-box[i].positions[1][1])/m1;
                            float x2=box[i].positions[2][0]+((j+0.5)-box[i].positions[2][1])/m2;
                            // sort their x coordinates from small to big
                            float array[4]={px,x0,x1,x2};
                            float inter=0;
                            for (int n=0;n<3;n++)  //sort array from small to big
                            {
                                for (int m=1;m<4-n;m++)
                                {
                                    if (array[n]>array[n+m])
                                    {
                                        inter=array[n+m];
                                        array[n+m]=array[n];
                                        array[n]=inter;
                                    }
                                }
                            }
                            float lx,rx;            //x coordinates for L and R points
                            float llambda,rlambda;  //coefficients for interpolating color of L and R points
                            vec4 lnormal,rnormal;   //interpolating results for L and R points
                            for (int n=1;n<3;n++)
                            {
                                if (array[n]==px) //find which element is point P, and left point is L, right point is R
                                {
                                    if (array[n-1]==x0) //left point for px is x0
                                    {
                                        lx=x0;
                                        llambda=(py-box[i].positions[0][1])/(box[i].positions[1][1]-box[i].positions[0][1]);
                                        lnormal=llambda*box[i].normals[1]+(1-llambda)*box[i].normals[0];
                                        if (array[n+1]==x1)
                                        {
                                            rx=x1;
                                            rlambda=(py-box[i].positions[1][1])/(box[i].positions[2][1]-box[i].positions[1][1]);
                                            rnormal=rlambda*box[i].normals[2]+(1-rlambda)*box[i].normals[1];
                                        }
                                        else if (array[n+1]==x2)
                                        {
                                            rx=x2;
                                            rlambda=(py-box[i].positions[2][1])/(box[i].positions[0][1]-box[i].positions[2][1]);
                                            rnormal=rlambda*box[i].normals[0]+(1-rlambda)*box[i].normals[2];
                                        }
                                    }
                                    else if (array[n-1]==x1)
                                    {
                                        lx=x1;
                                        llambda=(py-box[i].positions[1][1])/(box[i].positions[2][1]-box[i].positions[1][1]);
                                        lnormal=llambda*box[i].normals[2]+(1-llambda)*box[i].normals[1];
                                        if (array[n+1]==x0)
                                        {
                                            rx=x0;
                                            rlambda=(py-box[i].positions[0][1])/(box[i].positions[1][1]-box[i].positions[0][1]);
                                            rnormal=rlambda*box[i].normals[1]+(1-rlambda)*box[i].normals[0];
                                        }
                                        else if (array[n+1]==x2)
                                        {
                                            rx=x2;
                                            rlambda=(py-box[i].positions[2][1])/(box[i].positions[0][1]-box[i].positions[2][1]);
                                            rnormal=rlambda*box[i].normals[0]+(1-rlambda)*box[i].normals[2];
                                        }
                                    }
                                    else if (array[n-1]==x2)
                                    {
                                        lx=x2;
                                        llambda=(py-box[i].positions[2][1])/(box[i].positions[0][1]-box[i].positions[2][1]);
                                        lnormal=llambda*box[i].normals[0]+(1-llambda)*box[i].normals[2];
                                        if (array[n+1]==x0)
                                        {
                                            rx=x0;
                                            rlambda=(py-box[i].positions[0][1])/(box[i].positions[1][1]-box[i].positions[0][1]);
                                            rnormal=rlambda*box[i].normals[1]+(1-rlambda)*box[i].normals[0];
                                        }
                                        else if (array[n+1]==x1)
                                        {
                                            rx=x1;
                                            rlambda=(py-box[i].positions[1][1])/(box[i].positions[2][1]-box[i].positions[1][1]);
                                            rnormal=rlambda*box[i].normals[2]+(1-rlambda)*box[i].normals[1];
                                        }
                                    }
                                    break;
                                }
                            }
                            float lambda=(px-lx)/(rx-lx);
                            vec4 pnormal=lambda*rnormal+(1-lambda)*lnormal;
                            img->data[j*img->w+k].r=(pnormal[0]+1)*255/2;
                            img->data[j*img->w+k].g=(pnormal[1]+1)*255/2;
                            img->data[j*img->w+k].b=(pnormal[2]+1)*255/2;
                        }
                        else if (strcmp(argv,"--norm_gouraud_z")==0)
                        {
                            float m0=(box[i].positions[0][1]-box[i].positions[1][1])/(box[i].positions[0][0]-box[i].positions[1][0]);
                            float m1=(box[i].positions[1][1]-box[i].positions[2][1])/(box[i].positions[1][0]-box[i].positions[2][0]);
                            float m2=(box[i].positions[2][1]-box[i].positions[0][1])/(box[i].positions[2][0]-box[i].positions[0][0]);
                            // three intersection points
                            float x0=box[i].positions[0][0]+((j+0.5)-box[i].positions[0][1])/m0;
                            float x1=box[i].positions[1][0]+((j+0.5)-box[i].positions[1][1])/m1;
                            float x2=box[i].positions[2][0]+((j+0.5)-box[i].positions[2][1])/m2;
                            // sort their x coordinates from small to big
                            float array[4]={px,x0,x1,x2};
                            float inter=0;
                            for (int n=0;n<3;n++) //sort array from small to big
                            {
                                for (int m=1;m<4-n;m++)
                                {
                                    if (array[n]>array[n+m])
                                    {
                                        inter=array[n+m];
                                        array[n+m]=array[n];
                                        array[n]=inter;
                                    }
                                }
                            }
                            float lx,rx;            //x coordinates for L and R points
                            float lpz,rpz;          //perspective correct z coordinates for L and R
                            float llambda,rlambda;  //coefficients for interpolating color of L and R points
                            vec4 lnormal,rnormal;   //interpolating results for L and R points
                            for (int n=1;n<3;n++)
                            {
                                if (array[n]==px) //find which element is point P, and left point is L, right point is R
                                {
                                    if (array[n-1]==x0) //left point for px is x0
                                    {
                                        lx=x0;
                                        llambda=(py-box[i].positions[0][1])/(box[i].positions[1][1]-box[i].positions[0][1]);
                                        lpz=1/(llambda/box[i].positions[1][2]+(1-llambda)/box[i].positions[0][2]);  //correct z coord for L
                                        lnormal=llambda*box[i].normals[1]+(1-llambda)*box[i].normals[0];
                                        if (array[n+1]==x1)
                                        {
                                            rx=x1;
                                            rlambda=(py-box[i].positions[1][1])/(box[i].positions[2][1]-box[i].positions[1][1]);
                                            rpz=1/(rlambda/box[i].positions[2][2]+(1-rlambda)/box[i].positions[1][2]); //correct z coord for R
                                            rnormal=rlambda*box[i].normals[2]+(1-rlambda)*box[i].normals[1];
                                        }
                                        else if (array[n+1]==x2)
                                        {
                                            rx=x2;
                                            rlambda=(py-box[i].positions[2][1])/(box[i].positions[0][1]-box[i].positions[2][1]);
                                            rpz=1/(rlambda/box[i].positions[0][2]+(1-rlambda)/box[i].positions[2][2]); //correct z coord for R
                                            rnormal=rlambda*box[i].normals[0]+(1-rlambda)*box[i].normals[2];
                                        }
                                    }
                                    else if (array[n-1]==x1)
                                    {
                                        lx=x1;
                                        llambda=(py-box[i].positions[1][1])/(box[i].positions[2][1]-box[i].positions[1][1]);
                                        lpz=1/(llambda/box[i].positions[2][2]+(1-llambda)/box[i].positions[1][2]);
                                        lnormal=llambda*box[i].normals[2]+(1-llambda)*box[i].normals[1];
                                        if (array[n+1]==x0)
                                        {
                                            rx=x0;
                                            rlambda=(py-box[i].positions[0][1])/(box[i].positions[1][1]-box[i].positions[0][1]);
                                            rpz=1/(rlambda/box[i].positions[1][2]+(1-rlambda)/box[i].positions[0][2]);
                                            rnormal=rlambda*box[i].normals[1]+(1-rlambda)*box[i].normals[0];
                                        }
                                        else if (array[n+1]==x2)
                                        {
                                            rx=x2;
                                            rlambda=(py-box[i].positions[2][1])/(box[i].positions[0][1]-box[i].positions[2][1]);
                                            rpz=1/(rlambda/box[i].positions[0][2]+(1-rlambda)/box[i].positions[2][2]);
                                            rnormal=rlambda*box[i].normals[0]+(1-rlambda)*box[i].normals[2];
                                        }
                                    }
                                    else if (array[n-1]==x2)
                                    {
                                        lx=x2;
                                        llambda=(py-box[i].positions[2][1])/(box[i].positions[0][1]-box[i].positions[2][1]);
                                        lpz=1/(llambda/box[i].positions[0][2]+(1-llambda)/box[i].positions[2][2]);
                                        lnormal=llambda*box[i].normals[0]+(1-llambda)*box[i].normals[2];
                                        if (array[n+1]==x0)
                                        {
                                            rx=x0;
                                            rlambda=(py-box[i].positions[0][1])/(box[i].positions[1][1]-box[i].positions[0][1]);
                                            rpz=1/(rlambda/box[i].positions[1][2]+(1-rlambda)/box[i].positions[0][2]);
                                            rnormal=rlambda*box[i].normals[1]+(1-rlambda)*box[i].normals[0];
                                        }
                                        else if (array[n+1]==x1)
                                        {
                                            rx=x1;
                                            rlambda=(py-box[i].positions[1][1])/(box[i].positions[2][1]-box[i].positions[1][1]);
                                            rpz=1/(rlambda/box[i].positions[2][2]+(1-rlambda)/box[i].positions[1][2]);
                                            rnormal=rlambda*box[i].normals[2]+(1-rlambda)*box[i].normals[1];
                                        }
                                    }
                                    break;
                                }
                            }
                            float lambda=(px-lx)/(rx-lx);
                            vec4 pnormal=pz*(lambda*rnormal/rpz+(1-lambda)*lnormal/lpz);  //perspective correct normal
                            img->data[j*img->w+k].r=(pnormal[0]+1)*255/2;
                            img->data[j*img->w+k].g=(pnormal[1]+1)*255/2;
                            img->data[j*img->w+k].b=(pnormal[2]+1)*255/2;
                        }
                        else if (strcmp(argv,"--lambertian")==0)
                        {
                            float m0=(box[i].positions[0][1]-box[i].positions[1][1])/(box[i].positions[0][0]-box[i].positions[1][0]);
                            float m1=(box[i].positions[1][1]-box[i].positions[2][1])/(box[i].positions[1][0]-box[i].positions[2][0]);
                            float m2=(box[i].positions[2][1]-box[i].positions[0][1])/(box[i].positions[2][0]-box[i].positions[0][0]);
                            // three intersection points
                            float x0=box[i].positions[0][0]+((j+0.5)-box[i].positions[0][1])/m0;
                            float x1=box[i].positions[1][0]+((j+0.5)-box[i].positions[1][1])/m1;
                            float x2=box[i].positions[2][0]+((j+0.5)-box[i].positions[2][1])/m2;
                            // sort their x coordinates from small to big
                            float array[4]={px,x0,x1,x2};
                            float inter=0;
                            for (int n=0;n<3;n++) //sort array from small to big
                            {
                                for (int m=1;m<4-n;m++)
                                {
                                    if (array[n]>array[n+m])
                                    {
                                        inter=array[n+m];
                                        array[n+m]=array[n];
                                        array[n]=inter;
                                    }
                                }
                            }
                            vec4 lv2,rv2,fv2;
                            float ltheta,rtheta,ftheta;
                            float ldot,rdot,fdot;
                            float lx,rx;            //x coordinates for L and R points
                            float llambda,rlambda;  //coefficients for interpolating color of L and R points
                            vec4 lnormal,rnormal;   //interpolating results for L and R points
                            for (int n=1;n<3;n++)
                            {
                                if (array[n]==px) //find which element is point P, and left point is L, right point is R
                                {
                                    if (array[n-1]==x0) //left point for px is x0
                                    {
                                        lx=x0;
                                        llambda=(py-box[i].positions[0][1])/(box[i].positions[1][1]-box[i].positions[0][1]);
                                        box[i].normals[0].norm();
                                        box[i].normals[1].norm();
                                        ldot=dot(box[i].normals[0],box[i].normals[1]);
                                        if (ldot < 0)
                                        {
                                            box[i].normals[1]=-1*box[i].normals[1];
                                            ldot = -ldot;
                                        }


                                        ldot=max((float)-1,min(ldot,(float)1));
                                        ltheta=acos(ldot)*llambda;
                                        lv2=box[i].normals[1]-box[i].normals[0]*ldot;
                                        lv2.norm();
                                        lnormal=box[i].normals[0]*cos(ltheta)+lv2*sin(ltheta);
                                        if (ldot>0.9995)
                                        {
                                            lnormal=llambda*box[i].normals[1]+(1-llambda)*box[i].normals[0];
                                        }
                                        if (array[n+1]==x1)
                                        {
                                            rx=x1;
                                            rlambda=(py-box[i].positions[1][1])/(box[i].positions[2][1]-box[i].positions[1][1]);
                                            box[i].normals[1].norm();
                                            box[i].normals[2].norm();
                                            rdot=dot(box[i].normals[1],box[i].normals[2]);
                                            if (rdot < 0)
                                            {
                                                box[i].normals[2]=-1*box[i].normals[2];
                                                rdot = -rdot;
                                            }


                                            rdot=max((float)-1,min(rdot,(float)1));
                                            rtheta=acos(rdot)*rlambda;
                                            rv2=box[i].normals[2]-box[i].normals[1]*rdot;
                                            rv2.norm();
                                            rnormal=box[i].normals[1]*cos(rtheta)+rv2*sin(rtheta);
                                            if (rdot>0.9995)
                                            {
                                                rnormal=rlambda*box[i].normals[2]+(1-rlambda)*box[i].normals[1];
                                            }
                                        }
                                        else if (array[n+1]==x2)
                                        {
                                            rx=x2;
                                            rlambda=(py-box[i].positions[2][1])/(box[i].positions[0][1]-box[i].positions[2][1]);
                                            box[i].normals[2].norm();
                                            box[i].normals[0].norm();
                                            rdot=dot(box[i].normals[2],box[i].normals[0]);
                                            if (rdot < 0)
                                            {
                                                box[i].normals[0]=-1*box[i].normals[0];
                                                rdot = -rdot;
                                            }


                                            rdot=max((float)-1,min(rdot,(float)1));
                                            rtheta=acos(rdot)*rlambda;
                                            rv2=box[i].normals[0]-box[i].normals[2]*rdot;
                                            rv2.norm();
                                            rnormal=box[i].normals[2]*cos(rtheta)+rv2*sin(rtheta);
                                            if (rdot>0.9995)
                                            {
                                                rnormal=rlambda*box[i].normals[0]+(1-rlambda)*box[i].normals[2];
                                            }
                                        }
                                    }
                                    else if (array[n-1]==x1)
                                    {
                                        lx=x1;
                                        llambda=(py-box[i].positions[1][1])/(box[i].positions[2][1]-box[i].positions[1][1]);
                                        box[i].normals[1].norm();
                                        box[i].normals[2].norm();
                                        ldot=dot(box[i].normals[1],box[i].normals[2]);
                                        if (ldot < 0)
                                        {
                                            box[i].normals[2]=-1*box[i].normals[2];
                                            ldot = -ldot;
                                        }


                                        ldot=max((float)-1,min(ldot,(float)1));
                                        ltheta=acos(ldot)*llambda;
                                        lv2=box[i].normals[2]-box[i].normals[1]*ldot;
                                        lv2.norm();
                                        lnormal=box[i].normals[1]*cos(ltheta)+lv2*sin(ltheta);
                                        if (ldot>0.9995)
                                        {
                                            lnormal=llambda*box[i].normals[2]+(1-llambda)*box[i].normals[1];
                                        }
                                        if (array[n+1]==x0)
                                        {
                                            rx=x0;
                                            rlambda=(py-box[i].positions[0][1])/(box[i].positions[1][1]-box[i].positions[0][1]);
                                            box[i].normals[0].norm();
                                            box[i].normals[1].norm();
                                            rdot=dot(box[i].normals[0],box[i].normals[1]);
                                            if (rdot < 0)
                                            {
                                                box[i].normals[1]=-1*box[i].normals[1];
                                                rdot = -rdot;
                                            }


                                            rdot=max((float)-1,min(rdot,(float)1));
                                            rtheta=acos(rdot)*rlambda;
                                            rv2=box[i].normals[1]-box[i].normals[0]*rdot;
                                            rv2.norm();
                                            rnormal=box[i].normals[0]*cos(rtheta)+rv2*sin(rtheta);
                                            if (rdot>0.9995)
                                            {
                                                rnormal=rlambda*box[i].normals[1]+(1-rlambda)*box[i].normals[0];
                                            }
                                        }
                                        else if (array[n+1]==x2)
                                        {
                                            rx=x2;
                                            rlambda=(py-box[i].positions[2][1])/(box[i].positions[0][1]-box[i].positions[2][1]);
                                            box[i].normals[2].norm();
                                            box[i].normals[0].norm();
                                            rdot=dot(box[i].normals[2],box[i].normals[0]);
                                            if (rdot < 0)
                                            {
                                                box[i].normals[0]=-1*box[i].normals[0];
                                                rdot = -rdot;
                                            }


                                            rdot=max((float)-1,min(rdot,(float)1));
                                            rtheta=acos(rdot)*rlambda;
                                            rv2=box[i].normals[0]-box[i].normals[2]*rdot;
                                            rv2.norm();
                                            rnormal=box[i].normals[2]*cos(rtheta)+rv2*sin(rtheta);
                                            if (rdot>0.9995)
                                            {
                                                rnormal=rlambda*box[i].normals[0]+(1-rlambda)*box[i].normals[2];
                                            }
                                        }
                                    }
                                    else if (array[n-1]==x2)
                                    {
                                        lx=x2;
                                        llambda=(py-box[i].positions[2][1])/(box[i].positions[0][1]-box[i].positions[2][1]);
                                        box[i].normals[2].norm();
                                        box[i].normals[0].norm();
                                        ldot=dot(box[i].normals[2],box[i].normals[0]);
                                        if (ldot < 0)
                                        {
                                            box[i].normals[0]=-1*box[i].normals[0];
                                            ldot = -ldot;
                                        }


                                        ldot=max((float)-1,min(ldot,(float)1));
                                        ltheta=acos(ldot)*llambda;
                                        lv2=box[i].normals[0]-box[i].normals[2]*ldot;
                                        lv2.norm();
                                        lnormal=box[i].normals[2]*cos(ltheta)+lv2*sin(ltheta);
                                        if (ldot>0.9995)
                                        {
                                            lnormal=llambda*box[i].normals[0]+(1-llambda)*box[i].normals[2];
                                        }
                                        if (array[n+1]==x0)
                                        {
                                            rx=x0;
                                            rlambda=(py-box[i].positions[0][1])/(box[i].positions[1][1]-box[i].positions[0][1]);
                                            box[i].normals[0].norm();
                                            box[i].normals[1].norm();
                                            rdot=dot(box[i].normals[0],box[i].normals[1]);
                                            if (rdot < 0)
                                            {
                                                box[i].normals[1]=-1*box[i].normals[1];
                                                rdot = -rdot;
                                            }

                                            rdot=max((float)-1,min(rdot,(float)1));
                                            rtheta=acos(rdot)*rlambda;
                                            rv2=box[i].normals[1]-box[i].normals[0]*rdot;
                                            rv2.norm();
                                            rnormal=box[i].normals[0]*cos(rtheta)+rv2*sin(rtheta);
                                            if (rdot>0.9995)
                                            {
                                                rnormal=rlambda*box[i].normals[1]+(1-rlambda)*box[i].normals[0];
                                            }
                                        }
                                        else if (array[n+1]==x1)
                                        {
                                            rx=x1;
                                            rlambda=(py-box[i].positions[1][1])/(box[i].positions[2][1]-box[i].positions[1][1]);
                                            box[i].normals[1].norm();
                                            box[i].normals[2].norm();
                                            rdot=dot(box[i].normals[1],box[i].normals[2]);
                                            if (rdot < 0)
                                            {
                                                box[i].normals[2]=-1*box[i].normals[2];
                                                rdot = -rdot;
                                            }
                                            rdot=max((float)-1,min(rdot,(float)1));
                                            rtheta=acos(rdot)*rlambda;
                                            rv2=box[i].normals[2]-box[i].normals[1]*rdot;
                                            rv2.norm();
                                            rnormal=box[i].normals[1]*cos(rtheta)+rv2*sin(rtheta);
                                            if (rdot>0.9995)
                                            {
                                                rnormal=rlambda*box[i].normals[2]+(1-rlambda)*box[i].normals[1];
                                            }
                                        }
                                    }
                                    break;
                                }
                            }
                            float lambda=(px-lx)/(rx-lx);
                            lnormal.norm();
                            rnormal.norm();
                            fdot=dot(lnormal,rnormal);
                            if (fdot < 0)
                            {
                                lnormal=-1*lnormal;
                                fdot = -fdot;
                            }
                            fdot=max((float)-1,min(fdot,(float)1));
                            ftheta=acos(fdot)*lambda;
                            fv2=rnormal-lnormal*fdot;
                            fv2.norm();
                            vec4 pnormal=lnormal*cos(ftheta)+fv2*sin(ftheta);
                            pnormal.norm();
                            img->data[j*img->w+k].r=(pnormal[0]+1)*255/2;
                            img->data[j*img->w+k].g=(pnormal[1]+1)*255/2;
                            img->data[j*img->w+k].b=(pnormal[2]+1)*255/2;
                        }
                        else
                        {
                            cout << "Wrong Parameter Input!" << endl;
                            assert(0);
                        }
                    }
                }
            }
        }
    }
    return img;
}

img_t *new_img(int w, int h) {
    assert(w > 0);
    assert(h > 0);

    img_t *img = (img_t *) malloc(sizeof(img_t));

    img->w = w;
    img->h = h;

    img->data = (pixel_t *) malloc(w * h * sizeof(pixel_t));

    memset(img->data, 0, w * h * sizeof(pixel_t));

    return img;
}

void destroy_img(img_t **img) {
    free((*img)->data);
    (*img)->data = NULL;
    free(*img);
    *img = NULL;
}

void write_ppm(const img_t *img, const char *fname) {
    assert(img != NULL); // crash if img is NULL
    assert(fname != NULL); // crash if fname is NULL

    FILE *f = fopen(fname, "wb"); // open fname for writing in binary mode; clobbers any existing file
    assert(f != NULL); // crash if file did not open

    fprintf(f, "P6\n%d %d 255\n", img->w, img->h); // write the image header
    fwrite(img->data, img->w * img->h, 3, f); // write the image data

    fclose(f); // close the file
}
