// The mat4 class is a 4x4 matrix containing a floating point number in each cell.
// Each mat4 is column-major, meaning the [] operator should return the Nth column of the matrix.
//
// OpenGL expects matrices to be column-major, so you will run into serious problems later
//  in the semester if your matrix library does using column-major.


#include <iostream>
#include "vec4.h"
#include "mat4.h"
#include <math.h>
#include <assert.h>

using namespace std;
///----------------------------------------------------------------------
/// Constructors
///----------------------------------------------------------------------

/// Default Constructor.  Initialize to identity matrix.
mat4::mat4()
{
    vec4 veca(1,0,0,0);
    vec4 vecb(0,1,0,0);
    vec4 vecc(0,0,1,0);
    vec4 vecd(0,0,0,1);
    data[0]=veca;
    data[1]=vecb;
    data[2]=vecc;
    data[3]=vecd;
}

/// Initializes the diagonal values of the matrix to diag. All other values are 0.
mat4::mat4(float diag)
{
    vec4 veca(diag,0,0,0);
    vec4 vecb(0,diag,0,0);
    vec4 vecc(0,0,diag,0);
    vec4 vecd(0,0,0,diag);
    data[0]=veca;
    data[1]=vecb;
    data[2]=vecc;
    data[3]=vecd;
}

/// Initializes matrix with each vector representing a column in the matrix
mat4::mat4(const vec4 &col0, const vec4 &col1, const vec4 &col2, const vec4& col3)
{
    data[0]=col0;
    data[1]=col1;
    data[2]=col2;
    data[3]=col3;
}

mat4::mat4(const mat4 &m2) // copy constructor; copies values of m2 into this
{
    *this=m2;
}

///----------------------------------------------------------------------
/// Getters
///----------------------------------------------------------------------

/// Returns the values of the column at the index
/// Does NOT check the array bound because doing so is slow
const vec4 &mat4:: operator[](unsigned int index) const
{
    return data[index];
}

/// Returns a reference to the column at the index
/// Does NOT check the array bound because doing so is slow
vec4 &mat4:: operator[](unsigned int index)
{
    return data[index];
}


/// Returns the values of the column at the index
/// DOES check the array bound because doing so is slow
const vec4 &mat4:: operator()(unsigned int index) const
{
    if (index>3)
    {
        cout<<"ATTENTION! Index out of matrix!!!"<<endl;
        assert(0);
        return data[index];
    }
    else
    {
        return data[index];
    }
}

/// Returns a reference to the column at the index
/// DOES check the array bound because doing so is slow
vec4 &mat4:: operator()(unsigned int index)
{
    if (index>3)
    {
        cout<<"ATTENTION! Index out of matrix!!!"<<endl;
        assert(0);
        return data[index];
    }
    else
    {
        return data[index];
    }
}

///----------------------------------------------------------------------
/// Static Initializers
///----------------------------------------------------------------------

/// Creates a 3-D rotation matrix.
/// Takes an angle in degrees and an axis represented by its xyz components, and outputs a 4x4 rotation matrix
mat4 mat4:: rot(float angle, float x, float y, float z)
{
    mat4 mat;
    mat[0][0]=cos(angle)+x*x*(1-cos(angle));
    mat[0][1]=z*sin(angle)+x*y*(1-cos(angle));
    mat[0][2]=-y*sin(angle)+x*z*(1-cos(angle));
    mat[0][3]=0;
    mat[1][0]=-z*sin(angle)+x*y*(1-cos(angle));
    mat[1][1]=cos(angle)+y*y*(1-cos(angle));
    mat[1][2]=x*sin(angle)+y*z*(1-cos(angle));
    mat[1][3]=0;
    mat[2][0]=y*sin(angle)+x*z*(1-cos(angle));
    mat[2][1]=-x*sin(angle)+y*z*(1-cos(angle));
    mat[2][2]=cos(angle)+z*z*(1-cos(angle));
    mat[2][3]=0;
    vec4 vec(0,0,0,1);
    mat[3]=vec;
    return mat;
}

/// Takes an xyz displacement and outputs a 4x4 translation matrix
mat4 mat4:: trans(float x, float y, float z)
{
    mat4 mat(1);
    mat[3][0]=x;
    mat[3][1]=y;
    mat[3][2]=z;
    return mat;
}

/// Takes an xyz scale and outputs a 4x4 scale matrix
mat4 mat4:: scale(float x, float y, float z)
{
    mat4 mat(1);
    mat[0][0]=x;
    mat[1][1]=y;
    mat[2][2]=z;
    return mat;
}


///----------------------------------------------------------------------
/// Operator Functions
///----------------------------------------------------------------------

/// Assign m2 to this and return this
mat4 &mat4:: operator=(const mat4 &m2)
{
    this->data[0]=m2[0];
    this->data[1]=m2[1];
    this->data[2]=m2[2];
    this->data[3]=m2[3];
    return *this;
}

/// Test for equality
bool mat4:: operator==(const mat4 &m2) const
{
    if (data[0]==m2[0] && data[1]==m2[1] && data[2]==m2[2] && data[3]==m2[3])
    {
        return true;
    }
    else
    {
        return false;
    }
}

/// Test for inequality
bool mat4:: operator!=(const mat4 &m2) const
{
    {
        if (data[0]==m2[0] && data[1]==m2[1] && data[2]==m2[2] && data[3]==m2[3])
        {
            return false;
        }
        else
        {
            return true;
        }
    }
}

/// Element-wise arithmetic
/// e.g. += adds the elements of m2 to this and returns this (like regular +=)
///      +  returns a new matrix whose elements are the sums of this and v2
mat4 &mat4:: operator+=(const mat4 &m2)
{
    this->data[0]+=m2[0];
    this->data[1]+=m2[1];
    this->data[2]+=m2[2];
    this->data[3]+=m2[3];
    return *this;
}

mat4 &mat4:: operator-=(const mat4 &m2)
{
    this->data[0]-=m2[0];
    this->data[1]-=m2[1];
    this->data[2]-=m2[2];
    this->data[3]-=m2[3];
    return *this;
}

mat4 &mat4:: operator*=(float c)                 // multiplication by a scalar
{
    this->data[0]=this->data[0]*c;
    this->data[1]=this->data[1]*c;
    this->data[2]=this->data[2]*c;
    this->data[3]=this->data[3]*c;
    return *this;
}

mat4 &mat4:: operator/=(float c)                 // division by a scalar
{
    this->data[0]=this->data[0]/c;
    this->data[1]=this->data[1]/c;
    this->data[2]=this->data[2]/c;
    this->data[3]=this->data[3]/c;
    return *this;
}

mat4  mat4:: operator+(const mat4 &m2) const
{
    mat4 mat;
    mat[0]=data[0]+m2[0];
    mat[1]=data[1]+m2[1];
    mat[2]=data[2]+m2[2];
    mat[3]=data[3]+m2[3];
    return mat;
}
mat4  mat4:: operator-(const mat4 &m2) const
{
    mat4 mat;
    mat[0]=data[0]-m2[0];
    mat[1]=data[1]-m2[1];
    mat[2]=data[2]-m2[2];
    mat[3]=data[3]-m2[3];
    return mat;
}
mat4  mat4:: operator*(float c) const             // multiplication by a scalar
{
    mat4 mat;
    mat[0]=data[0]*c;
    mat[1]=data[1]*c;
    mat[2]=data[2]*c;
    mat[3]=data[3]*c;
    return mat;
}
mat4  mat4:: operator/(float c) const             // division by a scalar
{
    mat4 mat;
    mat[0]=data[0]/c;
    mat[1]=data[1]/c;
    mat[2]=data[2]/c;
    mat[3]=data[3]/c;
    return mat;
}

/// Matrix multiplication (m1 * m2)
mat4 mat4:: operator*(const mat4 &m2) const
{
    mat4 mat;
    mat4 m1;
    m1=*this;
    m1=m1.transpose();
    vec4 veca(dot(m1[0],m2[0]),dot(m1[0],m2[1]),dot(m1[0],m2[2]),dot(m1[0],m2[3]));
    vec4 vecb(dot(m1[1],m2[0]),dot(m1[1],m2[1]),dot(m1[1],m2[2]),dot(m1[1],m2[3]));
    vec4 vecc(dot(m1[2],m2[0]),dot(m1[2],m2[1]),dot(m1[2],m2[2]),dot(m1[2],m2[3]));
    vec4 vecd(dot(m1[3],m2[0]),dot(m1[3],m2[1]),dot(m1[3],m2[2]),dot(m1[3],m2[3]));
    mat[0]=veca;
    mat[1]=vecb;
    mat[2]=vecc;
    mat[3]=vecd;
    return mat.transpose();
}

/// Matrix/vector multiplication (m * v)
/// Assume v is a column vector (ie. a 4x1 matrix)
vec4 mat4:: operator*(const vec4 &v) const
{
    mat4 m1;
    m1=*this;
    m1=m1.transpose();
    vec4 vec(dot(m1[0],v),dot(m1[1],v),dot(m1[2],v),dot(m1[3],v));
    return vec;
}

///----------------------------------------------------------------------
/// Matrix Operations
///----------------------------------------------------------------------

/// Returns the transpose of the input matrix (v_ij == v_ji)
mat4 mat4:: transpose() const
{
    mat4 mat;
    for (int i=0;i<4;i++)
    {
        for (int j=0;j<4;j++)
        {
            mat[i][j]=this->data[j][i];
        }
    }
    return mat;
}

/// Returns a column of the input matrix
const vec4 &mat4:: col(unsigned int index) const  // const version
{
    return data[index];
}
vec4 &mat4:: col(unsigned int index)             // non-const version
{
    return data[index];
}

///----------------------------------------------------------------------
/// Other Functions (not part of the mat4 class)
///----------------------------------------------------------------------

/// Scalar multiplication (c * m)
mat4 operator*(float c, const mat4 &m)
{
    mat4 mat;
    mat[0]=m[0]*c;
    mat[1]=m[1]*c;
    mat[2]=m[2]*c;
    mat[3]=m[3]*c;
    return mat;
}

/// Vector/matrix multiplication (v * m)
/// Assume v is a row vector (ie. a 1x4 matrix)
vec4 operator*(const vec4 &v, const mat4 &m)
{
    vec4 vec;
    vec[0]=dot(v,m[0]);
    vec[1]=dot(v,m[1]);
    vec[2]=dot(v,m[2]);
    vec[3]=dot(v,m[3]);
    return vec;
}

/// Prints the matrix to a stream in a nice format
std::ostream &operator<<(std::ostream &o, const mat4 &m)
{
    mat4 mt=m.transpose();
    o << mt[0] <<endl;
    o << mt[1] <<endl;
    o << mt[2] <<endl;
    o << mt[3] <<endl;
    return o;
}

