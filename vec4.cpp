#include <iostream>
#include <math.h>
#include "vec4.h"
#include <assert.h>
using namespace std;

///----------------------------------------------------------------------
/// Constructors
///----------------------------------------------------------------------

vec4::vec4() // initialize the vector to (0, 0, 0, 0)
{
    data[0]=0;
    data[1]=0;
    data[2]=0;
    data[3]=0;
}

vec4::vec4(float x, float y, float z, float w)
{
    data[0]=x;
    data[1]=y;
    data[2]=z;
    data[3]=w;
}

vec4::vec4(const vec4 &v2) // copy constructor, copies values of v2 into this
{
    this->data[0]=v2[0];
    this->data[1]=v2[1];
    this->data[2]=v2[2];
    this->data[3]=v2[3];
}

///----------------------------------------------------------------------
/// Getters/Setters
///----------------------------------------------------------------------

/// Returns the value at index
/// Does NOT check the array bound because doing so is slow
float vec4::operator [](unsigned int index) const
{
    return data[index];
}

/// Returns a reference to the value at index
/// Does NOT check the array bound because doing so is slow
float &vec4::operator[](unsigned int index)
{
    return data[index];
}

/// Returns the value at index
/// DOES check the array bound with assert even though is slow
float vec4:: operator()(unsigned int index) const
{
    if (index>3)
    {
        cout<<"ATTENTION, Index out of array!!!"<<endl;
        assert(0);
        return data[index];
    }
    else
    {
        return data[index];
    }
}

/// Returns a reference to the value at index
/// DOES check the array bound with assert even though is slow
float &vec4:: operator()(unsigned int index)
{
    if (index>3)
    {
        cout<<"ATTENTION, Index out of array!!!"<<endl;
        assert(0);
        return data[index];
    }
    else
    {
        return this->data[index];
    }
}

///----------------------------------------------------------------------
/// Operator Methods
///----------------------------------------------------------------------

/// Assign v2 to this and return a reference to this
vec4 &vec4:: operator=(const vec4 &v2)
{
    this->data[0]=v2[0];
    this->data[1]=v2[1];
    this->data[2]=v2[2];
    this->data[3]=v2[3];
    return *this;
}

/// Test for equality
bool vec4:: operator==(const vec4 &v2) const	   //Component-wise comparison
{
    if (data[0]==v2[0] && data[1]==v2[1] && data[2]==v2[2] && data[3]==v2[3])
    {
        return true;
    }
    else
    {
        return false;
    }
}

/// Test for inequality
bool vec4:: operator!=(const vec4 &v2) const	   //Component-wise comparison
{
    if (data[0]==v2[0] && data[1]==v2[1] && data[2]==v2[2] && data[3]==v2[3])
    {
        return false;
    }
    else
    {
        return true;
    }
}

/// Arithmetic:
/// e.g. += adds v2 to this and return this (like regular +=)
///      +  returns a new vector that is sum of this and v2
vec4 &vec4:: operator+=(const vec4 &v2)
{
    this->data[0]=this->data[0]+v2[0];
    this->data[1]=this->data[1]+v2[1];
    this->data[2]=this->data[2]+v2[2];
    this->data[3]=this->data[3]+v2[3];
    return *this;
}
vec4 &vec4:: operator-=(const vec4 &v2)
{
    this->data[0]=this->data[0]-v2[0];
    this->data[1]=this->data[1]-v2[1];
    this->data[2]=this->data[2]-v2[2];
    this->data[3]=this->data[3]-v2[3];
    return *this;
}
vec4 &vec4:: operator*=(float c)                 // multiplication by a scalar
{
    this->data[0]=this->data[0]*c;
    this->data[1]=this->data[1]*c;
    this->data[2]=this->data[2]*c;
    this->data[3]=this->data[3]*c;
    return *this;
}

vec4 &vec4:: operator/=(float c)                 // division by a scalar
{
    this->data[0]=this->data[0]/c;
    this->data[1]=this->data[1]/c;
    this->data[2]=this->data[2]/c;
    this->data[3]=this->data[3]/c;
    return *this;
}

vec4 vec4:: operator+(const vec4 &v2) const
{
    vec4 vec;
    vec[0]=data[0]+v2[0];
    vec[1]=data[1]+v2[1];
    vec[2]=data[2]+v2[2];
    vec[3]=data[3]+v2[3];
    return vec;
}

vec4 vec4:: operator-(const vec4 &v2) const
{
    vec4 vec;
    vec[0]=data[0]-v2[0];
    vec[1]=data[1]-v2[1];
    vec[2]=data[2]-v2[2];
    vec[3]=data[3]-v2[3];
    return vec;
}
vec4 vec4:: operator*(float c) const             // multiplication by a scalar
{
    vec4 vec;
    vec[0]=data[0]*c;
    vec[1]=data[1]*c;
    vec[2]=data[2]*c;
    vec[3]=data[3]*c;
    return vec;
}
vec4 vec4:: operator/(float c) const             // division by a scalar
{
    vec4 vec;
    vec[0]=data[0]/c;
    vec[1]=data[1]/c;
    vec[2]=data[2]/c;
    vec[3]=data[3]/c;
    return vec;
}

///----------------------------------------------------------------------
/// Other Methods
///----------------------------------------------------------------------


/// return a new vec4 that is a normalized (unit-length) version of this one
vec4 vec4:: normalize() const
{
    vec4 vec;
    vec=*this;
    float a=length(vec);
    vec[0]=vec[0]/a;
    vec[1]=vec[1]/a;
    vec[2]=vec[2]/a;
    vec[3]=vec[3]/a;
    return vec;
}

/// noralize this vector in place
void vec4:: norm()
{
    float a=sqrt(data[0]*data[0]+data[1]*data[1]+data[2]*data[2]+data[3]*data[3]);
    data[0]=data[0]/a;
    data[1]=data[1]/a;
    data[2]=data[2]/a;
    data[3]=data[3]/a;
}


///----------------------------------------------------------------------
/// Other Functions (not part of the vec4 class)
///----------------------------------------------------------------------

/// Returns the geometric length of the input vector
float length(const vec4 &v)
{
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
}

/// Dot Product
float dot(const vec4 &v1, const vec4 &v2)
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3];
}

/// Cross Product
vec4 cross(const vec4 &v1, const vec4 &v2)      //Compute the result of v1 x v2 using only their X, Y, and Z elements.
{                                               //In other words, treat v1 and v2 as 3D vectors, not 4D vectors.
    vec4 vec;                                   //The fourth element of the resultant vector should be 0.
    vec[0]=v1[1]*v2[2]-v1[2]*v2[1];
    vec[1]=v1[2]*v2[0]-v1[0]*v2[2];
    vec[2]=v1[0]*v2[1]-v1[1]*v2[0];
    vec[3]=0;
    return vec;
}

/// Scalar Multiplication (c * v)
vec4 operator*(float c, const vec4 &v)
{
    vec4 vec;
    vec[0]=c*v[0];
    vec[1]=c*v[1];
    vec[2]=c*v[2];
    vec[3]=c*v[3];
    return vec;
}

/// Prints the vector to a stream in a nice format for integration with cout
std::ostream &operator<<(std::ostream &o, const vec4 &v)
{
    o<<"["<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<"]";
    return o;
}
