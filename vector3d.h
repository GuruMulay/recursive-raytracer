//===================================================================================
// Name        : vector3d.h
// Created on  : 15 October 2016
// Description : Class for 3D vector operations
//===================================================================================
#ifndef VECTOR3D_H_INCLUDED
#define VECTOR3D_H_INCLUDED

#include <iostream>
#include <math.h>

//#include "vector3d.cpp"

template <class T>
class Vector3D{
public:
    T x, y, z;
    T magnitude;

    Vector3D(){
        x=0, y=0, z=0;
    };

    //Vector3D(T, T, T);
	Vector3D(T x_val, T y_val, T z_val){
		x = x_val;
		y = y_val;
		z = z_val;
	} // Vector3D constructor

  // const std::string empty = std::string();
	virtual ~Vector3D(); // destructor Vector3D

	// General operations -----------------------------------------------------------------------------
    bool isZero(); // check if a zero vector
	T getMagnitude(); // magnitude/norm of the vector
	void normalize(); // normalize the vector
    bool isEqual(const Vector3D<T> &v_a) const; // comparison of two vectors

    // Dot and Cross Products of two vectors ----------------------------------------------------------
	T dotProduct(Vector3D<T> &v_a); // dot product of the vector with v_a
	Vector3D<T> crossProduct(Vector3D<T> &v_a); // cross product of the vector with v_a
	Vector3D<T> pairwiseProduct(Vector3D<T> &v_a); // pairwise product of the vector with v_a

	// Assignment Operation -----------------------------------------
    Vector3D<T>& operator=(const Vector3D<T> &v_a);

	// Arithmatic Operations on Vectors ---------------------------------------------------------------
    Vector3D<T> operator+(const Vector3D<T> &v_a) const; // addition of the vector with v_a
	Vector3D<T> operator-(const Vector3D<T> &v_a) const; // subtraction of v_a from the vector
	Vector3D<T> operator*(const T alpha) const; // multiplication by a scalar alpha
	Vector3D<T> operator/(const T alpha) const; // division by a scalar alpha

    // print vector's x, y, z components --------------------------------------------------------------
    void printVector3D(std::string str);
};

template<class T>
bool Vector3D<T>::isZero(){
    return (x==0 && y==0 && z==0);
} // check if a zero vector

template<class T>
T Vector3D<T>::getMagnitude(){
    magnitude = sqrt(x*x + y*y + z*z);
	return magnitude;
} // magnitude/norm of the vector

template<class T>
void Vector3D<T>::normalize(){
    T norm = sqrt(x*x + y*y + z*z);
	x = x/norm;
	y = y/norm;
	z = z/norm;
} // normalize the vector

template<class T>
T Vector3D<T>::dotProduct(Vector3D<T> &v_a){
    if (this->getMagnitude()==0 || v_a.getMagnitude()==0){std::cout << "------------------ WARNING: v1 or v2 is a zero vector! -------------------" << std::endl;}
	return (x*v_a.x + y*v_a.y + z*v_a.z);
} // dot product of the vector with v_a

template<class T>
Vector3D<T> Vector3D<T>::crossProduct(Vector3D<T> &v_a){
	if (this->getMagnitude()==0 || v_a.getMagnitude()==0){std::cout << "------------------ WARNING: v1 or v2 is a zero vector! -------------------" << std::endl;}
	return Vector3D<T>(y*v_a.z - z*v_a.y, z*v_a.x - x*v_a.z, x*v_a.y - y*v_a.x);
} // cross product of the vector with v_a

template<class T>
Vector3D<T> Vector3D<T>::pairwiseProduct(Vector3D<T> &v_a){
	return Vector3D<T>(x*v_a.x, y*v_a.y, z*v_a.z);
} // pairwise product of the vector with v_a

template<class T>
bool Vector3D<T>::isEqual(const Vector3D<T> &v_a) const{
    return (x==v_a.x && y==v_a.y && z==v_a.z);
} // comparison of two vectors


// Assignment Operation
template<class T>
Vector3D<T>& Vector3D<T>::operator=(const Vector3D<T> &v_a){
    this->x = v_a.x;
    this->y = v_a.y;
    this->z = v_a.z;
    return *this;
}

template<class T>
Vector3D<T> Vector3D<T>::operator+(const Vector3D<T> &v_a) const{
    return Vector3D<T>(x+v_a.x, y+v_a.y, z+v_a.z);
} // addition of the vector with v_a

template<class T>
Vector3D<T> Vector3D<T>::operator-(const Vector3D<T> &v_a) const{
    return Vector3D<T>(x-v_a.x, y-v_a.y, z-v_a.z);
} // subtraction of v_a from the vector

template<class T>
Vector3D<T> Vector3D<T>::operator*(const T alpha) const{
    return Vector3D<T>(alpha*x, alpha*y, alpha*z);
} // multiplication by a scalar alpha

template<class T>
Vector3D<T> Vector3D<T>::operator/(const T alpha) const{
    if (alpha==0){std::cout << "------------------ WARNING: Division by 0! alpha is 0! -------------------" << std::endl;}
    return Vector3D<T>(x/alpha, y/alpha, z/alpha);
} // division by a scalar alpha

template<class T>
void Vector3D<T>::printVector3D(std::string str){
    //~ std::cout <<  << std::endl;
    std::cout << str << " " << x << " " << y << " " << z << " (x, y, z) " << std::endl;
} // print vector's x, y, z components

template<class T>
Vector3D<T>::~Vector3D(){
} // Vector3D object destructor


#endif // VECTOR3D_H_INCLUDED


//===================================================================================
//Referred Library  : http://www.gamedev.net/ topic/ 487050-c-vector-class/
//             : http://web.engr.oregonstate.edu/ ~tgd/classes/ 533/code/lib/ vector.h
//===================================================================================
