#ifndef POINT_H_INCLUDED
#define POINT_H_INCLUDED

#include <iostream>
#include "vector3d.h"

// Point class =========================================================
class Point{
public:
    double x, y, z;
    Point ();
    Point(double, double, double);
    double dist_0(void) { return sqrt(x*x + y*y + z*z); }
    void printPoint();
    Vector3D<double>* getUnitVectorAToB(Point, Point);

};

Point::Point(){
    x=0,y=0,z=0;
}
Point::Point(double a, double b, double c) : x(a), y(b), z(c){}

Vector3D<double>* Point::getUnitVectorAToB(Point A, Point B){
    Vector3D<double>* ptr_vec;
    ptr_vec = new Vector3D<double>(B.x-A.x, B.y-A.y, B.z-A.z);
    ptr_vec->normalize();
    return ptr_vec;
}

void Point::printPoint(){
    std::cout << "Point is ------------->" << std::endl;
    std::cout << x << " " << y << " " << z << " <----- x, y, z" << std::endl;
} // constructor for printing point


#endif // POINT_H_INCLUDED
