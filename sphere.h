#ifndef SPHERE_H_INCLUDED
#define SPHERE_H_INCLUDED

#include <iostream>
#include "vector3d.h"

class Sphere{
    public:
        double radius;
        Vector3D<double> center;  // center vector from origin to source point
        //~ Vector3D<double> rgb;    // rgb value vector r = x; g = y; b = z

        Vector3D<double> ka;    // ambient
        Vector3D<double> kd;    // diffuse
        Vector3D<double> ks;    // specular (phong)
        Vector3D<double> kr;    // attenuation multiplier!
        double rho;

        Sphere();
        Sphere(double, double, double, double, double, double, double);

        void printSphereAttributes();
        ~Sphere();
};

Sphere::Sphere(){
    radius=1, center.x=1, center.y=1, center.z=1, ka.x=1, ka.y=1, ka.z=1;
}

//Sphere::Sphere(double a, double b, double c, double d, double e, double f, double g) : radius(a), center.x(b), center.y(c), center.z(d), rgb.x(e), rgb.y(f), rgb.z(g){}
Sphere::Sphere(double rVal, double posX, double posY, double posZ, double R, double G, double B){
    radius = rVal;
    center.x = posX;
    center.y = posY;
    center.z = posZ;
    ka.x = R;
    ka.y = G;
    ka.z = B;
}

void Sphere::printSphereAttributes(){
    std::cout << "Sphere Attributes -------------------------------------------------->" << std::endl;
    std::cout << center.x << " " << center.y << " " << center.z << " <----- center coordinates" << std::endl;
    //~ std::cout << rgb.x << " " << rgb.y << " " << rgb.z << " <----- RGB values in rgb vector" << std::endl;
    std::cout << ka.x << " " << ka.y << " " << ka.z << " <----- Ambient" << std::endl;
    std::cout << kd.x << " " << kd.y << " " << kd.z << " <----- Diffuse" << std::endl;
    std::cout << ks.x << " " << ks.y << " " << ks.z << " <----- Specular" << std::endl;
    std::cout << kr.x << " " << kr.y << " " << kr.z << " <----- Attenuation Multiplier" << std::endl;
    std::cout << radius << " <----- radius of the sphere" << std::endl;
    std::cout << rho << " <----- rho of the sphere" << std::endl;
    std::cout << "End ----------------------------------------------------------------" << '\n' << std::endl;

} // constructor for printing Sphere attributes

Sphere::~Sphere(){
} // Sphere object destructor

#endif // SPHERE_H_INCLUDED
