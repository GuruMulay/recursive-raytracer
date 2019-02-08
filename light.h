#ifndef LIGHT_H_INCLUDED
#define LIGHT_H_INCLUDED

#include <iostream>
#include "vector3d.h"

class Light{
    public:
        double w;
        Vector3D<double> position;       // position vector from origin to source point
        Vector3D<double> rgb;    // rgb value vector r = x; g = y; b = z

//////////        // store the x, y, z values of unit rays through the pixel plane
//////////        Matrix2D<double>* ptr_rayMatX;
//////////        Matrix2D<double>* ptr_rayMatY;
//////////        Matrix2D<double>* ptr_rayMatZ;

        Light();
        Light(double, double, double, double, double, double, double);

        void printLightAttributes();
        ~Light();
};

Light::Light(){
    w=1, position.x=1, position.y=1, position.z=1, rgb.x=1, rgb.y=0.5, rgb.z=0.5;
}

//Light::Light(double a, double b, double c, double d, double e, double f, double g) : w(a), position.x(b), position.y(c), position.z(d), rgb.x(e), rgb.y(f), rgb.z(g){}
Light::Light(double wVal, double posX, double posY, double posZ, double R, double G, double B) {
    w = wVal;
    position.x = posX;
    position.y = posY;
    position.z = posZ;
    rgb.x = R;
    rgb.y = G;
    rgb.z = B;
}


void Light::printLightAttributes(){
    std::cout << "Light Attributes -------------------------------------------------->" << std::endl;
    std::cout << position.x << " " << position.y << " " << position.z << " <----- position coordinates" << std::endl;
    std::cout << rgb.x << " " << rgb.y << " " << rgb.z << " <----- RGB values in rgb vector" << std::endl;
    std::cout << w << " <----- w of position vector" << std::endl;
    std::cout << "End ----------------------------------------------------------------" << '\n' << std::endl;

} // constructor for printing Light attributes

Light::~Light(){
}//Light object destructor


#endif // LIGHT_H_INCLUDED
