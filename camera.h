#ifndef CAMERA_H_INCLUDED
#define CAMERA_H_INCLUDED

#include <iostream>
#include <iomanip>  // std::setprecision(2)

#include <fstream>
//#include <sstream>
//#include <string>

#include "vector3d.h"
#include "matrix2d.h"
#include "ray.h"

class Camera{
    public:
//        double focal_pt_x, focal_pt_y, focal_pt_z;
//        double lookat_pt_x, lookat_pt_y, lookat_pt_z;
//        double up_x, up_y, up_z;

        double focalLength;
        Vector3D<double> eye;       // vector from origin to eye point
        Vector3D<double> lookat;    // vector from origin to lookat point
        Vector3D<double> up;        // must be a unit vector pointing in up direction
        //Vector3D<double>* eye_ptr;

        //ImagePlane* ptr_img;

        double left, bottom, right, top;
        unsigned res_u, res_v;  // res_u along u (horizontal column)

        //Ray* ray;
        // store the x, y, z values of unit rays through the pixel plane
        Matrix2D<double>* ptr_rayMatX;
        Matrix2D<double>* ptr_rayMatY;
        Matrix2D<double>* ptr_rayMatZ;
//        Matrix2D<double> matX;
//        Matrix2D<double> matY;
//        Matrix2D<double> matZ;
        // store the x, y, z values of pixel plane
        Matrix2D<double>* ptr_pixMatX;
        Matrix2D<double>* ptr_pixMatY;
        Matrix2D<double>* ptr_pixMatZ;

        // stores the depth values for each pixel
        Matrix2D<double>* ptr_depthMat;

        Vector3D<double> center_vector; // vector from origin to center of image plane
        // UNIT vectors u, v, w in the rotation matrix for the camera
        Vector3D<double> u_vector;
        Vector3D<double> v_vector;
        Vector3D<double> w_vector;

        Vector3D<double> pixel_pt; // vector from origin to pixel point
        Vector3D<double> d_vector; // UNIT vector from eye point to pixel point

        //int res_u, res_v;  // created a ptr to image class instead
        Camera();
        Camera(double, double, double, double, double, unsigned, unsigned);

        void getWVector();
        void getUVector();
        void getVVector();

        void getImagePlaneCenterVector();
        void getPixeledImagePlane();
        void printCameraAttributes();
        ~Camera();
};
Camera::Camera(){
    left=-1, bottom=-1, right=1, top=1, focalLength=2, res_u=256, res_v=256;
}
Camera::Camera(double a, double b, double c, double d, double e, unsigned ru, unsigned rv) : left(a), bottom(b), right(c), top(d), focalLength(e), res_u(ru), res_v(rv){}

void Camera::getWVector(){
    w_vector = eye - lookat; //vector subtraction
    w_vector.normalize();
//    std::cout << "+++++++++ Printing w vector +++++++++++" << std::endl;
//    w_vector.printVector3D();
}

void Camera::getUVector(){
    up.normalize();
    // special condition when up vector == w vector
    if (up.isEqual(w_vector) || up.isEqual(w_vector*(-1))){
        std:: cout << "Be Careful! UP and W point in the same direction !!!!!Equal!!!!!" << '\n' << std::endl;
        //w_vector.x = w_vector.x - 0.001; // changes w a little bit!
        w_vector.y = w_vector.y - 0.001;
    }
    u_vector = up.crossProduct(w_vector);
    u_vector.normalize();
//    std::cout << "+++++++++ Printing u vector +++++++++++" << std::endl;
//    u_vector.printVector3D();
}

void Camera::getVVector(){
    v_vector = w_vector.crossProduct(u_vector);
//    std::cout << "+++++++++ Printing v vector +++++++++++" << std::endl;
//    v_vector.printVector3D();
}

void Camera::getImagePlaneCenterVector(){
    center_vector = eye + ((w_vector)*(-focalLength));
//    std::cout << "+++++++++ Printing center vector +++++++++++" << std::endl;
//    center_vector.printVector3D();
} // get the center vector for the image plane

void Camera::getPixeledImagePlane(){

    double pxU, pyV;

    ptr_pixMatX = new Matrix2D<double>(res_v, res_u, 1);
    ptr_pixMatY = new Matrix2D<double>(res_v, res_u, 1);
    ptr_pixMatZ = new Matrix2D<double>(res_v, res_u, 1);

    ptr_rayMatX = new Matrix2D<double>(res_v, res_u, 1);
    ptr_rayMatY = new Matrix2D<double>(res_v, res_u, 1);
    ptr_rayMatZ = new Matrix2D<double>(res_v, res_u, 1);

//    ptr_matX->printMatrix();
//    ptr_matX->setElement(res_v-1, res_u-1, 55);
//    ptr_matX->printMatrix();

    //~ std::ofstream opstr1, opstr2;
    //~ opstr1.open("imagePlanePixelCoordinates.txt", std::ofstream::out);
    //~ opstr2.open("imagePlaneUnitRay.txt", std::ofstream::out);

    //opstr << "coordinates of image plane pixels";

    // started from top->bottom and left->right
    for (unsigned j=0; j<res_v; j++){
        pyV = (double)j/(res_v - 1) * (top - bottom) + bottom;  // at least one (int) should be typecasted to (double)
        for (unsigned i=0; i<res_u; i++){
            pxU = (double)i/(res_u - 1) * (right - left) + left;
            pixel_pt = center_vector + (u_vector * pxU) + (v_vector * pyV); // the original vectors for pixel in the image plane
            //~ opstr1 << "(" << std::setprecision(3) << std::fixed << pixel_pt.x << "," << std::setprecision(3) << std::fixed << pixel_pt.y << "," << std::setprecision(3) << std::fixed << pixel_pt.z << ") ";

            // Store these pixel positions in matrices corresponding to x, y, and z values
            // j rows i columns
            ptr_pixMatX->setElement(j, i, pixel_pt.x);
            ptr_pixMatY->setElement(j, i, pixel_pt.y);
            ptr_pixMatZ->setElement(j, i, pixel_pt.z);

            d_vector = pixel_pt - eye; // It's not unit yet as required
            d_vector.normalize();
            //~ opstr2 << "(" << std::setprecision(3) << std::fixed << d_vector.x << "," << std::setprecision(3) << std::fixed << d_vector.y << "," << std::setprecision(3) << std::fixed << d_vector.z << ") ";

            // Store these unit vectors in matrices corresponding to x, y, and z values of unit vectors values dx, dy, and dz
            ptr_rayMatX->setElement(j, i, d_vector.x);
            ptr_rayMatY->setElement(j, i, d_vector.y);
            ptr_rayMatZ->setElement(j, i, d_vector.z);
            // Stored x, y, and z values in individual matrices
        }
        //~ opstr1 << '\n';
        //~ opstr2 << '\n';
    }
    //~ opstr1.close(); // close the text file
    //~ opstr2.close();

    //ptr_matZ->printMatrix(); // sanity check
    //std::cout << " .......... Printing the final pixel_pt vector .......... " << std::endl;
    //pixel_pt.printVector3D(); // sanity check
}

void Camera::printCameraAttributes(){
    std::cout << '\n' << "Camera Attributes -------------------------------------------------->" << std::endl;
    //std::cout << focal_pt_x << " " << focal_pt_y << " " << focal_pt_z << " <----- focal point" << std::endl;
    //std::cout << eye_ptr << " <----- eye pointer address" << std::endl;
    //std::cout << lookat_pt_x << " " << lookat_pt_y << " " << lookat_pt_z << " <----- look at point" << std::endl;
    std::cout << eye.x << " " << eye.y << " " << eye.z << " <----- focal point" << std::endl;
    std::cout << lookat.x << " " << lookat.y << " " << lookat.z << " <----- look at point" << std::endl;
    std::cout << up.x << " " << up.y << " " << up.z << " <----- up vector" << std::endl;
    std::cout << focalLength << " <----- focal length" << std::endl;
    std::cout << left << " " << bottom << " " << right << " " << top << " <----- bounds" << std::endl;
    std::cout << res_u << " " << res_v << " <----- resolution" << '\n' << std::endl;

} // constructor for printing camera attributes

Camera::~Camera(){
}//Camera object destructor


#endif // CAMERA_H_INCLUDED
