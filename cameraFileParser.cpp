#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "cameraFileParser.h"
//#include "camera.h"

void getFocalLength(std::string inputFileName, double &focalLength){
    //std::cout << inputFileName << " <-----" << std::endl;

    std::string line;
    std::string dValue ("d");
    std::string term1; // to store d word
    //double d = 5; // to store d value

    std::ifstream cameraFile(inputFileName.c_str());
    std::size_t check = std::string::npos;
    //std::cout << check << " <----- check " << '\n';

    while (check == std::string::npos)
        {
            getline(cameraFile, line);
            check = line.find(dValue);
            //std::cout << check << '\n';
            //std::cout << line << '\n';
            //std::cin.get();
        } // while1

    //std::cout << line << " FOUND here " << '\n';
    std::stringstream ss(line);
    ss >> term1 >> focalLength;
    //std::cout << focalLength << " <----- focalLength" << std::endl;

    cameraFile.close();

    //return focalLength;
}

void getBounds(std::string inputFileName, double &left, double &bottom, double &right, double &top){
    std::string line;
    std::string bounds ("bounds");
    std::string term1; // to store bounds word

    std::ifstream cameraFile(inputFileName.c_str());
    std::size_t check = std::string::npos;

    while (check == std::string::npos)
        {
            getline(cameraFile, line);
            check = line.find(bounds);
        } // while1

    std::stringstream ss(line);
    ss >> term1 >> left >> bottom >> right >> top;

    cameraFile.close();
}

void getResolution(std::string inputFileName, int &res_u, int &res_v){
    std::string line;
    std::string res ("res");
    std::string term1; // to store res word

    std::ifstream cameraFile(inputFileName.c_str());
    std::size_t check = std::string::npos;

    while (check == std::string::npos)
        {
            getline(cameraFile, line);
            check = line.find(res);
        } // while1

    std::stringstream ss(line);
    ss >> term1 >> res_u >> res_v;

    cameraFile.close();
}

void getFocalPoint(std::string inputFileName, double &focal_pt_x, double &focal_pt_y, double &focal_pt_z){
    std::string line;
    std::string eye ("eye");
    std::string term1; // to store eye word

    std::ifstream cameraFile(inputFileName.c_str());
    std::size_t check = std::string::npos;

    while (check == std::string::npos)
        {
            getline(cameraFile, line);
            check = line.find(eye);
        } // while1

    std::stringstream ss(line);
    ss >> term1 >> focal_pt_x >> focal_pt_y >> focal_pt_z;

    cameraFile.close();
}

void getLookAtPoint(std::string inputFileName, double &lookat_pt_x, double &lookat_pt_y, double &lookat_pt_z){
    std::string line;
    std::string look ("look");
    std::string term1; // to store look word

    std::ifstream cameraFile(inputFileName.c_str());
    std::size_t check = std::string::npos;

    while (check == std::string::npos)
        {
            getline(cameraFile, line);
            check = line.find(look);
        } // while1

    std::stringstream ss(line);
    ss >> term1 >> lookat_pt_x >> lookat_pt_y >> lookat_pt_z;

    cameraFile.close();
}

void getUpVector(std::string inputFileName, double &up_x, double &up_y, double &up_z){
    std::string line;
    std::string up ("up");
    std::string term1; // to store up word

    std::ifstream cameraFile(inputFileName.c_str());
    std::size_t check = std::string::npos;

    while (check == std::string::npos)
        {
            getline(cameraFile, line);
            check = line.find(up);
        } // while1

    std::stringstream ss(line);
    ss >> term1 >> up_x >> up_y >> up_z;

    cameraFile.close();
}

Camera* getCameraAttributes(std::string inputFileName){

    Camera* ptrCamera;
    //ptrCamera = new Camera(); // default initialization
    Vector3D<double> eye;
    Vector3D<double> lookat;
    Vector3D<double> up;
    double focalLength;

    double left, bottom, right, top;
    int res_u, res_v;

    //double focal_pt_x, focal_pt_y, focal_pt_z;
    //Vector3D<double>* focal_ptr = &eye;
    //ptrCamera->eye_ptr = &eye;
    //double lookat_pt_x, lookat_pt_y, lookat_pt_z;
    //double up_x, up_y, up_z;

    std::ifstream cameraFile(inputFileName.c_str());

    if (cameraFile.is_open())
        {
            getFocalLength(inputFileName, focalLength);

            getBounds(inputFileName, left, bottom, right, top);
            //std::cout << left << " " << bottom << " " << right << " " << top << " <----- bounds" << std::endl;
            getResolution(inputFileName, res_u, res_v);
            //std::cout << res_u << " " << res_v << " <----- resolution" << std::endl;

            getFocalPoint(inputFileName, eye.x, eye.y, eye.z);
            //std::cout << eye.x << " eyex " << eye.y << " eyey " << eye.z  << " eyez " << std::endl;

            getLookAtPoint(inputFileName, lookat.x, lookat.y, lookat.z);
            //std::cout << lookat.x << " lookat_pt_x " << lookat.y << " lookat_pt_y " << lookat.z  << " lookat_pt_z " << std::endl;

            getUpVector(inputFileName, up.x, up.y, up.z);

            ptrCamera = new Camera(left, bottom, right, top, focalLength, res_u, res_v);

            ptrCamera->eye.x = eye.x;
            ptrCamera->eye.y = eye.y;
            ptrCamera->eye.z = eye.z;

            ptrCamera->lookat.x = lookat.x;
            ptrCamera->lookat.y = lookat.y;
            ptrCamera->lookat.z = lookat.z;

            ptrCamera->up.x = up.x;
            ptrCamera->up.y = up.y;
            ptrCamera->up.z = up.z;

//            std::cout << ptrCamera->eye.x << " eyex " << ptrCamera->eye.y << " eyey " << ptrCamera->eye.z  << " eyez" << std::endl;
//            std::cout << lookat.x << " lookat_pt_x " << lookat.y << " lookat_pt_y " << lookat.z  << " lookat_pt_z " << std::endl;
//            std::cout << up.x << " up x " << up.y << " up y " << up.z  << " up z " << std::endl;

            //std::cout << focal_ptr << " eye address " << std::endl;

        }//if

    else
        {
            std::cout << "Unable to open " << inputFileName << std::endl;
            return 0;
        }

    return ptrCamera;
}
