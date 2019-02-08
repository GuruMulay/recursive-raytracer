#ifndef CAMERAFILEPARSER_H_INCLUDED
#define CAMERAFILEPARSER_H_INCLUDED

#include "camera.h"
//#include "testMakefile.h"
//#include "cameraFileParser.cpp"

void getFocalLength(std::string inputFileName, double &focalLength);
void getBounds(std::string inputFileName, double &left, double &bottom, double &right, double &top);
void getResolution(std::string inputFileName, int &res_u, int &res_v);
void getFocalPoint(std::string inputFileName, double &focal_pt_x, double &focal_pt_y, double &focal_pt_z);
void getLookAtPoint(std::string inputFileName, double &lookat_pt_x, double &lookat_pt_y, double &lookat_pt_z);
void getUpVector(std::string inputFileName, double &up_x, double &up_y, double &up_z);
Camera* getCameraAttributes(std::string inputFileName);


#endif // CAMERAFILEPARSER_H_INCLUDED
