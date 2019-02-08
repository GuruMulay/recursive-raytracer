//============================================================================
// Name        : raytracer.cpp
// Created on  : 10 November 2016
// Author      : @guru
// Version     : 4.15
// Copyright   : Copyright
// Description : Program for Light Tracing
//============================================================================

#include <iostream>
#include <math.h>
#include <fstream>  // for parsing text files
#include <sstream>  // for parsing text files
#include <string>
#include <iomanip>  // std::setprecision(2)
#include <limits>   // to get infinity!

#include <vector>
#include <ctime> // for time calculation
//#include <cstdint>  // for pointer address

#include "vector3d.h"
#include "matrix2d.h"
#include "ray.h"
#include "camera.h"
#include "light.h"

#include "model.h"
#include "plyobject.h"
#include "sphere.h"
#include "point.h"

#include "cameraFileParser.h"
#include "cameraFileParser.cpp"
#include "sceneFileParser.h"
#include "sceneFileParser.cpp"

#include <map>
#include <sstream>
#include <cstdlib>

//using namespace std;

#define CN1 92
#define CN2 104

#define RN1 388
#define RN2 388

int getVertices(std::string inputFileName){
    //std::cout << inputFileName << "<-----" << std::endl;
    std::string line;
    std::string eleVert ("element vertex");
    std::string term1, term2;
    int nVertices = 0;

    std::ifstream plyFile(inputFileName.c_str());
    std::size_t check = std::string::npos;
    //std::cout << check << '\n';

    if (plyFile.is_open())
    {
        while (check == std::string::npos)
            {
            getline(plyFile, line);
            check = line.find(eleVert);
            //std::cout << check << '\n';
            //std::cout << line << '\n';
            //std::cin.get();
            } // while1

        //std::cout << line << " FOUND here ";
        std::stringstream ss(line);
        ss >> term1 >> term2 >> nVertices;
        //std::cout << nVertices << " <-----" << std::endl;

        plyFile.close();
    }//if

    else
        {
            std::cout << "Unable to open " << inputFileName << std::endl;
            return 0;
        }
    return nVertices;
}

int getFaces(std::string inputFileName){
    //std::cout << inputFileName << "<-----" << std::endl;
    std::string line;
    std::string eleFace ("element face");
    std::string term1, term2;
    int nFaces = 0;

    std::ifstream plyFile(inputFileName.c_str());
    std::size_t check = std::string::npos;
    //std::cout << check << '\n';

    if (plyFile.is_open())
    {
        while (check == std::string::npos)
            {
            getline(plyFile, line);
            check = line.find(eleFace);
            //std::cout << check << '\n';
            //std::cout << line << '\n';
            //std::cin.get();
            } // while1

        //std::cout << line << " FOUND here ";
        std::stringstream ss(line);
        ss >> term1 >> term2 >> nFaces;
        //std::cout << nVertices << " <-----" << std::endl;

        plyFile.close();
    }//if

    else
        {
            std::cout << "Unable to open " << inputFileName << std::endl;
            return 0;
        }
    return nFaces;
}

bool is_file_present(const std::string& inputFileName){
    std::ifstream inFile(inputFileName.c_str());
    return inFile.good();
}

bool is_ply(const std::string& inputFileName){
    std::string line;
    std::ifstream inFile(inputFileName.c_str());
    if (inFile.is_open())
        {
            getline(inFile, line);
            if(line == "ply")
                {
                    //std::cout << "It is a ply file!" << std::endl;
                    return true;
                }
            inFile.close();
        }//if
    else
        {
            std::cout << "Unable to open " << inputFileName << std::endl;
            return false;
        }
}

PlyObject* processPLY (std::string inputFileName, int nVert, int nFace){
    PlyObject* ptrPlyObject;
    ptrPlyObject = new PlyObject(4, nVert, nFace, 5); // initialize the PlyObject | 4 is for [x, y, z, 1] => homogeneous coordinates

    //std::cout << ptrPlyObject->n_vertices << " <------- nVerts " << ptrPlyObject->n_faces << " <------- nFaces" << std::endl;
    ptrPlyObject->getPlyObjectVerticesFaceAttributes(inputFileName);

    return ptrPlyObject;
} // ply processing assignment 3


//------------------------------------------------------------------------------------------------------------------------------------
//IMPORTANT NOTE: determinant M1, M2, and M3 can be calculated using the same expression used for determinant M with following changes:
//For M1: replace vector b with vector l (pixeL)
//For M2: replace vector c with vector l (pixeL)
//For M3: replace vector d with vector (a-l) i.e., (a-pixeL)
//------------------------------------------------------------------------------------------------------------------------------------
inline double calculateDeterminantM(const Vector3D<double> &a, const Vector3D<double> &b, const Vector3D<double> &c, const Vector3D<double> &d){
    double detM;
    // using simplified expression for the determinant
    detM = d.z*(c.y*(b.x-a.x) + a.y*(c.x-b.x) + b.y*(a.x-c.x)) +
           d.y*(c.z*(a.x-b.x) + a.z*(b.x-c.x) + b.z*(c.x-a.x)) +
           d.x*(c.z*(b.y-a.y) + a.z*(c.y-b.y) + b.z*(a.y-c.y));

    //det = c1.x*(c2.y*c3.z - c2.z*c3.y) - c2.x*(c1.y*c3.z - c1.z*c3.y) + c3.x*(c1.y*c2.z - c1.z*c2.y);
    //std::cout << detM << " <----- inside function " << std::endl;
    return detM;
}

inline void calculateDeterminantMAndM3(const Vector3D<double> &a, const Vector3D<double> &b, const Vector3D<double> &c, const Vector3D<double> &d, const Vector3D<double> &l, double &detM, double &detM3){
    double D1, D2, D3;
    double t1, t2, t3;
    t1 = (b.x-a.x);
    t2 = (c.x-b.x);
    t3 = (a.x-c.x);

    D1 =  c.y*t1 + a.y*t2 + b.y*t3;
    D2 = -c.z*t1 - a.z*t2 - b.z*t3;
    D3 = (c.z*(b.y-a.y) + a.z*(c.y-b.y) + b.z*(a.y-c.y));

//////////    D1 = (c.y*(b.x-a.x) + a.y*(c.x-b.x) + b.y*(a.x-c.x));
//////////    D2 = (c.z*(a.x-b.x) + a.z*(b.x-c.x) + b.z*(c.x-a.x));
//////////    D3 = (c.z*(b.y-a.y) + a.z*(c.y-b.y) + b.z*(a.y-c.y));

    detM =   d.z*D1 + d.y*D2 + d.x*D3;
    detM3 = (a.z-l.z)*D1 + (a.y-l.y)*D2 + (a.x-l.x)*D3;
}

inline void calculateDeterminantM1AndM2(const Vector3D<double> &a, const Vector3D<double> &b, const Vector3D<double> &c, const Vector3D<double> &d, const Vector3D<double> &l, double &detM1, double &detM2){
    double D1, D2, D3;
    double t1, t2, t3;
    t1 = (a.y-l.y);
    t2 = (a.z-l.z);
    t3 = (a.x-l.x);

    D1 = ( (d.z*(t1)) - (d.y*(t2)) );
    D2 = ( (d.z*(t3)) - (d.x*(t2)) );
    D3 = ( (d.y*(t3)) - (d.x*(t1)) );

//////////    D1 = ( (d.z*(a.y-l.y)) - (d.y*(a.z-l.z)) );
//////////    D2 = ( (d.z*(a.x-l.x)) - (d.x*(a.z-l.z)) );
//////////    D3 = ( (d.y*(a.x-l.x)) - (d.x*(a.y-l.y)) );

    detM2 =   (a.x-b.x)*D1 - (a.y-b.y)*D2 + (a.z-b.z)*D3;
    detM1 = -((a.x-c.x)*D1 - (a.y-c.y)*D2 + (a.z-c.z)*D3); // minus for column swap compensation
}


void initializeArrayElements(double arr[], int arrSize, double initValue){
    //std::cout << "----------- printing the array of t for a single pixel -----------" << std::endl;
    for(int i=0; i<arrSize; i++){
        arr[i] = initValue;
    }
}

void printArrayElements(double arr[], int arrSize){
    std::cout << "----------- printing the array of t for a single pixel -----------" << std::endl;
    for(int i=0; i<arrSize; i++){
        std::cout << arr[i] << ", ";
    }
    std::cout << '\n' << std::endl;
}

inline double getMinPositiveElement(const double arr[], const int arrSize){
    double tMin = INFINITY; // wrong-> needs to be a negative number!!!!! because for negative t, we write 239, 239, 239 in the final step
    for(int i=0; i<arrSize; i++){
        //if(tMin>arr[i] && arr[i]>0){tMin = arr[i];} // find minimum positive array element
        if(arr[i]>0 && arr[i]<tMin){tMin = arr[i];}
    }
    return tMin;
}


void calculateDepthValues(PlyObject* ptr_ply1, Camera* ptr_cm){
    //std::cout << ptr_ply1 << " <--- ptr_ply1 " << ptr_cm << " <--- ptr_cm" << std::endl; // sanity check
/*  checks
    ptr_ply1->ptr_faceMat->printMatrix();
    ptr_ply1->ptr_verticesMat->printMatrix();
    ptr_cm->ptr_pixMatY->printMatrix();
    ptr_cm->ptr_rayMatY->printMatrix();
    ptr_ply1->ptr_faceMat->setElement(1,2,1005); // changes reflected in main
*/
    //uintptr_t addressVertMat = ptr_ply1->ptr_verticesMat;   //, addressFaceMat;
    ptr_cm->ptr_depthMat = new Matrix2D<double>(ptr_cm->res_v, ptr_cm->res_u, -1); //initialize to -1 so that we could distinguish from required positive values
    //ptr_cm->ptr_depthMat->printMatrix();

    Vector3D<double> a, b, c, d, pixeL; // store temporary a, b, c, d
    unsigned v1, v2, v3;

    double tOfAPixel[ptr_ply1->n_faces]; //stores the t values for a single pixel after scanning all the faces; finally we pick minimum t>0
    initializeArrayElements(tOfAPixel, ptr_ply1->n_faces, -1.1-10); // set it to negative and not zero since object could be touching the image plane
    double tMin;

    int row = 0;
    int column = 0;
    int face = 0;
    double detM=0, detM1=0, detM2=0, detM3=0;
    //IMPORATANT CONSTANT EPSILON = 10^-10;
    double eps = 0.0000000001;
    std::cout << eps << " <<<----------------------------- epsilon value for [zero] comparison " << std::endl;

    double beta, gamma, tDepth;

    double cl, v;

    /*
    //-------------------------------- this section sets the values for a, b, c, d vectors -----------------------------//
    d.x = ptr_cm->ptr_rayMatX->getElement(0,0);
    d.y = ptr_cm->ptr_rayMatY->getElement(0,0);
    d.z = ptr_cm->ptr_rayMatZ->getElement(0,0);
    std::cout << " ----------------- printing d vector after setting ----------------------------- " << std::endl;
    d.printVector3D();

    pixeL.x = ptr_cm->ptr_pixMatX->getElement(0,0);
    pixeL.y = ptr_cm->ptr_pixMatY->getElement(0,0);
    pixeL.z = ptr_cm->ptr_pixMatZ->getElement(0,0);
    std::cout << " ----------------- printing pixeL vector after setting ----------------------------- " << std::endl;
    pixeL.printVector3D();

    v1 = ptr_ply1->ptr_faceMat->getElement(0,1); // 0th row 1st column
    v2 = ptr_ply1->ptr_faceMat->getElement(0,2);
    v3 = ptr_ply1->ptr_faceMat->getElement(0,3);
    std::cout << v1 << " " << v2 << " " << v3 << " <---------------------------------- printing vi values " << std::endl;

    a.x = ptr_ply1->ptr_verticesMat->getElement(v1,0);
    a.y = ptr_ply1->ptr_verticesMat->getElement(v1,1);
    a.z = ptr_ply1->ptr_verticesMat->getElement(v1,2);
    std::cout << " ----------------- printing a vector after setting --------------------------- " << std::endl;
    a.printVector3D();

    b.x = ptr_ply1->ptr_verticesMat->getElement(v2,0);
    b.y = ptr_ply1->ptr_verticesMat->getElement(v2,1);
    b.z = ptr_ply1->ptr_verticesMat->getElement(v2,2);
    std::cout << " ----------------- printing b vector after setting --------------------------- " << std::endl;
    b.printVector3D();

    c.x = ptr_ply1->ptr_verticesMat->getElement(v3,0);
    c.y = ptr_ply1->ptr_verticesMat->getElement(v3,1);
    c.z = ptr_ply1->ptr_verticesMat->getElement(v3,2);
    std::cout << " ----------------- printing c vector after setting --------------------------- " << std::endl;
    c.printVector3D();

    //std::cout << ptr_ply1->ptr_faceMat->getElement(1,0)<< " <----------------- face element in function" << std::endl;

    std::cout << " <----------------- printing determinant values -------------------------- " << std::endl;
    detM = calculateDeterminantM(a, b, c, d);
    std::cout << detM << " <------------ detM " << std::endl;

    detM1 = calculateDeterminantM(a, pixeL, c, d);                //NOTE: b is replaced by pixeL
    std::cout << detM1/detM << " <------------ detM1/detM " << std::endl;

    detM2 = calculateDeterminantM(a, b, pixeL, d);                //NOTE: c is replaced by pixeL
    std::cout << detM2/detM << " <------------ detM2/detM " << std::endl;

    detM3 = calculateDeterminantM(a, b, c, (a-pixeL));            //NOTE: d is replaced by a-pixeL
    std::cout << detM3/detM << " <------------ detM3/detM " << std::endl;
    //-----------------------------------------------------------------------------------------------------------------//

    printArrayElements(tOfAPixel, ptr_ply1->n_faces); //stores the t values for a single pixel after scanning all the faces; finally we pick minimum t>0
    tMin = getMinPositiveElement(tOfAPixel, ptr_ply1->n_faces);
    std::cout << tMin << " <------------ tMin for that pixel " << std::endl;

    //-----------------------------------------------------------------------------------------------------------------//
    */

    int sp = 0;
    double dt;
    double total_time = 0;
    int start_sec_matrix, stop_sec_matrix;
    start_sec_matrix = clock();


//////////    std::ofstream op1;
//////////    op1.open("tvals.txt", std::ofstream::out);

    // for every pixel and ray through that pixel::::: (starting from left-bottom corner)
    for(row=0; row < ptr_cm->res_v; row++){ // ptr_cm->res_u
        for(column=0; column < ptr_cm->res_u; column++){ // ptr_cm->res_v

            // ---------------------------- Firstly, get values of d and pixeL vectors for Pixel Number (row, column) ----- //

            d.x = ptr_cm->ptr_rayMatX->getElement(row,column);
            d.y = ptr_cm->ptr_rayMatY->getElement(row,column);
            d.z = ptr_cm->ptr_rayMatZ->getElement(row,column);
////            std::cout << " ----------------- printing d vector after setting ----------------------------- " << std::endl;
////            d.printVector3D();

            pixeL.x = ptr_cm->ptr_pixMatX->getElement(row,column);
            pixeL.y = ptr_cm->ptr_pixMatY->getElement(row,column);
            pixeL.z = ptr_cm->ptr_pixMatZ->getElement(row,column);
////            std::cout << " ----------------- printing pixeL vector after setting ----------------------------- " << std::endl;
////            pixeL.printVector3D();

/*------------ Sphere Exit (doen't work for ellelltri.ply) ----------------------------------------------------------------
            cl = (ptr_ply1->object_center - pixeL).getMagnitude();
            v = (ptr_ply1->object_center - pixeL).dotProduct(d);
            //std::cout << " cl and v and r " << cl << ", " << v << ", " << (ptr_ply1->radiusOfSphere) << std::endl;
            //int ch = getchar();

            if((ptr_ply1->radiusOfSphere)*(ptr_ply1->radiusOfSphere) < (cl*cl-v*v)){
                tOfAPixel[face] = -2.55;
                sp++;
                continue;
            }
*/ //----------------------------------------------------------------------------------------------------------------------
            // ------------------------------------------------------------------------------------------------------ //
            // For every face::::: Find beta, gamma, and tDepth values
            for(face=0; face<ptr_ply1->n_faces; face++){

////                std::cout << " >>>>>>>>>>>>>>>>>>> for face number " << face << " <<<<<<<<<<<<<<<<< " << std::endl;
                ptr_ply1->ptr_faceMat->getElementXYZ(face, v1, v2, v3); // 0th row 1st column
//                v2 = ptr_ply1->ptr_faceMat->getElement(face, 2);
//                v3 = ptr_ply1->ptr_faceMat->getElement(face, 3);
////                std::cout << v1 << " " << v2 << " " << v3 << " <---------------------------------- printing vi values " << std::endl;

                ptr_ply1->ptr_verticesMat->getElementXYZTriangle(v1, v2, v3, a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z);
                //////////ptr_ply1->ptr_faceMat->getElementXYZTriangle(*(ptr_ply1->ptr_verticeMat), face, a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z);
//                ptr_ply1->ptr_verticesMat->getElement1(v1,1, a.y);
//                ptr_ply1->ptr_verticesMat->getElement1(v1,2, a.z);
////                std::cout << " ----------------- printing a vector after setting --------------------------- " << std::endl;
////                a.printVector3D();

//////////                ptr_ply1->ptr_verticesMat->getElementXYZ(v2,0, b.x, b.y, b.z);
//                ptr_ply1->ptr_verticesMat->getElement1(v2,1, b.y);
//                ptr_ply1->ptr_verticesMat->getElement1(v2,2, b.z);
////                std::cout << " ----------------- printing b vector after setting --------------------------- " << std::endl;
////                b.printVector3D();

//////////                ptr_ply1->ptr_verticesMat->getElementXYZ(v3,0, c.x, c.y, c.z);
//                ptr_ply1->ptr_verticesMat->getElement1(v3,1, c.y);
//                ptr_ply1->ptr_verticesMat->getElement1(v3,2, c.z);
////                std::cout << " ----------------- printing c vector after setting --------------------------- " << std::endl;
////                c.printVector3D();

////                std::cout << " <----------------- printing determinant values -------------------------- " << std::endl;
                calculateDeterminantMAndM3(a, b, c, d, (pixeL), detM, detM3);
//////////                op1 << detM << '\n';
////                std::cout << detM << " <------------ detM " << std::endl;
                // if determinant is zero then the ray does not intersect the face => there is no solution; hence set tOfAPixel[face] to -2 and DONT break! (but no need to calculate beta, gamma, and tDepth)
                if(fabs(detM) < eps){
                    ///////tOfAPixel[face] = -2; // removed to save time; tDepth is already initialized to a negative value
////                    std::cout << " set t value to -2 (detM = 0) and continue to next iteration i.e., the next face! " << std::endl;
                    continue;
                } // if detM == 0

                //detM3 = calculateDeterminantM(a, b, c, (a-pixeL));            //NOTE: d is replaced by a-pixeL
                tDepth = detM3/detM;
////                std::cout << detM3/detM << " <------------ tDepth = detM3/detM " << std::endl;

//////////                if(tDepth<=0){op1 << tDepth << ", " << row << ", " << column << ", " << face << ", " << a.x << ", " << b.x << ", " << c.x << '\n';
//////////                    op1 << v1 << ", " << v2 << ", " << v3 << '\n';
//////////                    op1 << d.x << ", " << d.y << ", " << d.z << ", " << pixeL.x << ", " << pixeL.y << ", " << pixeL.z << ", " << '\n';
//////////                }

                if(tDepth<=0){
                    ///////tOfAPixel[face] = tDepth; // removed to save time; tDepth is already initialized to a negative value
////                    std::cout << " set t value to tDepth and continue to next iteration i.e., the next face! " << std::endl;
                    continue;
                } // if tDepth is negative

                calculateDeterminantM1AndM2(a, b, c, d, (pixeL), detM1, detM2);
                //detM1 = calculateDeterminantM(a, pixeL, c, d);                //NOTE: b is replaced by pixeL
                beta = detM1/detM;
////                std::cout << detM1/detM << " <------------ beta = detM1/detM " << std::endl;
                if(beta<0){
                    tOfAPixel[face] = -3; // Here, tDepth is positive already. But beta is negative => No intersection! So we gotta neglect this positive tDepth by setting it to -3.
////                    std::cout << " set t value to -3 (beta negative) and continue to  the next face! " << std::endl;
                    continue;
                } // if beta is negative

                //detM2 = calculateDeterminantM(a, b, pixeL, d);                //NOTE: c is replaced by pixeL
                gamma = detM2/detM;
////                std::cout << detM2/detM << " <------------ gamma = detM2/detM " << std::endl;
                if(gamma<0){
                    tOfAPixel[face] = -5; // Here, tDepth is positive already. But gamma is negative => No intersection! So we gotta neglect this positive tDepth by setting it to -5.
////                    std::cout << " set t value to -5 (gamma negative) and continue to  the next face! " << std::endl;
                    continue;
                } // if gamma is negative

                //Final step of sieving (EARLY EXIT/RETURN strategy), check if (beta+gamma) <= 1
                //////////////////////// CHANGE this for square shaped faces!!!!!!!!!!
                if((beta+gamma)>1){
                    tOfAPixel[face] = -7; // Here, tDepth is positive already. But (beta+gamma) > 1) => No intersection! So we gotta neglect this positive tDepth by setting it to -7.
////                    std::cout << " set t value to -7 (gamma negative) and continue to  the next face! " << std::endl;
                    continue;
                } // if (beta+gamma) > 1

                // When everything from pdf notes (ray-triangle intersection) is satisfied then assign tDepth which is positive to tOfAPixel[face] array element
                if(tDepth>=0){
                    tOfAPixel[face] = tDepth;
                    //std::cout << tDepth << " this is final step tDepth value; it has to be positive! " << std::endl;
                }

            } // for loop face

////            printArrayElements(tOfAPixel, ptr_ply1->n_faces); // stores the t values for a single pixel after scanning all the faces; finally we pick minimum t>0 minimum [for the closest first face that is hit by the ray]
            tMin = getMinPositiveElement(tOfAPixel, ptr_ply1->n_faces);
////            std::cout << tMin << " <------------ tMin for pixel number -> " << row << ", " << column << std::endl;

            ptr_cm->ptr_depthMat->setElement(row, column, tMin);
        } // for loop column

    } // for loop row

    std::cout << "sphere exits -----> " << sp << std::endl;
    stop_sec_matrix = clock();
    dt  = (stop_sec_matrix-start_sec_matrix)/double(CLOCKS_PER_SEC)*1;
    std::cout << "time for depth calulation: " << dt << std::endl;
    //std::cout << "time for matrix assignment: " << total_time << std::endl;
    //std::cout << ptr_ply1->ptr_faceMat << " <----- pointer address"  << std::endl;
    //check depthMat after scanning all the pixels
    //std::cout << "+++++++++++++++++ Checking after traversing all faces for every pixel +++++++++++++++" << std::endl;
    //ptr_cm->ptr_depthMat->printMatrix();
//////////    op1.close();
} // Depth Matrix Generation


void generatePpmImage(Camera* ptr_cm, std::string outputFilePpm){
    std::cout << "<<<<<<<<<<<<<<<<<<< Inside generatePpmImage function >>>>>>>>>>>>>>>>>>>>>>> " << std::endl;
    int row=0, column=0;
    double tMax = -1;
    double tMin = INFINITY;

    int redVal=0, greenVal=0, blueVal=0;
    double tRatio=0;

    // for loops to replace INFINITY with a negative number!!!!!
    for(row=0; row<ptr_cm->res_v; row++){
        for(column=0; column<ptr_cm->res_u; column++){
            if(ptr_cm->ptr_depthMat->getElement(row,column) == INFINITY){ptr_cm->ptr_depthMat->setElement(row,column,-5.55);}
        }
    }

    //-------------------------------------------------------------------------------------------------------------------------
    // find tmax and tmin
    for(row=0; row<ptr_cm->res_v; row++){
        for(column=0; column<ptr_cm->res_u; column++){
            if(tMax < ptr_cm->ptr_depthMat->getElement(row,column)){tMax = ptr_cm->ptr_depthMat->getElement(row,column);}
        }
    }

    for(row=0; row<ptr_cm->res_v; row++){
        for(column=0; column<ptr_cm->res_u; column++){
            if(ptr_cm->ptr_depthMat->getElement(row,column) < tMin && ptr_cm->ptr_depthMat->getElement(row,column) >= 0){tMin = ptr_cm->ptr_depthMat->getElement(row,column);}
        } // ">= 0 " is IMPORTANT tMin could possibly be 0
    }

    std::cout << "tMin = " << tMin << ", tMax = " << tMax << std::endl;

    //-------------------------------------------------------------------------------------------------------------------------
    // FINAL FINAL STEP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    std::ofstream op_ppm;
    op_ppm.open(outputFilePpm.c_str(), std::ofstream::out);

    //PPM HEADER -----------------------------------------------------
    op_ppm << "P3" << '\n';
    op_ppm << ptr_cm->res_u << " " << ptr_cm->res_v << " " << "255" << '\n';

    //PPM PIXEL RGB VALUES -------------------------------------------
//// inverts about center: both axes -u, -v
////    for(row=ptr_cm->res_u-1; row>=0; row--){
////        for(column=ptr_cm->res_v-1; column>=0; column--){

//// column iteration reversal => -v => flip vertically ///////////////
//// row iteration reversal => -u => flip horizontally ////////////////
//////////    for(column=ptr_cm->res_v-1; column>=0; column--){
//////////        for(row=0; row<ptr_cm->res_u; row++){


    //********************* rows scanned in inverted manner before printing because they were generated from bottom->top
    //********************* columns scanned in left to right manner
    for(row=ptr_cm->res_v-1; row>=0; row--){
        for(column=0; column<ptr_cm->res_u; column++){

            if(ptr_cm->ptr_depthMat->getElement(row,column) < 0){
                op_ppm << "239" << " " << "239" << " " << "239" << " ";
            }

            else{
                tRatio = 2 * (ptr_cm->ptr_depthMat->getElement(row,column) - tMin) / (tMax - tMin);
                redVal = std::max((double)0, (255 * (1 - tRatio))); // typecasted 0 from int to double to match both the arguments for std::max() function
                blueVal = std::max((double)0, (255 * (tRatio - 1)));
                greenVal = 255 - blueVal - redVal;

                op_ppm << redVal << " " << greenVal << " " << blueVal << " ";
            }

        } //for column

        op_ppm << '\n'; // newline after each row in PPM file

    } // for row

    op_ppm << '\n'; // to match the given PPM files
    op_ppm.close();

//////////    ptr_cm->ptr_depthMat->printMatrix(); // sanity check

} // Generate PPM Output Image


void generatePpmImageSphere(Camera* ptr_cm, std::string outputFilePpm){
    // duplicate of generatePpmImage() ???
    std::cout << "<<<<<<<<<<<<<<<<<<< Inside generatePpmImage for Sphere function >>>>>>>>>>>>>>>>>>>>>>> " << std::endl;
    int row=0, column=0;
    double tMax = -1;
    double tMin = INFINITY;

    int redVal=0, greenVal=0, blueVal=0;
    double tRatio=0;

    // for loops to replace INFINITY with a negative number!!!!!
    for(row=0; row<ptr_cm->res_v; row++){
        for(column=0; column<ptr_cm->res_u; column++){
            if(ptr_cm->ptr_depthMat->getElement(row,column) == INFINITY){ptr_cm->ptr_depthMat->setElement(row,column,-5.55);}
        }
    }

    //-------------------------------------------------------------------------------------------------------------------------
    // find tmax and tmin
    for(row=0; row<ptr_cm->res_v; row++){
        for(column=0; column<ptr_cm->res_u; column++){
            if(tMax < ptr_cm->ptr_depthMat->getElement(row,column)){tMax = ptr_cm->ptr_depthMat->getElement(row,column);}
        }
    }

    for(row=0; row<ptr_cm->res_v; row++){
        for(column=0; column<ptr_cm->res_u; column++){
            if(ptr_cm->ptr_depthMat->getElement(row,column) < tMin && ptr_cm->ptr_depthMat->getElement(row,column) >= 0){tMin = ptr_cm->ptr_depthMat->getElement(row,column);}
        } // ">= 0 " is IMPORTANT tMin could possibly be 0
    }

    std::cout << "tMin = " << tMin << ", tMax = " << tMax << std::endl;

    //-------------------------------------------------------------------------------------------------------------------------
    // FINAL FINAL STEP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    std::ofstream op_ppm;
    op_ppm.open(outputFilePpm.c_str(), std::ofstream::out);

    //PPM HEADER -----------------------------------------------------
    op_ppm << "P3" << '\n';
    op_ppm << ptr_cm->res_u << " " << ptr_cm->res_v << " " << "255" << '\n';

    //PPM PIXEL RGB VALUES -------------------------------------------
//// inverts about center: both axes -u, -v
////    for(row=ptr_cm->res_u-1; row>=0; row--){
////        for(column=ptr_cm->res_v-1; column>=0; column--){

//// column iteration reversal => -v => flip vertically ///////////////
//// row iteration reversal => -u => flip horizontally ////////////////
//////////    for(column=ptr_cm->res_v-1; column>=0; column--){
//////////        for(row=0; row<ptr_cm->res_u; row++){


    //********************* rows scanned in inverted manner before printing because they were generated from bottom->top
    //********************* columns scanned in left to right manner
    for(row=ptr_cm->res_v-1; row>=0; row--){
        for(column=0; column<ptr_cm->res_u; column++){

            if(ptr_cm->ptr_depthMat->getElement(row,column) < 0){
                op_ppm << "239" << " " << "239" << " " << "239" << " ";
            }

            else{
                tRatio = 2 * (ptr_cm->ptr_depthMat->getElement(row,column) - tMin) / (tMax - tMin);
                redVal = std::max((double)0, (255 * (1 - tRatio))); // typecasted 0 from int to double to match both the arguments for std::max() function
                blueVal = std::max((double)0, (255 * (tRatio - 1)));
                greenVal = 255 - blueVal - redVal;

                op_ppm << redVal << " " << greenVal << " " << blueVal << " ";
            }

        } //for column

        op_ppm << '\n'; // newline after each row in PPM file

    } // for row

    op_ppm << '\n'; // to match the given PPM files
    op_ppm.close();
} // Generate PPM Output Image


Sphere* rayToWhichSphere(Ray* ray, Sphere** ptr_sphArray, unsigned& nSpheres, int& intersectSphere){
    //std::cout << ray << " ptrRay inside function" <<std::endl;

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ray->ptr_bestSph = NULL; // SAVED, FIANLLY!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ray->bestT_sphere = INFINITY; // SAVED MY DAY, FIANLLY! // SAVED ANOTHER DAY, AGAIN by putting this outside the for loop for spheres
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for (unsigned n=0; n<nSpheres; n++){
        //ray->d.printVector3D();
        // intersect values: 2 => default value 1 => return 1 from ray->raySphereIntersectionTest meaning intersected
        // 0 => return 0 from ray->raySphereIntersectionTest meaning no intersection
        intersectSphere = ray->raySphereIntersectionTest(ptr_sphArray[n]);
        //std::cout << intersect << " ? " << n << " <---------- intersection ? sphere number" << std::endl;
    }

    return ray->ptr_bestSph;
}


Model* rayToWhichModel(Ray* ray, Model** ptr_modelArray, unsigned& nModels, int& intersectModel){
    //std::cout << ray << " ptrRay inside function" <<std::endl;

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ray->ptr_bestModel = NULL; // SAVED, FIANLLY!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ray->bestT_model = INFINITY; // SAVED MY DAY, FIANLLY! // SAVED ANOTHER DAY, AGAIN by putting this outside the for loop for spheres
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for (unsigned n=0; n<nModels; n++){
        //ray->d.printVector3D();
        // intersect values: 2 => default value 1 => return 1 from ray->raySphereIntersectionTest meaning intersected
        // 0 => return 0 from ray->raySphereIntersectionTest meaning no intersection
        intersectModel = ray->rayModelIntersectionTest(ptr_modelArray[n], false);
        // std::cout << intersectModel << " ? " << n << " <---------- intersection (1==True; 0==False) ? model number" << std::endl;
    }

    return ray->ptr_bestModel;
}


void calculateRaySphereRGB(Camera* ptr_cm, Sphere* ptr_sph, unsigned& nSpheres, Light** ptr_lightArray, unsigned& nLightSources, Vector3D<double> &intersection, Vector3D<double> &rgbColor){
    Vector3D<double> surfaceNormal, ambient, vectL, vectR, vectV;
    ambient.x = 0;
    ambient.y = 0;
    ambient.z = 0;
    //intersection .printVector3D();

    surfaceNormal = intersection - ptr_sph->center;
    surfaceNormal.normalize();
    //surfaceNormal.printVector3D();
    //ptr_sph->rgb.printVector3D();

    // *************** Amibient reflection *********************
    rgbColor = ambient.pairwiseProduct(ptr_sph->ka);
    //rgbColor.printVector3D();

    for (unsigned n=0; n<nLightSources; n++){
        vectL = ptr_lightArray[n]->position - intersection;
        vectL.normalize(); // toL
        //vectL.printVector3D();

        if (surfaceNormal.dotProduct(vectL) > 0){
            // **************************** Diffuse reflection **************************************
            rgbColor = rgbColor + ( ptr_sph->kd.pairwiseProduct(ptr_lightArray[n]->rgb) ) * surfaceNormal.dotProduct(vectL);

            vectV = ptr_cm->eye - intersection;
            vectV.normalize(); // toC

            vectR = (surfaceNormal*2*surfaceNormal.dotProduct(vectL)) - vectL;
            // **************************** Specular reflection (Phong) *****************************
            rgbColor = rgbColor + ( ptr_sph->ka.pairwiseProduct(ptr_lightArray[n]->rgb) ) * pow(vectV.dotProduct(vectR), 16);
        }
    }

}

bool shadow(Vector3D<double>& pointOfIntersection, Vector3D<double>& lightLocation, unsigned& nSpheres, Sphere** ptr_sphArray, unsigned& nModels, Model** ptr_modelArray, bool debugFlag){
    
    Vector3D<double> pToL, pToLNormalized;
    double projLD;
    pToL = lightLocation - pointOfIntersection;
    pToLNormalized = lightLocation - pointOfIntersection;
    pToLNormalized.normalize();
    Ray pToLRay(pointOfIntersection.x, pointOfIntersection.y, pointOfIntersection.z, pToLNormalized.x, pToLNormalized.y, pToLNormalized.z);
    Ray* ptr_pToLRay = &pToLRay;

    projLD = fabs(pToL.dotProduct(ptr_pToLRay->d));  // gives distance between p and l
    //projLD = pToL.getMagnitude();  // gives distance between p and l

    // cs410 A4 2017 SHADOWS ----------------------------------------
    // check if toL goes through any other objects. If there is an object between light source and point of intersection
    // then it means there should be a shadow. Do not add light source's contribution in this case.

    bool ifShadow = false;
    bool sphereModelInbetween = false;
    int toLIntersectSphere = 2;
    int toLIntersectModel = 2;

	ptr_pToLRay->eye.x = pointOfIntersection.x;
	ptr_pToLRay->eye.y = pointOfIntersection.y;
	ptr_pToLRay->eye.z = pointOfIntersection.z;
	
	// -------------------------------------------------------
    ptr_pToLRay->d.x = pToLNormalized.x;
    ptr_pToLRay->d.y = pToLNormalized.y;
    ptr_pToLRay->d.z = pToLNormalized.z;
	
	
    if (nSpheres != 0) {
        for (unsigned n=0; n<nSpheres; n++){
            toLIntersectSphere = ptr_pToLRay->raySphereIntersectionTest(ptr_sphArray[n]);
            if (toLIntersectSphere == 1){
                ifShadow = true;
                break;
            }
        }
    }

	ptr_pToLRay->eye.x = pointOfIntersection.x;
	ptr_pToLRay->eye.y = pointOfIntersection.y;
	ptr_pToLRay->eye.z = pointOfIntersection.z;
	
	// -------------------------------------------------------
    ptr_pToLRay->d.x = pToLNormalized.x;
    ptr_pToLRay->d.y = pToLNormalized.y;
    ptr_pToLRay->d.z = pToLNormalized.z;
	
    if (nModels != 0) {
        for (unsigned n=0; n<nModels; n++){
            toLIntersectModel = ptr_pToLRay->rayModelIntersectionTest(ptr_modelArray[n], debugFlag);
			if (debugFlag) {
				std::cout << "toLIntersectModel=" << toLIntersectModel << std::endl;
			}
            if (toLIntersectModel == 1){
                ifShadow = true;
                break;
            }
        }
    }

    if ((ptr_pToLRay->bestT_sphere > 0 && ptr_pToLRay->bestT_sphere < projLD) || (ptr_pToLRay->bestT_model > 0 && ptr_pToLRay->bestT_model < projLD)){
        sphereModelInbetween = true;
    }

    if (ifShadow && sphereModelInbetween){
        return true;
    }

	if (debugFlag) {
		//std::cout << 512-row-1 << " " << column <<  " row, column ================================" << std::endl;
		std::cout << "shadow flags " << " ifShadow=" << ifShadow  << " sphereModelInbetween=" << sphereModelInbetween  << " toLIntersectSphere=" << toLIntersectSphere  << " toLIntersectModel=" << toLIntersectModel << std::endl;
		std::cout << "projLD=" << projLD << std::endl;
	}
	
	//delete ptr_pToLRay;
	//delete pToLRay;
    return false;

}


void rayTracer(Ray* ptrRay, int& row, int& column, Sphere** ptr_sphArray, unsigned& nSpheres, int& intersectSphere, Model** ptr_modelArray, unsigned& nModels, int& intersectModel, Light** ptr_lightArray, unsigned& nLightSources, Vector3D<double>& ambient, Vector3D<double>& accumRGB, Vector3D<double>& attenuationMultiplier, unsigned reflectionLevel){

	//~ std::cout << " ***** inside rayTracer ***** Printing initial attenuationMultiplier vector " << std::endl;
	//~ attenuationMultiplier.printVector3D();
	//~ std::cout << reflectionLevel << " refletionLevel" << std::endl;
    bool isSurfaceRgbDone;

	unsigned level;
	unsigned objectType; // 1 => sphere; 2 => model
	level = reflectionLevel;
	Vector3D<double> dInv; // opposite vector of d (direction of the ray)
	Vector3D<double> newD; // new value of vector d
	Vector3D<double> attenuationMultiplierNew;

	Vector3D<double> rgbColor;
	Vector3D<double> surfaceNormal, vectL, vectR, vectV;
  Vector3D<double> toEye;  // same as vectV

	Sphere* ptrSphTest = 0;  // Initialize to zero. Very Important! (in case sphere is not present)
	Model* ptrModelTest = 0; // Initialize to zero. Very Important! (in case model is not present)
	//~ int intersectSphere = 2;
	
	bool debugFlag = false;
	
    //~ accumRGB.printVector3D();
    //~ std::cout << nSpheres << " nSpheres" << std::endl;
    if (nSpheres != 0) {
		ptrSphTest = rayToWhichSphere(ptrRay, ptr_sphArray, nSpheres, intersectSphere);
	}
	//~ else {
		//~ accumRGB.x = 0;
		//~ accumRGB.y = 0;
		//~ accumRGB.z = 0;
	//~ }

	if (nModels != 0) {
		ptrModelTest = rayToWhichModel(ptrRay, ptr_modelArray, nModels, intersectModel);
	}
	//~ else {
		//~ accumRGB.x = 0;
		//~ accumRGB.y = 0;
		//~ accumRGB.z = 0;
	//~ }

    objectType = ptrRay->findBestIntersectionObject();
    //~ std::cout << ptrSphTest << " ptrSphTest"  << " | " << ptrRay->objectType << " ptrRay->objectType" << std::endl;
	
     /*if (row >= 90 && row < 93 && column >= 139 &&  column < 142) {
		 std::cout << row << " " << column <<  " row, column" << std::endl;
		 //std::cout << ptrSphTest << " ptrSphTest" << std::endl;
		 std::cout << ptrModelTest << " ptrModelTest" << std::endl;
		 surfaceNormal.printVector3D("surface normal");
	 }*/
	
	if (row >= 512-RN2 && row <= 512-RN1 && column >= CN1 && column <= CN2) {
		debugFlag = true;
	}
	
	if (ptrSphTest && ptrRay->objectType == 1){  // ptrSphTest will be 0 if there are no spheres; as it's initialized to 0.

		surfaceNormal = ptrRay->best_intersection_sphere - ptrRay->ptr_bestSph->center;
		surfaceNormal.normalize();
		// surfaceNormal.printVector3D("Sphere: surface normal = ");
		//ptrRay->ptr_bestSph->center.printVector3D();
		//std::cout << ptrRay << " ptrRay after best find" <<std::endl;
		//std::cin.get();
		//surfaceNormal.printVector3D();
		//ptr_sph->rgb.printVector3D();

		// *************** Amibient reflection *********************
		//ptrRay->ptr_bestSph->ka.printVector3D();
		rgbColor = ambient.pairwiseProduct(ptrRay->ptr_bestSph->ka);
		//std::cout << "Printing initial RGB vector " << std::endl;
		//rgbColor.printVector3D();

		for (unsigned n=0; n<nLightSources; n++){
			vectL = ptr_lightArray[n]->position - ptrRay->best_intersection_sphere;
			vectL.normalize(); // toL
			///std::cout << "Printing best intersection vector ---> " << n << std::endl;
			///ptrRay->best_intersection.printVector3D();
			//ptr_lightArray[n]->position.printVector3D();
			//vectL.printVector3D();

            bool hasShadow = false;
            // // cs410 A4 2017 SHADOWS ----------------------------------------
            // // check if toL goes through any other objects. If there is an object between light source and point of intersection
            // // then it means there should be a shadow. Do not add light source's contribution in this case.

            // // create a ray of toL vector
            // bool ifShadow = false;
            // int toLIntersectSphere = 2;
            // int toLIntersectModel = 2;
            // Ray toLRay(ptrRay->best_intersection_sphere.x, ptrRay->best_intersection_sphere.y, ptrRay->best_intersection_sphere.z, vectL.x, vectL.y, vectL.z);
            // Ray* ptr_toLRay = &toLRay;
            // // check for intersection with other objects
            // if (nSpheres != 0) {
            //     for (unsigned n=0; n<nSpheres; n++){
            //         // intersect values: 2 => default value 1 => return 1 from ray->raySphereIntersectionTest meaning intersected
            //         // 0 => return 0 from ray->raySphereIntersectionTest meaning no intersection
            //         toLIntersectSphere = ptr_toLRay->raySphereIntersectionTest(ptr_sphArray[n]);
            //         if (toLIntersectSphere == 1){
            //             ifShadow = true;
            //             break;
            //         }
            //     }
            // }
            // if (nModels != 0) {
            //     for (unsigned n=0; n<nModels; n++){
            //         toLIntersectModel = ptr_toLRay->rayModelIntersectionTest(ptr_modelArray[n]);
            //         if (toLIntersectModel == 1){
            //             ifShadow = true;
            //             break;
            //         }
            //     }
            // }

            hasShadow = shadow(ptrRay->best_intersection_sphere, ptr_lightArray[n]->position, nSpheres, ptr_sphArray, nModels, ptr_modelArray, debugFlag);


			if (surfaceNormal.dotProduct(vectL) > 0.000001 && hasShadow == false){
				// **************************** Diffuse reflection **************************************
				rgbColor = rgbColor + ( ptrRay->ptr_bestSph->kd.pairwiseProduct(ptr_lightArray[n]->rgb) ) * surfaceNormal.dotProduct(vectL);

				vectV = ptrRay->eye - ptrRay->best_intersection_sphere;
				//vectV = ptr_cm->eye - ray->best_intersection;
				vectV.normalize(); // toC

				vectR = (surfaceNormal*2*surfaceNormal.dotProduct(vectL)) - vectL; // spR
				// **************************** Specular reflection (Phong) *****************************
				rgbColor = rgbColor + ( ptrRay->ptr_bestSph->ks.pairwiseProduct(ptr_lightArray[n]->rgb) ) * pow(vectV.dotProduct(vectR), 19);
			} // if

            hasShadow = false;
		} // for lights

		accumRGB = accumRGB + attenuationMultiplier.pairwiseProduct(rgbColor);
		//~ accumRGB = accumRGB + attenuationMultiplier.pairwiseProduct(rgbColor);
		//~ accumRGB.x = accumRGB.x + attenuationMultiplier.x*rgbColor.x;
		//~ accumRGB.y = accumRGB.y + attenuationMultiplier.y*rgbColor.y;
		//~ accumRGB.z = accumRGB.z + attenuationMultiplier.z*rgbColor.z;

		//~ std::cout << "Printing initial accumRGB vector " << std::endl;
		//~ accumRGB.printVector3D();
		//~ if (row > 12 && row < 14 && column > 0 &&  column < 7) {
			//~ std::cout << row << " " << column <<  " row, column " << reflectionLevel << " reflectionLevel " << std::endl;
			//~ std::cout << ptrSphTest << " ptrSphTest" << std::endl;
			//~ surfaceNormal.printVector3D();
			//~ ptrRay->printRayAttributes();
		//~ }

		if (reflectionLevel > 0){
			//~ std::cout << reflectionLevel << " reflectionLevel " << std::endl;
			dInv = ptrRay->d*(-1); // Uinv // to accomodate for the direction convention of incident ray
			newD = (surfaceNormal*2*surfaceNormal.dotProduct(dInv)) - dInv; // spR i.e., R which is newD
			// Normalize only if it's not a zero vector (case: surface normal matches with the dInv vector)
			if (!(newD.x == 0 && newD.y == 0 && newD.z == 0)){
				newD.normalize(); // refR // Important!
			}

			// ==========================================================
			//~ //Ray RNew(ptrRay->best_intersection_sphere.x, ptrRay->best_intersection_sphere.y, ptrRay->best_intersection_sphere.z, newD.x, newD.y, newD.z);
			//~ //Ray* ptrRayNew = &RNew;

			ptrRay->eye.x = ptrRay->best_intersection_sphere.x;
			ptrRay->eye.y = ptrRay->best_intersection_sphere.y;
			ptrRay->eye.z = ptrRay->best_intersection_sphere.z;

			ptrRay->d.x = newD.x;
			ptrRay->d.y = newD.y;
			ptrRay->d.z = newD.z;

			//~ if (row > 12 && row < 14 && column > 0 &&  column < 7) {
				//~ std::cout << reflectionLevel << " reflectionLevel " << std::endl;
				//~ std::cout << row << " " << column <<  " row, column" << std::endl;
				//~ std::cout << ptrSphTest << " ptrSphTest" << std::endl;
				//~ dInv.printVector3D();
				//~ newD.printVector3D();
				//~ ptrRay->eye.printVector3D();
				//~ ptrRay->d.printVector3D();
				//~ ptrRay->printRayAttributes();
			//~ }

			// NOTE this attenuation step: IMPORTANT!!!!!
			attenuationMultiplierNew = ptrRay->ptr_bestSph->kr.pairwiseProduct(attenuationMultiplier);
			//~ attenuationMultiplier.x = ptrRay->ptr_bestSph->kr.x*attenuationMultiplier.x;
			//~ attenuationMultiplier.y = ptrRay->ptr_bestSph->kr.y*attenuationMultiplier.y;
			//~ attenuationMultiplier.z = ptrRay->ptr_bestSph->kr.z*attenuationMultiplier.z;

			rayTracer(ptrRay, row, column, ptr_sphArray, nSpheres, intersectSphere, ptr_modelArray, nModels, intersectModel, ptr_lightArray, nLightSources, ambient, accumRGB, attenuationMultiplierNew, reflectionLevel-1);
			// ==========================================================
		} // if level > 0

	} // if (ptrSphTest)


	else if (ptrModelTest && ptrRay->objectType == 2){
		// *************** Amibient reflection *********************
        //~ ptrRay->ptr_bestModel->ka.printVector3D();
        unsigned bestFaceMaterialIndex = ptrRay->ptr_bestModel->facelist_material_index[(unsigned)ptrRay->bestFaceIndex];
        rgbColor = ambient.pairwiseProduct(ptrRay->ptr_bestModel->ka[bestFaceMaterialIndex]);  // mtlchange!
        //~ std::cout << "Printing initial RGB vector " << std::endl;
        //~ rgbColor.printVector3D();

        // 2017
        // std::cout << "before ptrRay->bestNormal: " << std::endl;
        // ptrRay->bestNormal.printVector3D();
        toEye = ptrRay->eye - ptrRay->best_intersection_model;
        
        // commented to check if the vertices A B C follow a specific order for PA4 2018
        // If I comment this check then the shadows does not work. So keep this check!
        if (toEye.dotProduct(ptrRay->bestNormal) < 0){
          ptrRay->bestNormal = ptrRay->bestNormal*(-1);
          // std::cout << "after ptrRay->bestNormal: " << std::endl;
          // ptrRay->bestNormal.printVector3D();
        }
        
        //ptrRay->bestNormal.printVector3D("Model:: surface normal = ");

		if (row >= 512-RN2 && row <= 512-RN1 && column >= CN1 && column <= CN2) {
				std::cout << "\n\n" << 512-row-1 << " " << column <<  " row, column ================================" << std::endl;
		}
		
        for (unsigned n=0; n<nLightSources; n++){
            vectL = ptr_lightArray[n]->position - ptrRay->best_intersection_model;
            vectL.normalize(); // toL
            //~ std::cout << "Printing best intersection vector ---> " << n << std::endl;
            //~ ptrRay->best_intersection_model.printVector3D();
            //ptr_lightArray[n]->position.printVector3D();
            //////////vectL.printVector3D();

            bool hasShadow = false;
            // // cs410 A4 2017 SHADOWS ----------------------------------------
            // // check if toL goes through any other objects. If there is an object between light source and point of intersection
            // // then it means there should be a shadow. Do not add light source's contribution in this case.

            // // create a ray of toL vector
            // bool ifShadow = false;
            // Ray toLRay(ptrRay->best_intersection_model.x, ptrRay->best_intersection_model.y, ptrRay->best_intersection_model.z, vectL.x, vectL.y, vectL.z);
            // Ray* ptr_toLRay = &toLRay;
            // int toLIntersectSphere = 2;
            // int toLIntersectModel = 2;
            // // check for intersection with other objects
            // if (nSpheres != 0) {
            //     for (unsigned n=0; n<nSpheres; n++){
            //         // intersect values: 2 => default value 1 => return 1 from ray->raySphereIntersectionTest meaning intersected
            //         // 0 => return 0 from ray->raySphereIntersectionTest meaning no intersection
            //         toLIntersectSphere = ptr_toLRay->raySphereIntersectionTest(ptr_sphArray[n]);
            //         if (toLIntersectSphere == 1){
            //             ifShadow = true;
            //             break;
            //         }
            //     }
            // }
            // if (nModels != 0) {
            //     for (unsigned n=0; n<nModels; n++){
            //         toLIntersectModel = ptr_toLRay->rayModelIntersectionTest(ptr_modelArray[n]);
            //         if (toLIntersectModel == 1){
            //             ifShadow = true;
            //             break;
            //         }
            //     }
            // }

            hasShadow = shadow(ptrRay->best_intersection_model, ptr_lightArray[n]->position, nSpheres, ptr_sphArray, nModels, ptr_modelArray, debugFlag);
			
			if (row >= 512-RN2 && row <= 512-RN1 && column >= CN1 && column <= CN2) {
				//std::cout << row << " " << column <<  " row, column ================================" << std::endl;
				std::cout << n << " nth light source" << std::endl;
				ptrRay->bestNormal.printVector3D("ptrRay->bestNormal");
				vectL.printVector3D("vectL");
				std::cout << ptrModelTest << " ptrModelTest" << std::endl;
				std::cout << hasShadow << " hasShadow" << std::endl;
				std::cout << (ptrRay->bestNormal.dotProduct(vectL) > 0.000001) << " dotp>0?" << std::endl;
				if (ptrRay->bestNormal.dotProduct(vectL) > 0.000001 && hasShadow == false){
				std::cout << " Diffuse YES =========================================================" << std::endl;}
			}

			if (ptrRay->bestNormal.dotProduct(vectL) > 0.000001 && hasShadow == false){
				//std::cout << " Diffuse YES =========================================" << std::endl;
				// **************************** Diffuse reflection **************************************
				rgbColor = rgbColor + ( ptrRay->ptr_bestModel->kd[bestFaceMaterialIndex].pairwiseProduct(ptr_lightArray[n]->rgb) ) * ptrRay->bestNormal.dotProduct(vectL); // mtlchange!
				//~ std::cout << "Printing diffuse RGB vector " << n << std::endl;
				//~ rgbColor.printVector3D();

				vectV = ptrRay->eye - ptrRay->best_intersection_model;
				//////////vectV = pt8_cm->eye - ray->best_intersection;
				vectV.normalize(); // toC

				vectR = (ptrRay->bestNormal*2*ptrRay->bestNormal.dotProduct(vectL)) - vectL;
				//~ std::cout << "Printing r vector " << n << std::endl;
				//~ vectR.printVector3D();
				// **************************** Specular reflection (Phong) *****************************
				//~ std::cout << ptrRay->ptr_bestModel->Ns << " ptrRay->ptr_bestModel->Ns" << std::endl;
				//~ powerNs = ptrRay->ptr_bestModel->Ns;
				//~ std::cout << powerNs << " powerNs" << std::endl;

				/////IMPORTANT --------->// ????? ptrRay->ptr_bestModel->Ns needs to be typecasted to (int)! ????? otherwise it gives -nan -nan -nan as output!
				rgbColor = rgbColor + ( ptrRay->ptr_bestModel->ks[bestFaceMaterialIndex].pairwiseProduct(ptr_lightArray[n]->rgb) ) * pow((vectV.dotProduct(vectR)), (int)ptrRay->ptr_bestModel->Ns[bestFaceMaterialIndex]);  // mtlchange!
				//~ std::cout << "Printing specular RGB vector "  << n << std::endl;
				//~ rgbColor.printVector3D();
			  } // if
			
			if (row >= 512-RN2 && row <= 512-RN1 && column >= CN1 && column <= CN2) {
				rgbColor.printVector3D("rgbColor");
			}
              hasShadow = false;
		} // for lights

		/*if (row >= 512-186 && row <= 512-184 && column >= 297 && column <= 304) {
        //std::cout << "Printing final RGB vector " << std::endl;
			std::cout << "row, column ##### " << 512-row-1 << " " << column << std::endl;
			rgbColor.printVector3D("################### final RGB vector ");
		}*/
		if (row >= 512-RN2 && row <= 512-RN1 && column >= CN1 && column <= CN2) {
        //std::cout << "Printing final RGB vector " << std::endl;
			std::cout << "row, column ##### " << 512-row-1 << " " << column << std::endl;
			rgbColor.printVector3D("################### final RGB vector ");
		}
		//if(1){rgbColor.printVector3D("################### final RGB vector ");}


		accumRGB = accumRGB + attenuationMultiplier.pairwiseProduct(rgbColor);

		if (reflectionLevel > 0){
			//~ std::cout << reflectionLevel << " reflectionLevel " << std::endl;
			dInv = ptrRay->d*(-1); // Uinv // to accomodate for the direction convention of incident ray
			newD = (ptrRay->bestNormal*2*ptrRay->bestNormal.dotProduct(dInv)) - dInv; // spR i.e., R which is newD
			// Normalize only if it's not a zero vector (case: surface normal matches with the dInv vector)
			if (!(newD.x == 0 && newD.y == 0 && newD.z == 0)){
				newD.normalize(); // refR // Important!
			}

			ptrRay->eye.x = ptrRay->best_intersection_model.x;
			ptrRay->eye.y = ptrRay->best_intersection_model.y;
			ptrRay->eye.z = ptrRay->best_intersection_model.z;

			ptrRay->d.x = newD.x;
			ptrRay->d.y = newD.y;
			ptrRay->d.z = newD.z;

			//~ if (row > 12 && row < 14 && column > 0 &&  column < 7) {
				//~ std::cout << reflectionLevel << " reflectionLevel " << std::endl;
				//~ std::cout << row << " " << column <<  " row, column" << std::endl;
				//~ std::cout << ptrSphTest << " ptrSphTest" << std::endl;
				//~ dInv.printVector3D();
				//~ newD.printVector3D();
				//~ ptrRay->eye.printVector3D();
				//~ ptrRay->d.printVector3D();
				//~ ptrRay->printRayAttributes();
			//~ }

			// NOTE this attenuation step: IMPORTANT!!!!!
			attenuationMultiplierNew = ptrRay->ptr_bestModel->kr.pairwiseProduct(attenuationMultiplier);
			//~ attenuationMultiplier.x = ptrRay->ptr_bestSph->kr.x*attenuationMultiplier.x;
			//~ attenuationMultiplier.y = ptrRay->ptr_bestSph->kr.y*attenuationMultiplier.y;
			//~ attenuationMultiplier.z = ptrRay->ptr_bestSph->kr.z*attenuationMultiplier.z;

			rayTracer(ptrRay, row, column, ptr_sphArray, nSpheres, intersectSphere, ptr_modelArray, nModels, intersectModel, ptr_lightArray, nLightSources, ambient, accumRGB, attenuationMultiplierNew, reflectionLevel-1);
			// ==========================================================
		} // if level > 0

	} // if ptrModelTest else if

}


void rayTracerDepth(Camera* ptr_cm, Ray* ptrRay, int& row, int& column, Sphere** ptr_sphArray, unsigned& nSpheres, int& intersectSphere, Model** ptr_modelArray, unsigned& nModels, int& intersectModel, double& rgbDepth){
    // Go through every pixel and find the depth of at that pixel.

    unsigned objectType; // 1 => sphere; 2 => model

    // A2_2017
    Vector3D<double> pix_eye_diff;  // stores diff vector between eye and the pixel coordiante

    Sphere* ptrSphTest = 0;  // Initialize to zero. Very Important! (in case sphere is not present)
    Model* ptrModelTest = 0; // Initialize to zero. Very Important! (in case model is not present)
    //~ int intersectSphere = 2;

    //~ std::cout << nSpheres << " nSpheres" << std::endl;
    if (nSpheres != 0) {
        ptrSphTest = rayToWhichSphere(ptrRay, ptr_sphArray, nSpheres, intersectSphere);
    }

    if (nModels != 0) {
        ptrModelTest = rayToWhichModel(ptrRay, ptr_modelArray, nModels, intersectModel);
        // std::cout << ptrModelTest << " ptrModelTest (0==no object; other==that pointed object intersected)" << std::endl;
    }

    objectType = ptrRay->findBestIntersectionObject();
    // return 1 => sphere; return 2 => .obj object; return 0 => no intersection
    // std::cout << ptrRay->objectType << " ptrRay->objectType " << std::endl;

    //~ if (row > 5 && row < 9 && column > 5 &&  column < 9) {
        //~ std::cout << row << " " << column <<  " row, column" << std::endl;
        //~ std::cout << ptrSphTest << " ptrSphTest" << std::endl;
        //~ surfaceNormal.printVector3D();
    //~ }

    if (ptrSphTest && ptrRay->objectType == 1){  // ptrSphTest will be 0 if there are no spheres; as it's initialized to 0.
        // std::cout << "****************************calculating depth values for sphere" << std::endl;
        pix_eye_diff = ptrRay->pixel_xyz - ptrRay->eye;
        rgbDepth = ptrRay->bestT_sphere - pix_eye_diff.getMagnitude();

    } // if (ptrSphTest)


    else if (ptrModelTest && ptrRay->objectType == 2){
        // std::cout << "(((((((((((((((((((((((((((calculating depth values for model" << std::endl;
        // std::cout << "pixel_xyz, eye: " << std::endl;
        // ptrRay->pixel_xyz.printVector3D();
        // ptrRay->eye.printVector3D();
        pix_eye_diff = ptrRay->pixel_xyz - ptrRay->eye;
        // std::cout << "pix_eye_diff_magnitude, pix_eye_diff:" << pix_eye_diff.getMagnitude() << std::endl;
        // pix_eye_diff.printVector3D();
        rgbDepth = ptrRay->bestT_model - pix_eye_diff.getMagnitude();

        // std::cout << "Final depth value =========================" << rgbDepth << std::endl;

    } // if ptrModelTest else if
}


void renderImage(Camera* ptr_cm, Sphere** ptr_sphArray, unsigned& nSpheres, Light** ptr_lightArray, unsigned& nLightSources, Model** ptr_modelArray, unsigned& nModels, Vector3D<double>& ambient, unsigned recursionLevel, Vector3D<double>& rgbColor, std::string outputFilePpm){

    ptr_cm->ptr_depthMat = new Matrix2D<double>(ptr_cm->res_v, ptr_cm->res_u, -1); // initialize to -1 so that we could distinguish from required positive values

    Ray R(ptr_cm->eye.x, ptr_cm->eye.y, ptr_cm->eye.z, 0, 0, 0);
    Ray* ptrRay = &R;

    Sphere* ptrSphTest = 0;  // Initialize to zero. Very Important! (in case sphere is not present)
    Model* ptrModelTest = 0; // Initialize to zero. Very Important! (in case model is not present)
    double powerNs;

    // intersect values: 2 => default value 1 => return 1 from ray->raySphereIntersectionTest meaning intersected
    // 0 => return 0 from ray->raySphereIntersectionTest meaning no intersection
    int intersectSphere = 2;
    int intersectModel = 2;
    int sp0 = 0;
    int sp1 = 0;
    int sp00 = 0;
    int sp11 = 0;
    std::cout << " check ----------------------------------  "  << ptr_cm->eye.x << std::endl;
    ptrRay->printRayAttributes();

    //std::cin.get();
//    ray->eye.x = ptr_cm->eye.x;
//    ray->eye.y = ptr_cm->eye.y;
//    ray->eye.z = ptr_cm->eye.z;

    //std::cout << " check ----------------------------------  " << std::endl;

//    Vector3D<double> d, tImage, intersection; // d is the unit ray // tImage is from eye to pixel
//    double scalarC, scalarV, scalarD, ScalarTImage, scalarTEye;
//    scalarC = (ptr_sphArray[0]->center - ptr_cm->eye).getMagnitude();

	// A5 additions ---------->
	//~ unsigned reflectionLevel = 6;
	Vector3D<double> attenuationMultiplier(1.0,1.0,1.0); // NOTE: intitalized to all ones
	Vector3D<double> accumRGB;

    Vector3D<double> surfaceNormal, surfaceNormalModel, vectL, vectR, vectV; // ambient,
    //ambient.x = 0; ambient.y = 0; ambient.z = 0;
    std::cout << "Ambient RGB vector " << std::endl;
    ambient.printVector3D("");

    //~ std::ofstream op1;
    //~ op1.open("ptr_bool.txt", std::ofstream::out);

    //~ std::ofstream op2;
    //~ op2.open("model_checks.txt", std::ofstream::out);

    /////////////////////////////////////////////////////////////////////////
    std::ofstream op_ppm;
    op_ppm.open(outputFilePpm.c_str(), std::ofstream::out);

    // PPM HEADER -----------------------------------------------------
    op_ppm << "P3" << '\n';
    op_ppm << ptr_cm->res_u << " " << ptr_cm->res_v << " " << "255" << '\n';

    //std::cin.get();
    /////////////////////////////////////////////////////////////////////////
    std::cout << "Printing initial RGB vector " << std::endl;
    rgbColor.printVector3D("");

    for(int row=ptr_cm->res_v-1; row >=0 ; row--){
        for(int column=0; column < ptr_cm->res_u; column++){
          // std::cout << "pixel (r, c) -------> (" << ptr_cm->res_v-1 - row << ", " << column << ") " << std::endl;
//////////            std::cout << "RGB vector -----------> (r,c) -------> (" << ptr_cm->res_v-1 - row << ", " << column << ") "<< std::endl;
//////////            rgbColor.printVector3D();

			// -------------------------------------------------------
			// Looks redundant, but it ain't !!! Solved first issues in a5. Before, only first intersecting ray was having correct rgb, next wasn't even intersecting when it should have!
			ptrRay->eye.x = ptr_cm->eye.x;
			ptrRay->eye.y = ptr_cm->eye.y;
			ptrRay->eye.z = ptr_cm->eye.z;
			// -------------------------------------------------------

            ptrRay->d.x = ptr_cm->ptr_rayMatX->getElement(row,column);
            ptrRay->d.y = ptr_cm->ptr_rayMatY->getElement(row,column);
            ptrRay->d.z = ptr_cm->ptr_rayMatZ->getElement(row,column);

            //R.printRayAttributes();
            //ptrRay->printRayAttributes();
            //std::cout << ptrRay << " ptrRay" <<std::endl;
            //std::cin.get();

            if (nSpheres == 0 && nModels == 0) {
                rgbColor.x = 0; // 255;
                rgbColor.y = 0; // 255;
                rgbColor.z = 0; // 255;
                op_ppm << rgbColor.x << " " << rgbColor.y << " " << rgbColor.z << " ";
            } // if nObjects is zero

            else {

                if (nSpheres != 0 || nModels != 0) {
                    //~ ptrSphTest = rayToWhichSphere(ptrRay, ptr_sphArray, nSpheres, intersectSphere);
                    //ptrSphTest = rayToWhichSphere(ptrRay, ptr_sphArray, nSpheres, intersect);
                    //~ op1 << intersectSphere << " " << ptrRay->bestT_sphere << " (" << std::setprecision(2) << std::fixed << ptr_cm->res_v-1 - row << ", " << column << ") , (" << std::setprecision(3) << std::fixed << ptrRay->d.x << std::setprecision(3) << std::fixed << ", " << ptrRay->d.y <<  ", " << std::setprecision(3) << std::fixed << ptrRay->d.z << ")" << '\t';
                    //~ std:: cout << ptrSphTest << " <<<<<----- ptrSphTest "<< std::endl;
                //~ }
                //std::cout << ptrRay->bestT_sphere << " ptrRay->bestT_sphere" << std::endl;

                //~ if (nModels != 0) {
                    //~ ptrModelTest = rayToWhichModel(ptrRay, ptr_modelArray, nModels, intersectModel);
                    //~ op2 << intersectModel << " " << ptrRay->bestT_model << " " << ptrRay->ptr_bestModel << " (" << std::setprecision(2) << std::fixed << ptr_cm->res_v-1 - row << ", " << column << ") , (" << std::setprecision(3) << std::fixed << ptrRay->d.x << std::setprecision(3) << std::fixed << ", " << ptrRay->d.y <<  ", " << std::setprecision(3) << std::fixed << ptrRay->d.z << ")" << '\t';
                //~ }

                //std::cout << ptrRay->bestT_model << " ptrRay->bestT_model" << std::endl;


                //ptrRay->printRayAttributes();
                //ptrRay->bestT_sphere < ptrRay->bestT_model && ptrRay->bestT_sphere >= 0

                ///if (nSpheres != 0) {

                        //if (ptrSphTest && ptrRay->bestT_sphere < ptrRay->bestT_model && ptrRay->bestT_sphere >= 0){ // if not NULL
                        //~ if (ptrSphTest && ptrRay->bestT_sphere >= 0 && ptrRay->bestT_sphere < ptrRay->bestT_model){ // if not NULL
                            //~ sp1++;
                            //~ std::cout << ptr_cm->res_v-1 - row << ", " << column << "  <<<---- row, column" << std::endl;
                            //std::cout <<[* "Printing initial RGB vector ----- " << std::endl;
                            //rgbColor.printVector3D();
							//~ std::cout << row << " " << column <<  " row, column" << std::endl;
							rayTracer(ptrRay, row, column, ptr_sphArray, nSpheres, intersectSphere, ptr_modelArray, nModels, intersectModel, ptr_lightArray, nLightSources, ambient, accumRGB, attenuationMultiplier, recursionLevel);

                            // A5 additions ----------------------------
                            //~ accumRGB = accumRGB + ();


                            //

                            // following solved the problem with the specular reflection term of rgbColor
                            accumRGB.x = fmax(0, fmin(255, floor(255*accumRGB.x)));
                            accumRGB.y = fmax(0, fmin(255, floor(255*accumRGB.y)));
                            accumRGB.z = fmax(0, fmin(255, floor(255*accumRGB.z)));
                            op_ppm << accumRGB.x << " " << accumRGB.y << " " << accumRGB.z << " ";

                        //~ } // if ptrSphTest

//                        else {
//                            sp0++;
//                            rgbColor.x = 0; // 255;
//                            rgbColor.y = 0; // 255;
//                            rgbColor.z = 0; // 255;
//                            op_ppm << rgbColor.x << " " << rgbColor.y << " " << rgbColor.z << " ";
//                        }

			} // if nSpheres != 0 added newly a5
                ///} // if nSpheres is not zero

                ///if (nModels != 0){
// ------------------------
/*
                        else if (ptrModelTest && ptrRay->bestT_model >= 0 && ptrRay->bestT_model <= ptrRay->bestT_sphere){ // note <= in third condition
							// do models processing add it on to rgbColor calculated for
                            sp11++;
                            //~ std::cout << "inside sp11++ " << std::endl;

                            //std::cout << ptrRay << " ptrRay after best find" <<std::endl;
                            //std::cin.get();

                            // *************** Amibient reflection *********************
                            //~ ptrRay->ptr_bestModel->ka.printVector3D();
                            rgbColor = ambient.pairwiseProduct(ptrRay->ptr_bestModel->ka);
                            //~ std::cout << "Printing initial RGB vector " << std::endl;
                            //~ rgbColor.printVector3D();

                            for (unsigned n=0; n<nLightSources; n++){
                                vectL = ptr_lightArray[n]->position - ptrRay->best_intersection_model;
                                vectL.normalize(); // toL
                                //~ std::cout << "Printing best intersection vector ---> " << n << std::endl;
                                //~ ptrRay->best_intersection_model.printVector3D();
                                //ptr_lightArray[n]->position.printVector3D();
            //////////                    vectL.printVector3D();

                                if (ptrRay->bestNormal.dotProduct(vectL) > 0){
                                    // **************************** Diffuse reflection **************************************
                                    rgbColor = rgbColor + ( ptrRay->ptr_bestModel->kd.pairwiseProduct(ptr_lightArray[n]->rgb) ) * ptrRay->bestNormal.dotProduct(vectL);
									//~ std::cout << "Printing diffuse RGB vector " << n << std::endl;
									//~ rgbColor.printVector3D();

                                    vectV = ptr_cm->eye - ptrRay->best_intersection_model;
            //////////                        vectV = ptr_cm->eye - ray->best_intersection;
                                    vectV.normalize(); // toC

                                    //~ std::cout << "Printing v vector " << n << std::endl;
									//~ vectV.printVector3D();

                                    vectR = (ptrRay->bestNormal*2*ptrRay->bestNormal.dotProduct(vectL)) - vectL;
                                    //~ std::cout << "Printing r vector " << n << std::endl;
									//~ vectR.printVector3D();
                                    // **************************** Specular reflection (Phong) *****************************
                                    //~ std::cout << ptrRay->ptr_bestModel->Ns << " ptrRay->ptr_bestModel->Ns" << std::endl;
                                    //~ powerNs = ptrRay->ptr_bestModel->Ns;
                                    //~ std::cout << powerNs << " powerNs" << std::endl;

           /////IMPORTANT --------->// ????? ptrRay->ptr_bestModel->Ns needs to be typecasted to (int)! ????? otherwise it gives -nan -nan -nan as output!
                                    rgbColor = rgbColor + ( ptrRay->ptr_bestModel->ks.pairwiseProduct(ptr_lightArray[n]->rgb) ) * pow((vectV.dotProduct(vectR)), (int)ptrRay->ptr_bestModel->Ns);
                                    //~ std::cout << "Printing specular RGB vector "  << n << std::endl;
									//~ rgbColor.printVector3D();
                                } // if
                            } // for

                            //~ std::cout << "Printing final RGB vector " << std::endl;
                            //~ rgbColor.printVector3D();

                            // following solved the problem with the specular reflection term of rgbColor
                            rgbColor.x = fmax(0, fmin(255, floor(255*rgbColor.x)));
                            rgbColor.y = fmax(0, fmin(255, floor(255*rgbColor.y)));
                            rgbColor.z = fmax(0, fmin(255, floor(255*rgbColor.z)));
                            op_ppm << rgbColor.x << " " << rgbColor.y << " " << rgbColor.z << " ";

							//~ rgbColor.x = 128; // 255;
                            //~ rgbColor.y = 128; // 255;
                            //~ rgbColor.z = 128; // 255;
                            //~ op_ppm << rgbColor.x << " " << rgbColor.y << " " << rgbColor.z << " ";
                        } // if ptrModelTest else if


                ///} // if nModels != 0
*/
// ----------------

// have to do something about this else
/*
						else {
                            sp0++;
                            rgbColor.x = 0; // 255;
                            rgbColor.y = 0; // 255;
                            rgbColor.z = 0; // 255;
                            op_ppm << rgbColor.x << " " << rgbColor.y << " " << rgbColor.z << " ";
                        } // sphere model or else
*/
            } // else if nObjects (nSpheres or nModels) is not zero



/*
            //******************************************** check if center is infront of the image plane ******************************************
            scalarV = (ptr_sphArray[0]->center - ptr_cm->eye).dotProduct(d);
            //std::cout << " cl and v and r " << cl << ", " << v << ", " << (ptr_ply1->radiusOfSphere) << std::endl;
            //int ch = getchar();
            scalarD = (ptr_sphArray[0]->radius)*(ptr_sphArray[0]->radius) - (scalarC*scalarC-scalarV*scalarV); // d^2

            if(scalarD<0){
                //tOfAPixel[face] = -2.55;
                ptr_cm->ptr_depthMat->setElement(row, column, -3);
                rgbColor.x = 0; // 255;
                rgbColor.y = 0; // 255;
                rgbColor.z = 0; // 255;
            }

            else {
                // distance from intersection to eye
                scalarTEye = fabs(scalarV)-fabs(sqrt(scalarD)); // = (v-d)
                ptr_cm->ptr_depthMat->setElement(row, column, scalarTEye); // Q-E = (scalarV-scalarD)*d

                /*
                // changes tdepth values from eye to image plane ******************************************************************************
                tImage.x = ptr_cm->ptr_pixMatX->getElement(row, column) - ptr_cm->eye.x;
                tImage.y = ptr_cm->ptr_pixMatY->getElement(row, column) - ptr_cm->eye.y;
                tImage.z = ptr_cm->ptr_pixMatZ->getElement(row, column) - ptr_cm->eye.z;
                ScalarTImage = ptr_cm->ptr_depthMat->getElement(row, column) - tImage.getMagnitude();
                ptr_cm->ptr_depthMat->setElement(row, column, ScalarTImage);
                /

                // Find coordinates of the point of intersection
                intersection = ptr_cm->eye + d*scalarTEye; // d is the unit ray
                //intersection.printVector3D();


//////////                calculateRaySphereRGB(ptr_cm, ptr_sphArray[0], nSpheres, ptr_lightArray, nLightSources, intersection, rgbColor);
                //rgbColor = rgbColor*255;
//////////                // following solved the problem with the specular reflection term of rgbColor
//////////                rgbColor.x = fmax(0, fmin(255, floor(255*rgbColor.x)));
//////////                rgbColor.y = fmax(0, fmin(255, floor(255*rgbColor.y)));
//////////                rgbColor.z = fmax(0, fmin(255, floor(255*rgbColor.z)));
            } // else

*/
//////////        op1 << std::setprecision(3) << std::fixed << ptr_cm->ptr_depthMat->getElement(row, column) << '\t';
//////////            op_ppm << rgbColor.x << " " << rgbColor.y << " " << rgbColor.z << " ";

            rgbColor.x=0;
            rgbColor.y=0;
            rgbColor.z=0;

            // A5 additions ------------>
            accumRGB.x=0;
            accumRGB.y=0;
            accumRGB.z=0;

            ptrSphTest = NULL;
        } // for column

        //~ op1 << '\n';
        //~ op2 << '\n';
        op_ppm << '\n'; // newline after each row in PPM file

    } // for row

    std::cout << sp0 << " sp0 - no intersection " << sp1 << " sp1 intersected " << sp11 << " sp11 model intersects"<< std::endl;
    //ptr_cm->ptr_depthMat->printMatrix();
    //~ op1.close();
    //~ op2.close();
    op_ppm << '\n'; // to match the given PPM files
    op_ppm.close();
}



void renderDepthImage(Camera* ptr_cm, Sphere** ptr_sphArray, unsigned& nSpheres, Model** ptr_modelArray, unsigned& nModels, Vector3D<double>& rgbColor, std::string outputFilePpm){
    ptr_cm->ptr_depthMat = new Matrix2D<double>(ptr_cm->res_v, ptr_cm->res_u, -1.2345); // initialize to -1 so that we could distinguish from required positive values

    Ray R(ptr_cm->eye.x, ptr_cm->eye.y, ptr_cm->eye.z, 0, 0, 0);  // why dx, dy, dz = 0?
    Ray* ptrRay = &R;

    Sphere* ptrSphTest = 0;  // Initialize to zero. Very Important! (in case sphere is not present)
    Model* ptrModelTest = 0; // Initialize to zero. Very Important! (in case model is not present)

    // A1_2017
    double rgbDepth = -2.2345;  // a double storing the depth value at every pixel  // set it to a negative value // solves the "running" color issue

    int intersectSphere = 2;
    int intersectModel = 2;
    std::cout << " check eye vector ----------------------------------  "  << ptr_cm->eye.x << std::endl;
    ptrRay->printRayAttributes();

    //~ std::ofstream op1;
    //~ op1.open("ptr_bool.txt", std::ofstream::out);

    //~ std::ofstream op2;
    //~ op2.open("model_checks.txt", std::ofstream::out);

    /////////////////////////////////////////////////////////////////////////
    std::cout << "============================ Writing the output PPM file ===========================" << std::endl;
    std::ofstream op_ppm;
    op_ppm.open(outputFilePpm.c_str(), std::ofstream::out);
    // PPM HEADER -----------------------------------------------------
    op_ppm << "P3" << '\n';
    op_ppm << ptr_cm->res_u << " " << ptr_cm->res_v << " " << "255" << '\n';

    /////////////////////////////////////////////////////////////////////////
    std::cout << "Printing initial RGB vector " << std::endl;
    rgbColor.printVector3D("");


    if (nSpheres == 0 && nModels == 0) {
        std::cout << "No sphere or model present!!!!! " << std::endl;
        for(int row=ptr_cm->res_v-1; row >=0 ; row--){
            for(int column=0; column < ptr_cm->res_u; column++){
                op_ppm << 239 << " " << 239 << " " << 239 << " ";
            }  // for column
        }  // for row
    }  // if nSpheres and nModels == 0

    // Go through every pixel and store the depth values.
    else if (nSpheres != 0 || nModels != 0) {
        // -------------------------------------------------------
        for(int row=ptr_cm->res_v-1; row >=0 ; row--){
            for(int column=0; column < ptr_cm->res_u; column++){
                // std::cout << "(row, column) = (" << row << ", " << column << ")" << std::endl;

                // Looks redundant, but it ain't !!! Solved first issues in a5. Before, only first intersecting ray was having correct rgb, next wasn't even intersecting when it should have!
                ptrRay->eye.x = ptr_cm->eye.x;
                ptrRay->eye.y = ptr_cm->eye.y;
                ptrRay->eye.z = ptr_cm->eye.z;
                // std::cout << "Printing eye vector of the ray " << std::endl;
                // ptrRay->eye.printVector3D();  // correct for a2
                // -------------------------------------------------------

                ptrRay->d.x = ptr_cm->ptr_rayMatX->getElement(row,column);
                ptrRay->d.y = ptr_cm->ptr_rayMatY->getElement(row,column);
                ptrRay->d.z = ptr_cm->ptr_rayMatZ->getElement(row,column);
                // std::cout << "Printing d vector of the ray " << std::endl;
                // ptrRay->d.printVector3D();

                // A2_2017: added a xyz coordinate of the pixel (pixel_xyz) from which the ray passes.
                // Needed for calculating depth which will be bestT_* - (pixel_xyz - eye)
                // bestT_* is from eye to best intersection, we want depth from pixel location to best intersection
                // So, subtract from bestT_*, the distance between eye and pixel coordinate (pixel_xyz - eye)
                ptrRay->pixel_xyz.x = ptr_cm->ptr_pixMatX->getElement(row,column);
                ptrRay->pixel_xyz.y = ptr_cm->ptr_pixMatY->getElement(row,column);
                ptrRay->pixel_xyz.z = ptr_cm->ptr_pixMatZ->getElement(row,column);
                // std::cout << "Printing pixel_xyz vector of the ray " << std::endl;
                // ptrRay->pixel_xyz.printVector3D();

                // stores the depths in the ptr_cm->depthMat, but for that assign rbbDepth to the ptr_cm->depthMat(r,c)
                rayTracerDepth(ptr_cm, ptrRay, row, column, ptr_sphArray, nSpheres, intersectSphere, ptr_modelArray, nModels, intersectModel, rgbDepth);
                ptr_cm->ptr_depthMat->setElement(row, column, rgbDepth);
                rgbDepth = -3.2345;  // solves the "running" color issue

                // ptrSphTest = NULL;
                // + other pointer to Null?
            } // for column

        } // for row

    }  // else if

    // Printing the depth matrix
    // std::cout << "Printing the depth matrix ========================= " << std::endl;
    // ptr_cm->ptr_depthMat->printMatrix();

    // Go through depth matrix and store the rgb of that depth to a ppm file.
    // Apply the formula to convert the depth values to an rgb color.
    generatePpmImage(ptr_cm, outputFilePpm);

}


int main(int argc, char *argv[]){
    //std::cout << "CS410 Assignment 5 by Gururaj Mulay" << std::endl;
    //std::cin.get();

    int nVertices = 0;
    int nFaces = 0;
    unsigned nLightSources = 0;
    unsigned nSpheres = 0;
    unsigned nModels = 0;
    std::string modelName = "NULL";

    Vector3D<double> rgbColor;
    Vector3D<double> ambient;

    unsigned recursionLevel;

    double dt;
    double total_time = 0;
    int start_sec_matrix, stop_sec_matrix;
    start_sec_matrix = clock();

    std::map<std::string,int> modelmap;

// CMD Inputs 1 --------------------------------------------------------------------------------------------------
    std::cout << "Executable name is -----> " << argv[0] << std::endl;
    std::cout << "Output PPM name is -----> " << argv[2] << std::endl;
    std::string inputFileCamera(argv[1]);
    //std::string inputFilePly(argv[2]);
    std::string outputFilePpm(argv[2]); // A5_2016
    //std::string inputFileScene(argv[4]);
    //std::string outputFilePpm = "opImage.ppm";


// Camera File Processing --------------------------------------------------------------------------------------------------
    Camera* ptr_cm = getCameraAttributes(inputFileCamera);

    ptr_cm->printCameraAttributes();
    ptr_cm->getWVector();
    ptr_cm->getUVector();
    ptr_cm->getVVector();
    ptr_cm->getImagePlaneCenterVector();
    ptr_cm->getPixeledImagePlane();
    std::cout << ptr_cm->res_v << " res_v | " <<  ptr_cm->res_u << " res_u" << std::endl;
//    ptr_cm->ptr_rayMatY->printMatrix();

//-----------------------------------------------------------------------------------------------------------------------------
    // find the total light sources in the scene file
    getNumberOfLightSourcesSpheresModels(inputFileCamera, nLightSources, nSpheres, nModels);

    std::cout << nLightSources << " <--- number of light sources" << std::endl;
    std::cout << nSpheres << " <--- number of spheres" << std::endl;
    std::cout << nModels << " <--- number of models" << std::endl;

// A5_2016
    getAmbientRGBVector(inputFileCamera, ambient);
    std::cout << "ambient vector is -----------------------> " << std::endl;
    ambient.printVector3D("ambient!!!!!");

    getRecursionLevel(inputFileCamera, recursionLevel);
    std::cout << recursionLevel  << " recursionLevel " << std::endl;

    // Important: Created an array lightSource of pointers of type Light*; Array pointer Light** - a double poiner
    Light** lightSourceArray = new Light*[nLightSources];
    //std::cout << lightSourceArray[0] << " | " << lightSourceArray[1] << " <--- ptrs to light sources (at init)" << std::endl;
    // Now create an object for every instance of the light source (nLightSources times)
    for(unsigned lightNumber=0; lightNumber<nLightSources; lightNumber++){
        lightSourceArray[lightNumber] = getLightAttributes(inputFileCamera, lightNumber);
        //lightSourceArray[lightNumber]->printLightAttributes();
    }
    //std::cout << size(lightSourceArray) << " <<<<<<<<------------------- " << std::endl;
    //std::cout << lightSourceArray[0] << " | " << lightSourceArray[1] << " <--- ptrs to light source (after object creation)" << std::endl;

//

  Sphere** sphereArray = new Sphere*[nSpheres];
  for(unsigned sphereNumber=0; sphereNumber<nSpheres; sphereNumber++){
    sphereArray[sphereNumber] = getSphereAttributes(inputFileCamera, sphereNumber);
    sphereArray[sphereNumber]->printSphereAttributes();
  }

	Model** modelArray = new Model*[nModels];
	for(unsigned modelNumber=0; modelNumber<nModels; modelNumber++){
		modelArray[modelNumber] = parseModelLine(inputFileCamera, modelNumber, modelName);
		modelArray[modelNumber]->getNumberOfVerticesAndFaces(modelArray[modelNumber]->modelName);
        modelArray[modelNumber]->getModelCommentBlock(modelArray[modelNumber]->modelName);
		if(nModels != 0) {
			modelArray[modelNumber]->getModelVerticesFaceAttributes(modelArray[modelNumber]->modelName); // also gets the modelMaterialName
			modelArray[modelNumber]->getModelMaterialAttributes(modelArray[modelNumber]->modelMaterialName);
			modelArray[modelNumber]->printModelAttributes();

		} // if (nModels != 0)
	}

  // A2_2017 commented for A3
  // Sphere and Model Processing -------------------------------------------------------------------------------------------------------------------
  // renderDepthImage(ptr_cm, sphereArray, nSpheres, modelArray, nModels, rgbColor, outputFilePpm);
  //generatePpmImageSphere(ptr_cm, outputFilePpm, rgbColor);


// A5_2016
// Sphere and Model Processing -------------------------------------------------------------------------------------------------------------------
    std::cout << sphereArray[0] << " sphere 0 pointer" << std::endl;
    renderImage(ptr_cm, sphereArray, nSpheres, lightSourceArray, nLightSources, modelArray, nModels, ambient, recursionLevel, rgbColor, outputFilePpm);
    //generatePpmImageSphere(ptr_cm, outputFilePpm, rgbColor);

//    Ray r1(1,2,3,4,5,6);
//    r1.printRayAttributes();
//



// Depth Values Calculation ---------------------------------------------------------------------------------------------------
    // calculateDepthValues(ptr_ply1, ptr_cm);
    //std::cout << "printing depth matrix after processing ----------------------------------> " << std::endl;
    //ptr_cm->ptr_depthMat->printMatrix();



//------------------------------------------------------------------------------------------------------------
//FINAL STEP of PRINTING PPM IMAGE
    // generatePpmImage(ptr_cm, outputFilePpm);



//------------------------------------------------------------------------------------------------------------
/*
    ptr_ply1->findMinElementArray();
    ptr_ply1->findMaxElementArray();
    ptr_ply1->getObjectCenter();
    ptr_ply1->getSphereRadius();
*/


// General Testing --------------------------------------------------------------------------------------------------
/*
    Matrix2D<double> MM(3, 5, 2.5);
    MM.printMatrix();
    //std::cout << MM(2, 3) << " <----- MM(2,3) " << std::endl;

    double num = 100000;
    double num1 = INFINITY;
    double inf = std::numeric_limits<double>::infinity();
    std::cout << inf << " <----- INFINITY! " << std::endl;
    double num2 = inf;

    if (num1==inf){
        std::cout << num2 << " <----- num equal to INFINITY! " << std::endl;
    }
*/

/* A1_2017
    // create folder named after driver file
    std::string driverdir = inputFileCamera.substr(0, inputFileCamera.size() - 4);
    int dircreated = system(("mkdir -p ./" + driverdir).c_str());
    std::cout << "dir created ??? " << dircreated << std::endl;

    if(nModels != 0) {
      for(unsigned modelNumber=0; modelNumber<nModels; modelNumber++){
        // modelArray[modelNumber]->printModelAttributes();
        // modelArray[modelNumber]->printTransformedModelVerticesFaces();

        std::string model_name = modelArray[modelNumber]->modelName;
        std::string outputFileObj;

        std::cout << "Object is =============================>" << model_name << std::endl;
        std::cout << modelArray[modelNumber]->firstCommentBlock << std::endl;

        // add model to modelmap if not exists or increase the count if exists
        if (modelmap.count(model_name)>0){  // if exists
          // std::cout << " is an element of modelmap.\n";
          // std::cout << modelmap[model_name] << std::endl;
          modelmap[model_name] = modelmap[model_name] + 1;  // increase count by 1
        }  // if

        else{
          // std::cout << " is not an element of modelmap.\n";
          modelmap[model_name] = 1;  // intialize count to 1
        }  // else


        int count = modelmap[model_name] - 1;
        // std::string strCount = std::to_string(count);
        std::stringstream ss;
        ss << count;
        std::string strCount = ss.str();

        // TODO make it work for more than 10 objects.
        outputFileObj = model_name.substr(0, model_name.size() - 4) + "_mw0" + strCount + ".obj";
        // std::cout << "The output file name is: " << outputFileObj << std::endl;


        // Write output files
        std::ofstream op_obj;
        std::string outpath = "./" + driverdir + "/" + outputFileObj;
        std::cout << "output file path ===" << outpath << std::endl;
        op_obj.open(outpath.c_str(), std::ofstream::out);

        // PPM HEADER -----------------------------------------------------
        op_obj << modelArray[modelNumber]->firstCommentBlock.substr(0, modelArray[modelNumber]->firstCommentBlock.size() - 1) << '\n';

        double matElement;
        unsigned i, j;
        std::cout << "Matrix is ------------------> " << std::endl;

        for (i=0; i<modelArray[modelNumber]->ptr_verticesMat->rows; i++){
          std::cout << "v ";
          op_obj << "v ";
          for (j=0; j<modelArray[modelNumber]->ptr_verticesMat->columns; j++){
            matElement = modelArray[modelNumber]->ptr_verticesMat->getElement(i, j);
            if (j==2){
              std::cout << std::setprecision(7) << std::fixed << matElement;
              op_obj << std::setprecision(7) << std::fixed << matElement;
            }
            else {
              std::cout << std::setprecision(7) << std::fixed << matElement << " ";
              op_obj << std::setprecision(7) << std::fixed << matElement << " ";
            }

          }  // for j
          std::cout << std::endl;
          op_obj << std::endl;
        }  // for i

        for (i=0; i<modelArray[modelNumber]->ptr_faceMat->rows; i++){
          std::cout << "f ";
          op_obj << "f ";
          for (j=0; j<modelArray[modelNumber]->ptr_faceMat->columns; j++){
            matElement = modelArray[modelNumber]->ptr_faceMat->getElement(i, j);
            // std::cout << "j ===" << j << std::endl;
            if (j%3 == 0){  // v index
              std::cout << (int)matElement << "/";
              op_obj << (int)matElement << "/";
            }  // if v
            if (j%3 == 1){  // don't write the vt index
              std::cout << "/";
              op_obj << "/";
            }  // if vt
            if (j%3 == 2){  // vn index
              if (j==8){  // last element
                std::cout << (int)matElement;
                op_obj << (int)matElement;
              }
              else {
                std::cout << (int)matElement << " ";
                op_obj << (int)matElement << " ";
              }
            }  // if vn

          }  // for j
          std::cout << std::endl;
          op_obj << std::endl;
        }  // for i


      }  // for every model
    }  // if

    // write out the transformed model coordinates

*/ // A1


// Exit --------------------------------------------------------------------------------------------------
// delete ptr_ply1;

// A5_2016
    delete ptr_cm;

    for(unsigned lightNumber=0; lightNumber<nLightSources; lightNumber++){
        delete lightSourceArray[lightNumber];
    }
    delete[] lightSourceArray;

    for(unsigned sphereNumber=0; sphereNumber<nSpheres; sphereNumber++){
        delete sphereArray[sphereNumber];
    }
    delete[] sphereArray;
//

    stop_sec_matrix = clock();
    dt  = (stop_sec_matrix-start_sec_matrix)/double(CLOCKS_PER_SEC)*1;
    std::cout << "total time " << dt << std::endl;

    //~ for(unsigned modelNumber=0; modelNumber<nModels; modelNumber++){
        //~ delete modelArray[modelNumber];
    //~ }
    //~ delete[] modelArray;

    std::cout << "Exiting with return 0 ..... " << std::endl;
    return 0;
} //main


// NOTES:
// mtlchange!
// This program will not work as expected if the order of materials mentioned in .obj is different from what .mtl specifies.
