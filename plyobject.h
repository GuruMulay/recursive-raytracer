#ifndef PLYOBJECT_H_INCLUDED
#define PLYOBJECT_H_INCLUDED


#include <sstream>
#include <fstream>

#include <iostream>
#include <algorithm>

#include "vector3d.h"
#include "matrix2d.h"

class PlyObject{
public:
    int n_axes, n_vertices, n_faces, max_polygonSides;
    double **vertices;
    double *maxElementArray;
    double *minElementArray;
    double *meanArray;
    double *stdDeviationArray;

    // model transformation attributes
    Vector3D<double> modelTranslation;
    Matrix2D<double>* ptr_rotationMat;

    Vector3D<double> object_center;
    double radiusOfSphere;
    //int n_polygonSides = 4;

    Matrix2D<double>* ptr_verticesMat;
    Matrix2D<unsigned>* ptr_faceMat; //unsigned because its either number of vertices per face or vertex number

    //double **rotationMatrix;
    //double vertices[m][n];
    //bool isPLY = false, isFormatASCII = false;

    PlyObject(int, int, int, int);

    void getPlyObjectVerticesFaceAttributes(std::string);

    void assignMatrixValue(int, int, double, int);
    void printCoordinates();
    //void printRotationMatrix();
    void findMaxElementArray();
    void findMinElementArray();
    void findMeanArray();

    void getObjectCenter();
    void getSphereRadius();

    void findStdDeviationArray();
    void objectCentering();
    void objectRounding();
    void printBoundingBox();
    ~PlyObject(); //PlyObject destructor

public:
    void allocateMatrix(){
		vertices = new double*[n_axes];
		for(int i = 0; i < n_axes; i++)
            {
                vertices[i] = new double[n_vertices];
            }
	}
}; // PlyObject class

PlyObject::PlyObject(int nAxes, int nVertices, int nFaces, int max_polySides){
    n_axes = nAxes;
    n_vertices = nVertices;
    n_faces = nFaces;
    max_polygonSides = max_polySides; // solved error!

    allocateMatrix();
    for(int i = 0; i < n_axes; i++)
        {
            for(int j = 0; j < n_vertices; j++)
                {
                vertices[i][j] = 1.0;
                }
        }
}// constructor for matrix of vertices


void PlyObject::assignMatrixValue(int nAxes, int nVertices, double value, int flag){
	double newValue = value;
    vertices[nAxes][nVertices] = newValue;
    //std::cout << "assigned to vertex matrix! \n";
} // constructor to assign matrix value

void PlyObject::printCoordinates(){
	for(int i = 0; i < n_axes; i++)
        {
            std::cout << "axis " << i << " vertices of the PlyObject --->" << std::endl;
            std::cout << "" << std::endl;
            for(int j = 0; j < n_vertices; j++)
            {
                std::cout << vertices[i][j] << "  " << std::endl;
            }
            std::cout << "" << std::endl;
        }
} // constructor for printing the coordinates of the vertices

void PlyObject::findMaxElementArray(){
    maxElementArray = new double[n_axes-1];
    for(int j=0; j<n_axes-1; j++){maxElementArray[j] = ptr_verticesMat->getElement(0,j);} // initialize array

    for(int j=0; j<(n_axes-1); j++)
        {
            for(int i=1; i<n_vertices; i++)
                {
                    if (maxElementArray[j] < ptr_verticesMat->getElement(i,j)){
                        maxElementArray[j] = ptr_verticesMat->getElement(i,j);
                    }
                }
        }
     std::cout << "max element array is --- " << maxElementArray[0] << ", " << maxElementArray[1] << ", " << maxElementArray[2] << std::endl;
    //std::cout << "max element is --- " << maxElementArray[2] << std::endl;
} // constructor to find maximum element between the coordinates

void PlyObject::findMinElementArray(){
    minElementArray = new double[n_axes-1];
    for(int j=0; j<n_axes-1; j++){minElementArray[j] = ptr_verticesMat->getElement(0,j);} // initialize array
    //std::cout << "min element array is --- " << minElementArray[0] << ", " << minElementArray[1] << ", " << minElementArray[2] << std::endl;

    for(int j=0; j<(n_axes-1); j++)
        {
            for(int i=1; i<n_vertices; i++)
                {
                    if (minElementArray[j] > ptr_verticesMat->getElement(i,j)){
                        minElementArray[j] = ptr_verticesMat->getElement(i,j);
                    }
                }
        }
    std::cout << "min element array is --- " << minElementArray[0] << ", " << minElementArray[1] << ", " << minElementArray[2] << std::endl;
    //std::cout << "min element is --- " << minElementArray[0] << std::endl;
} // constructor to find minimum element between the coordinates


void PlyObject::getObjectCenter(){
    object_center.x = 0.5*(maxElementArray[0]+minElementArray[0]);
    object_center.y = 0.5*(maxElementArray[1]+minElementArray[1]);
    object_center.z = 0.5*(maxElementArray[2]+minElementArray[2]);
    std::cout << "center is --- " << object_center.x << ", " << object_center.y << ", " << object_center.z << std::endl;
}

void PlyObject::getSphereRadius(){
    radiusOfSphere = std::max((maxElementArray[0]-object_center.x), (maxElementArray[1]-object_center.y));
    radiusOfSphere = std::max(radiusOfSphere,(maxElementArray[2]-object_center.z));
    std::cout << "max radius is --- " << radiusOfSphere << std::endl;
}


 void PlyObject::printBoundingBox(){
 	std::cout << "Bounding Box: " << minElementArray[0] << " <= x <= " << maxElementArray[0] << ", " << minElementArray[1]
 	<< " <= y <= " << maxElementArray[1] << ", " << minElementArray[2] << " <= z <= " << maxElementArray[2] << std::endl;
 } // constructor for printing the bounding box


void PlyObject::findMeanArray(){
    meanArray = new double[n_axes-1];
    for(int i=0; i<n_axes-1; i++){meanArray[i] = 0;} // initialize array

        for (int i=0; i<(n_axes-1); i++)
        {
            for(int j=0; j<n_vertices; j++)
                {
                    meanArray[i] = meanArray[i] + vertices[i][j];
                }
            meanArray[i] = meanArray[i]/n_vertices;
        }//for
        std::cout << "Mean Vertex = (" << meanArray[0] << ", " << meanArray[1] << ", " << meanArray[2] << ")" << std::endl;
} // constructor to find mean of the coordinates

 void PlyObject::findStdDeviationArray(){
     stdDeviationArray = new double[n_axes-1];
     for(int i=0; i<n_axes-1; i++){stdDeviationArray[i] = 0;} // initialize array
////     double mean[n_axes-1] = {0};
////     double variance[n_axes-1] = {0};
     double mean[n_axes-1];
     double variance[n_axes-1];

         for (int i=0; i<(n_axes-1); i++)
         {
             for(int j=0; j<n_vertices; j++)
                 {
                     mean[i] = mean[i] + vertices[i][j];
                 }
             mean[i] = mean[i]/n_vertices;
             //std::cout << "mean is " << mean[i] << std::endl;

             for(int j=0; j<n_vertices; j++)
                 {
                     variance[i] = variance[i] + (vertices[i][j] - mean[i])*(vertices[i][j] - mean[i]);
                 }
             variance[i] = variance[i]/n_vertices;
             stdDeviationArray[i] = sqrt(variance[i]);
             //std::cout << "sd is " << stdDeviationArray[i] << std::endl;
         }//for
         std::cout << "Standard Deviations: x = " << stdDeviationArray[0] << ", y = " << stdDeviationArray[1] << ", z = " << stdDeviationArray[2] << std::endl;
 } // constructor to find standard deviation of the coordinates

 void PlyObject::objectCentering(){
     double rotationMatrix[4][4] = {
         {1.0,0,0,0},
         {0,1.0,0,0},
         {0,0,1.0,0},
         {0,0,0,1.0}};
////     double mean[n_axes-1] = {0}, sum[n_axes-1] = {0}, sd[n_axes-1] = {0};
     double mean[n_axes-1], sum[n_axes-1], sd[n_axes-1];
     meanArray = new double[n_axes-1];

     for (int i=0; i< (n_axes-1); i++)
         {
             for(int j=0; j<n_vertices; j++)
                 {
                     sum[i] = sum[i]+ vertices[i][j];
                 }
             mean[i] = sum[i]/n_vertices;
             meanArray[i] = mean[i];
             //std::cout << "mean is " << mean[i] << std::endl;

             for(int j=0; j<n_vertices; j++)
                 {
                     vertices[i][j] = vertices[i][j] - mean[i];
                 }//object centering
             mean[i] = 0;
             sum[i] = 0;
         }//for
 }// constructor for object centering

 void PlyObject::objectRounding(){
     double rotationMatrix[4][4] = {
         {1.0,0,0,0},
         {0,1.0,0,0},
         {0,0,1.0,0},
         {0,0,0,1.0}};
////     double mean[n_axes-1] = {0}, sum[n_axes-1] = {0}, sd[n_axes-1] = {0};
     double mean[n_axes-1], sum[n_axes-1], sd[n_axes-1];
     int i = 0, j = 0;
     stdDeviationArray = new double[n_axes-1];

     for (i=0; i< (n_axes-1); i++)
         {
             for(j=0; j<n_vertices; j++)
                 {
                     sum[i] = sum[i] + vertices[i][j];
                 }
             mean[i] = sum[i]/n_vertices;
             //std::cout << "mean is " << mean[i] << std::endl;

             for(j=0; j<n_vertices; j++)
                 {
                     sum[i] = sum[i] + (vertices[i][j] - mean[i])*(vertices[i][j] - mean[i]);
                 }
             sd[i] = sqrt(sum[i]/n_vertices);
             stdDeviationArray[i] = sd[i];
             //std::cout << "sd is " << sd[i] << std::endl;

             for(j=0; j<n_vertices; j++)
                 {
                     vertices[i][j] = vertices[i][j]/ sd[i];
                 }//object rounding
             mean[i] = 0;
             sum[i] = 0;
         }//for
 }// constructor for object rounding

void PlyObject::getPlyObjectVerticesFaceAttributes(std::string inputFileName){
    ptr_verticesMat = new Matrix2D<double>(this->n_vertices, this->n_axes, 1);
    //ptr_verticesMat->printMatrix();
    ptr_faceMat = new Matrix2D<unsigned>(this->n_faces, this->max_polygonSides, 0);
    //ptr_faceMat->printMatrix();

    std::string line;
    std::ifstream plyFile (inputFileName.c_str());

    if (plyFile.is_open())
    {
        do
            {
                getline(plyFile, line);
                //std::cout << line << '\n';
            }while (line != "end_header"); // while1
        //std::cout << "header ended --------------------. " << '\n';

        //--------------------- j X i : j rows, i columns ------------------------//
        int  j = 0;
        while (j!=this->n_vertices)
            {
                getline(plyFile, line);
                //std::cout << line << '\n';
                std::stringstream ss(line);
                double data[this->n_axes-1];

                for(int i=0; i<this->n_axes-1; i++)
                    {
                        ss >> data[i];
                        ptr_verticesMat->setElement(j, i, data[i]); // IMP indices are different: jth row, ith column *****
                        //this->assignMatrixValue(i, j, data[i], 1);
                        //std::cout << data[i] << " \t \t" << std::endl;
                    }//for
                j = j+1;
            }// while2
        //ptr_verticesMat->printMatrix();

        //--------------------- k X i : k rows, i columns ------------------------//
        int  k = 0;
        while (k!=this->n_faces)
            {
                getline(plyFile, line);
                //std::cout << line << '\n';
                std::stringstream ss(line);

                unsigned nVerts;
                ss >> nVerts; // read first column value
                //std::cout << nVerts << '\n';
                ptr_faceMat->setElement(k, 0, nVerts); // wrote the nVerts value in the 0th column
                unsigned data[nVerts]; // its number of vertices == unsigned

                for(int i=1; i<nVerts+1; i++)
                    {
                        ss >> data[i];
                        ptr_faceMat->setElement(k, i, data[i]); // IMP indices are different: kth row, ith column *****
                        //this->assignMatrixValue(i, j, data[i], 1);
                        //std::cout << data[i] << " \t \t" << std::endl;
                    }//for

                k = k+1;
            }// while3
        //ptr_faceMat->printMatrix();

        plyFile.close();
    }//if open

    else std::cout << "Unable to open file / File not found";
    //ptrPlyObject->printCoordinates();
    //std::cout << "------------------- after function call \n";
}

PlyObject::~PlyObject(){
	for(int i = 0; i < n_axes; i++)
        {
            delete [] vertices[i];
        }
	delete [] vertices;

}//Ply object destructor


#endif // PLYOBJECT_H_INCLUDED
