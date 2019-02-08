#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

#include <sstream>
#include <fstream>
#include <math.h>

#include <iostream>
#include <algorithm>

#include "vector3d.h"
#include "matrix2d.h"

#include <map>
#include <utility>

#define PI 3.14159265


class Model{
public:
    int n_axes, n_vertices, n_vns, n_faces, max_polygonSides;
	
	// smoothing cow
	int smoothingType = 1;  // 0: no smoothing; 1: smoothing a4 2018
	double smoothingAngleLimit = 22;  // angles between diff versions (coming from diff faces) of same vertex normal should be smaller than this angle
	double smoothingAngleCosineLimit = 0.92387953251;  // 0.86602; // 0.92718385456;  // cosine between two vn should be greater than this number
	
    double **vertices;
    double *maxElementArray;
    double *minElementArray;
    double *meanArray;
    double *stdDeviationArray;

    std::string modelName;
    std::string modelMaterialName;

    std::string firstCommentBlock;

    // a4 2017:  to store the materials given in the mtl file
    std::vector <std::string> materials;

    //Vector3D<double> rgb;    // rgb value vector r = x; g = y; b = z
    
    Vector3D<double> k_temp; // k temp
    // std::vector because it will store k coefficient for all the materials defined in the .mtl file
    std::vector <double> Ns; // specular exponent!
    std::vector <double> Ni; // 
    std::vector < Vector3D<double> > ka;    // ambient
    std::vector < Vector3D<double> > kd;    // diffuse
    std::vector < Vector3D<double> > ks;    // specular (phong)
    Vector3D<double> kr;    // attenuation multiplier!
    std::vector < Vector3D<double> > ke;    // emission coefficient
    std::vector <double> rho;  // 1 means only reflected light // 0 means only surface illumination

    // model transformation attributes
    Vector3D<double> preTransformVector;
    Vector3D<double> scaleVector;
    Vector3D<double> translationVector;
    Vector3D<double> rotationAxis;
    Vector3D<double> transformedVector;
    double rotationAngleDegrees;
    Vector3D<double> preTransformNormal;
    Vector3D<double> transformedNormal;

    Matrix2D<double>* ptr_rotationMat;
	Matrix2D<double>* ptr_rwMat;
	Matrix2D<double>* ptr_rwTMat;

    // for sphere fitting
    Vector3D<double> object_center;
    double radiusOfSphere;
    //int n_polygonSides = 4;

    Matrix2D<double>* ptr_verticesMat;
	//Matrix2D<double>* ptr_smoothNormalForVerticesMat;  // for every vertex store a smooth normal a4_2018
    Matrix2D<unsigned>* ptr_faceMat; // unsigned because it's either number of vertices per face or vertex number
    Matrix2D<double>* ptr_vnsMat;
    std::vector <unsigned> facelist_material_index;  // stores the 0-index of the material used for every face. Its length = n_faces
    //double **rotationMatrix;
    //double vertices[m][n];
    //bool isPLY = false, isFormatASCII = false;
	
	// a4 2018 multimap to store the list of faces for which the specific vertex was used
	std::multimap<unsigned, unsigned> vertex_to_face;
	//std::multimap<unsigned, unsigned>* ptr_vertex_to_face;
	std::pair<std::multimap<unsigned, unsigned>::iterator, std::multimap<unsigned, unsigned>::iterator> mmm;
	std::pair<std::multimap<unsigned, unsigned>::iterator, std::multimap<unsigned, unsigned>::iterator> mmm1;
	std::pair<std::multimap<unsigned, unsigned>::iterator, std::multimap<unsigned, unsigned>::iterator> mmm2;
	std::pair<std::multimap<unsigned, unsigned>::iterator, std::multimap<unsigned, unsigned>::iterator> mmm3;
	
	
    Model(int, int, int, int);

	void calculateRotationMatrix();
    void getNumberOfVerticesAndFaces(std::string);
    void getModelVerticesFaceAttributes(std::string);
    void getModelMaterialAttributes(std::string);
    void getModelCommentBlock(std::string);

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

    void printTransformedModelVerticesFaces();

    void printModelAttributes();

    ~Model(); //Model destructor

public:
    void allocateMatrix(){
		vertices = new double*[n_axes];
		for(int i = 0; i < n_axes; i++)
            {
                vertices[i] = new double[n_vertices];
            }
	}
}; // Model class

Model::Model(int nAxes, int nVertices, int nFaces, int max_polySides){
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


void Model::assignMatrixValue(int nAxes, int nVertices, double value, int flag){
	double newValue = value;
    vertices[nAxes][nVertices] = newValue;
    //std::cout << "assigned to vertex matrix! \n";
} // constructor to assign matrix value

void Model::printCoordinates(){
	for(int i = 0; i < n_axes; i++)
        {
            std::cout << "axis " << i << " vertices of the Model --->" << std::endl;
            std::cout << "" << std::endl;
            for(int j = 0; j < n_vertices; j++)
            {
                std::cout << vertices[i][j] << "  " << std::endl;
            }
            std::cout << "" << std::endl;
        }
} // constructor for printing the coordinates of the vertices

void Model::findMaxElementArray(){
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

void Model::findMinElementArray(){
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


void Model::getObjectCenter(){
    object_center.x = 0.5*(maxElementArray[0]+minElementArray[0]);
    object_center.y = 0.5*(maxElementArray[1]+minElementArray[1]);
    object_center.z = 0.5*(maxElementArray[2]+minElementArray[2]);
    std::cout << "center is --- " << object_center.x << ", " << object_center.y << ", " << object_center.z << std::endl;
}

void Model::getSphereRadius(){
    radiusOfSphere = std::max((maxElementArray[0]-object_center.x), (maxElementArray[1]-object_center.y));
    radiusOfSphere = std::max(radiusOfSphere,(maxElementArray[2]-object_center.z));
    std::cout << "max radius is --- " << radiusOfSphere << std::endl;
}

 void Model::printBoundingBox(){
 	std::cout << "Bounding Box: " << minElementArray[0] << " <= x <= " << maxElementArray[0] << ", " << minElementArray[1]
 	<< " <= y <= " << maxElementArray[1] << ", " << minElementArray[2] << " <= z <= " << maxElementArray[2] << std::endl;
 } // constructor for printing the bounding box


void Model::findMeanArray(){
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


 void Model::findStdDeviationArray(){
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

 void Model::objectCentering(){
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

 void Model::objectRounding(){
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


void Model::calculateRotationMatrix(){
	ptr_rotationMat = new Matrix2D<double>(3, 3, 0);
	ptr_rwMat = new Matrix2D<double>(3, 3, 0);
	ptr_rwTMat = new Matrix2D<double>(3, 3, 0);
	Vector3D<double> m, u, v;

	this->rotationAxis.normalize();
  std::cout << "normalized rotation axis:" << std::endl;
  this->rotationAxis.printVector3D("");

	// to make it random ->
	m.x = this->rotationAxis.x - this->rotationAxis.y - 5.235555;
	m.y = this->rotationAxis.y + this->rotationAxis.z - 1.1500000;
	m.z = this->rotationAxis.z - 2.5557777;
	m.normalize();
  std::cout << "normalized m (=w):" << std::endl;
  m.printVector3D("");

	u = this->rotationAxis.crossProduct(m);
  // u = m.crossProduct(this->rotationAxis);
  u.normalize();
  std::cout << "u vector:" << std::endl;
  u.printVector3D("");

	v = this->rotationAxis.crossProduct(u);  // w X u = v
  // v = u.crossProduct(this->rotationAxis);  // u X w = v
  v.normalize();
  std::cout << "v vector:" << std::endl;
  v.printVector3D("");

  // anti-clockwise rotation =>
  // | cos() -sin() |
  // | sin()  cos() |
	ptr_rotationMat->setElement(0, 0, cos(this->rotationAngleDegrees*PI/180));
	ptr_rotationMat->setElement(0, 1, -1*sin(this->rotationAngleDegrees*PI/180));
	ptr_rotationMat->setElement(0, 2, 0);

	ptr_rotationMat->setElement(1, 0, 1*sin(this->rotationAngleDegrees*PI/180));
	ptr_rotationMat->setElement(1, 1, cos(this->rotationAngleDegrees*PI/180));
	ptr_rotationMat->setElement(1, 2, 0);

	ptr_rotationMat->setElement(2, 0, 0);
	ptr_rotationMat->setElement(2, 1, 0);
	ptr_rotationMat->setElement(2, 2, 1);

  std::cout << "Rotation matrix is:" << std::endl;
	ptr_rotationMat->printMatrix();


	ptr_rwMat->setElement(0, 0, u.x);
	ptr_rwMat->setElement(0, 1, u.y);
	ptr_rwMat->setElement(0, 2, u.z);

	ptr_rwMat->setElement(1, 0, v.x);
	ptr_rwMat->setElement(1, 1, v.y);
	ptr_rwMat->setElement(1, 2, v.z);

	ptr_rwMat->setElement(2, 0, this->rotationAxis.x);
	ptr_rwMat->setElement(2, 1, this->rotationAxis.y);
	ptr_rwMat->setElement(2, 2, this->rotationAxis.z);

  std::cout << "rw matrix is:" << std::endl;
	ptr_rwMat->printMatrix();


	ptr_rwTMat->setElement(0, 0, u.x);
	ptr_rwTMat->setElement(1, 0, u.y);
	ptr_rwTMat->setElement(2, 0, u.z);

	ptr_rwTMat->setElement(0, 1, v.x);
	ptr_rwTMat->setElement(1, 1, v.y);
	ptr_rwTMat->setElement(2, 1, v.z);

	ptr_rwTMat->setElement(0, 2, this->rotationAxis.x);
	ptr_rwTMat->setElement(1, 2, this->rotationAxis.y);
	ptr_rwTMat->setElement(2, 2, this->rotationAxis.z);

  std::cout << "rwT matrix is:" << std::endl;
	ptr_rwTMat->printMatrix();

	(*ptr_rotationMat) = (*ptr_rotationMat)*(*ptr_rwMat);
	// ptr_rotationMat->printMatrix();
	(*ptr_rotationMat) = (*ptr_rwTMat)*(*ptr_rotationMat);

  std::cout << "The final rotation matrix is:" << std::endl;
	ptr_rotationMat->printMatrix();

	// testing ----------->>>>>
	Vector3D<double> A(1,0,0);
	Vector3D<double> B;

  std::cout << "Before rotation:" << std::endl;
  A.printVector3D("");
	B = (*ptr_rotationMat)*A;
  std::cout << "After rotation:" << std::endl;
	B.printVector3D("");

}


void Model::getNumberOfVerticesAndFaces(std::string inputFileName){
    std::string line;
    std::string searchKeyword_v("v "); // IMPORTANT: SPACE is necessary after the keyword
    std::string searchKeyword_vn ("vn "); // IMPORTANT: SPACE is necessary after the keyword
    std::string searchKeyword_f ("f "); // IMPORTANT: SPACE is necessary after the keyword
    std::ifstream sceneFile(inputFileName.c_str());

    this->n_vertices = 0;
    this->n_vns = 0;
    this->n_faces = 0;

    while (1)
        {
            getline(sceneFile, line);

            if (line.compare(0, searchKeyword_v.length(), searchKeyword_v) == 0){
                //std::cout << line << " light is present in the line " << std::endl;
                this->n_vertices++;
            }

            if (line.compare(0, searchKeyword_vn.length(), searchKeyword_vn) == 0){
                //std::cout << line << " sphere is present in the line " << std::endl;
                this->n_vns++;
            }

            if (line.compare(0, searchKeyword_f.length(), searchKeyword_f) == 0){
                //std::cout << line << " model is present in the line " << std::endl;
                this->n_faces++;
            }

            //std::cin.get();
            if(sceneFile.eof()){break;} // break if end of file
        } // while1
    sceneFile.close();

}


void Model::getModelVerticesFaceAttributes(std::string inputFileName){

    // CHANGE this later
    this->n_axes = 3;
    std::cout << this->n_axes << " n_axes (x, y, z) (test)" << std::endl;
    ptr_verticesMat = new Matrix2D<double>(this->n_vertices, this->n_axes, 1);
	//ptr_smoothNormalForVerticesMat = new Matrix2D<double>(this->n_vertices, this->n_axes, 1);
    ptr_vnsMat = new Matrix2D<double>(this->n_vns, this->n_axes, 1);
    //ptr_vnsMat->printMatrix();
    ptr_faceMat = new Matrix2D<unsigned>(this->n_faces, 9, 0); // 9 because 1/2/3 4/5/6 7/8/9
    //ptr_faceMat->printMatrix();

    std::string line;
    std::string term[4];
    std::string term1, term2, term3, term4; // to store the terms v, vn, or f
    std::string slash ("/");
    // std::cout << slash << " <------- slash" << std::endl;
	  calculateRotationMatrix();
    // // => A5_2016
    std::string searchKeyword_v("v "); // IMPORTANT: SPACE is necessary after the keyword
    std::string searchKeyword_vn ("vn "); // IMPORTANT: SPACE is necessary after the keyword
    std::string searchKeyword_f ("f "); // IMPORTANT: SPACE is necessary after the keyword
    std::string searchKeyword_mtl ("mtllib "); // IMPORTANT: SPACE is necessary after the keyword
    std::string searchKeyword_usemtl("usemtl "); // IMPORTANT: SPACE is necessary after the keyword
    std::ifstream sceneFile(inputFileName.c_str());

    //std::cout << " test 2" << std::endl;
    n_vertices = 0;
    n_vns = 0;
    n_faces = 0;

    unsigned current_mtl_index;  // 1 indexed instead of zero indexed

    while (1)
        {
            getline(sceneFile, line);
            std::stringstream ss(line);
            double data1[this->n_axes-1];  // for vertices
            double data2[this->n_axes-1];  // for vn
            double data3[this->n_axes-1];  // not used

            if (line.compare(0, searchKeyword_mtl.length(), searchKeyword_mtl) == 0){
                ss >> term1 >> this->modelMaterialName;
                std::cout << this->modelMaterialName << " this->modelMaterialName" << std::endl;
            }

            if (line.compare(0, searchKeyword_usemtl.length(), searchKeyword_usemtl) == 0){
                current_mtl_index = current_mtl_index + 1;  // 1 indexed!
            }

            if (line.compare(0, searchKeyword_v.length(), searchKeyword_v) == 0){
                //std::cout << " test v" << std::endl;
                //std::cout << line << " light is present in the line " << std::endl;

                ss >> term1 >> data1[0] >> data1[1] >> data1[2];
                //std::cout << term1 << "\t" << data1[0] << " \t" << data1[1] << " \t" << data1[2] << " \t \t" << std::endl;

                this->preTransformVector.x = data1[0];
                this->preTransformVector.y = data1[1];
                this->preTransformVector.z = data1[2];

                // ORDER is Rotate, Scale, and then Translate (A1_2017)
                // ANGLE AXIS ROTATION =======================================================================
                // std::cout << "Rotating data ..." << std::endl;
                this->transformedVector = (*ptr_rotationMat)*(this->preTransformVector);
                // ===========================================================================================

                // SCALING ===================================================================================
                // std::cout << "Scaling data ..." << std::endl;
                this->transformedVector.x = this->scaleVector.x*this->transformedVector.x;
                this->transformedVector.y = this->scaleVector.y*this->transformedVector.y;
                this->transformedVector.z = this->scaleVector.z*this->transformedVector.z;
                // ===========================================================================================

                // TRANSLATION ===============================================================================
                // std::cout << "Translating data ..." << std::endl;
                this->transformedVector.x = this->transformedVector.x + this->translationVector.x;
                this->transformedVector.y = this->transformedVector.y + this->translationVector.y;
                this->transformedVector.z = this->transformedVector.z + this->translationVector.z;
                // ===========================================================================================

                ptr_verticesMat->setElement(this->n_vertices, 0, this->transformedVector.x);
                ptr_verticesMat->setElement(this->n_vertices, 1, this->transformedVector.y);
                ptr_verticesMat->setElement(this->n_vertices, 2, this->transformedVector.z);

                this->n_vertices++;
            } // if for vs

            if (line.compare(0, searchKeyword_vn.length(), searchKeyword_vn) == 0){
                //std::cout << " test vn" << std::endl;
                //std::cout << line << " sphere is present in the line " << std::endl;
                ss >> term1 >> data2[0] >> data2[1] >> data2[2];
                //std::cout << term1 << " \t" << data2[0] << " \t" << data2[1] << " \t" << data2[2] << " \t \t" << std::endl;

                // ptr_vnsMat->setElement(this->n_vns, 0, data2[0]);
                // ptr_vnsMat->setElement(this->n_vns, 1, data2[1]);
                // ptr_vnsMat->setElement(this->n_vns, 2, data2[2]);

                // TODO vn rotation
                this->preTransformNormal.x = data2[0];
                this->preTransformNormal.y = data2[1];
                this->preTransformNormal.z = data2[2];
                std::cout << "this->preTransformNormal: " << this->preTransformNormal.x << ", " << this->preTransformNormal.y << ", " << this->preTransformNormal.z << " " << std::endl;
                this->transformedNormal = (*ptr_rotationMat)*(this->preTransformNormal);
                std::cout << "this->transformedNormal: " << this->transformedNormal.x << ", " << this->transformedNormal.y << ", " << this->transformedNormal.z << " " << std::endl;

                ptr_vnsMat->setElement(this->n_vns, 0, this->transformedNormal.x);
                ptr_vnsMat->setElement(this->n_vns, 1, this->transformedNormal.y);
                ptr_vnsMat->setElement(this->n_vns, 2, this->transformedNormal.z);
                // Done

                this->n_vns++;
            } // if for vns

            if (line.compare(0, searchKeyword_f.length(), searchKeyword_f) == 0){
                //std::cout << " test f" << std::endl;
                std::cout << line << " f is present in the line " << std::endl;
                ss >> term[0] >> term[1] >> term[2] >> term[3];
                //std::cout << term[0] << " \t" << term[1] << " \t" << term[2] << " \t" << term[3] << std::endl;

                // assign the material index that was read in inside current_mtl_index (1 indexed)
                facelist_material_index.push_back(current_mtl_index - 1);  //-1 because it is 1 indexed 

                for (unsigned t=1; t<4; t++){
                    std::size_t found1 = term[t].find(slash); // first occurence
                    std::size_t found2 = term[t].rfind(slash); // last occurence
                    //std::cout << found << " found" << std::endl;

                    if (found1!=std::string::npos && found2!=std::string::npos){
                        //std::cout << found1 << " found1 in the term" << std::endl;
                        //std::cout << found2 << " found2 in the term" << std::endl;

                        std::string num1 = term[t].substr(0, found1-0);
                        std::string num3 = term[t].substr(found2+1); // till the end
                        std::string num2;

                        if(((int)found2-(int)found1) == 1){
                            //std::cout << (int)found2-(int)found1 << " found2 - found1" << std::endl;
                            num2 = "0";
                        }

                        else {
                            num2 = term[t].substr(found1+1, found2-found1-1);
                        }

						std::cout << this->n_faces << " this->n_faces" << std::endl;  // starts with 0 (0-indexed) and face mat is also 0-indexed
                        std::cout << num1 << " <- double num1" << std::endl;
                        std::cout << num3 << " <- double num3" << std::endl;
                        std::cout << num2 << " <- double num2" << std::endl;

						unsigned number1 = (unsigned)atof(num1.c_str());
						unsigned number2 = (unsigned)atof(num2.c_str());
						unsigned number3 = (unsigned)atof(num3.c_str());
							
						// In the multimap, create a key with vertex index (number1) and value is the current index of face (this->n_faces)
						// IMPORTANT !!!!!!!!!!!!!! key (vertex index) is 1-indexed, values (face index) are 0-indexed
						vertex_to_face.insert(std::pair<unsigned, unsigned>(number1, (unsigned)this->n_faces));
						
						
						//std::cout << atof(num.c_str()) << " atof" << std::endl;
                        // for t=1; coulumns are: 0,1,2; for t=2; coulumns are: 3,4,5; for t=3; coulumns are: 6,7,8;
                        ptr_faceMat->setElement(this->n_faces, 3*t - 3, number1);  //0 for t=1
                        ptr_faceMat->setElement(this->n_faces, 3*t - 2, number2);  //1
                        ptr_faceMat->setElement(this->n_faces, 3*t - 1, number3);  //2
                    }

                } // for each term 1, 2, 3

                this->n_faces++;
            } // if for faces

            //std::cin.get();
            if(sceneFile.eof()){break;} // break if end of file
        } // while1
    sceneFile.close();
}


void Model::getModelMaterialAttributes(std::string inputFileName){

    std::cout << "reading model material properties for " << inputFileName << std::endl;
	// CHANGE THIS LATER  ...changed.
	this->kr.x = 1.0;  // 0.9;
	this->kr.y = 1.0;  // 0.9;
	this->kr.z = 1.0;  // 0.9;

    // this->rho = 0.9; 
    double rhoterm = 0.9;

    std::string line;
    std::string term1; // to store the terms Ka, Kd, Ks;
    double nterm;
    std::string mtlname;
    // int mtlcount = -1;  // to make indexing start from 0

    std::string searchKeyword_newmtl("newmtl "); // IMPORTANT: SPACE is necessary after the keyword
    std::string searchKeyword_Ns("Ns "); // IMPORTANT: SPACE is necessary after the keyword
	std::string searchKeyword_Ni("Ni "); // IMPORTANT: SPACE is necessary after the keyword
    std::string searchKeyword_Ka("Ka "); // IMPORTANT: SPACE is necessary after the keyword
    std::string searchKeyword_Kd ("Kd "); // IMPORTANT: SPACE is necessary after the keyword
    std::string searchKeyword_Ks ("Ks "); // IMPORTANT: SPACE is necessary after the keyword
    std::string searchKeyword_Ke ("Ke "); // IMPORTANT: SPACE is necessary after the keyword
    std::ifstream sceneFile(inputFileName.c_str());

    while (1)
        {
            getline(sceneFile, line);
            std::stringstream ss(line);

            if (line.compare(0, searchKeyword_newmtl.length(), searchKeyword_newmtl) == 0){
                ss >> term1 >> mtlname;
                this->materials.push_back(mtlname);
                // mtlcount = mtlcount + 1;
            }
                  
			if (line.compare(0, searchKeyword_Ns.length(), searchKeyword_Ns) == 0){
                ss >> term1 >> nterm;
                this->Ns.push_back(nterm);
            }

            if (line.compare(0, searchKeyword_Ka.length(), searchKeyword_Ka) == 0){
                ss >> term1 >> this->k_temp.x >> this->k_temp.y >> this->k_temp.z;
                this->ka.push_back(k_temp);
            }

			if (line.compare(0, searchKeyword_Kd.length(), searchKeyword_Kd) == 0){
                ss >> term1 >> this->k_temp.x >> this->k_temp.y >> this->k_temp.z;
                this->kd.push_back(k_temp);
            }

            if (line.compare(0, searchKeyword_Ks.length(), searchKeyword_Ks) == 0){
                ss >> term1 >> this->k_temp.x >> this->k_temp.y >> this->k_temp.z;
                this->ks.push_back(k_temp);
            }

            if (line.compare(0, searchKeyword_Ke.length(), searchKeyword_Ke) == 0){
                ss >> term1 >> this->k_temp.x >> this->k_temp.y >> this->k_temp.z;
                this->ke.push_back(k_temp);
            }

            if (line.compare(0, searchKeyword_Ni.length(), searchKeyword_Ni) == 0){
                ss >> term1 >> nterm;
                this->Ni.push_back(nterm);
            }

            this->rho.push_back(rhoterm);

            //std::cin.get();
            if(sceneFile.eof()){break;} // break if end of file
        } // while1
    sceneFile.close();
}


void Model::getModelCommentBlock(std::string inputFileName){
    // stores the intial comment block in a string
    std::string line;
    std::string searchKeyword_hash("#"); // searches for # representing a comment
    std::ifstream sceneFile(inputFileName.c_str());

    this->firstCommentBlock = "";
    while (1)
        {
            getline(sceneFile, line);

            if (line.compare(0, searchKeyword_hash.length(), searchKeyword_hash) == 0){
                //std::cout << line << " light is present in the line " << std::endl;
                this->firstCommentBlock = this->firstCommentBlock + line + "\n";
            }

            else {
              break;  // break immidiately after no # is found => first comment block ended
            }

            //std::cin.get();
            if(sceneFile.eof()){break;} // break if end of file

        } // while1
    sceneFile.close();

}


void Model::printTransformedModelVerticesFaces(){
  std::cout << "transformed Model Vertices -------------------------------------------------->" << std::endl;
  ptr_verticesMat->printMatrix();
  ptr_faceMat->printMatrix();

}


//void Model::printModelVertexToFacesMultimap(){	
//}


void Model::printModelAttributes(){
    std::cout << "Model Attributes -------------------------------------------------->" << std::endl;
    std::cout << modelName << " <----- model name " << modelMaterialName << " <----- model material name" << std::endl;
    std::cout << n_vertices << " <- nVertices " << n_faces << " <- nFaces " << n_vns << " <- nVns" << std::endl;
    std::cout << "n(rows, columns) vertMat = " << ptr_verticesMat->getN_Rows() << " " << ptr_verticesMat->getN_Columns() << " last element = " << std::endl; // << ptr_verticesMat->getElement(7, 2) << std::endl;
    std::cout << "n(rows, columns) faceMat = " << ptr_faceMat->getN_Rows() << " " << ptr_faceMat->getN_Columns() << " last element = " << std::endl; // << ptr_faceMat->getElement(11, 0) << std::endl;
    std::cout << "n(rows, columns) normMat = " << ptr_vnsMat->getN_Rows() << " " << ptr_vnsMat->getN_Columns() << " last element = " << std::endl; // << ptr_vnsMat->getElement(5, 1) << std::endl;

    std::cout << "vert mat ------------------------------------" << std::endl;
    ptr_verticesMat->printMatrix();
    std::cout << "face mat ------------------------------------" << std::endl;
    ptr_faceMat->printMatrix();
    std::cout << "vn mat ------------------------------------" << std::endl;
    ptr_vnsMat->printMatrix();
    std::cout << this->facelist_material_index.size() << " size(facelist_material_index) ------------------------------------" << std::endl;
    for (int m=0; m<this->facelist_material_index.size(); m++){
        std::cout << this->facelist_material_index[m] << " " << std::endl;
    }

	std::cout << "Model Vertex to Faces Multimap -------------------------------------------------->" << std::endl;
	for (unsigned v=0; v<this->n_vertices+1; v++){
		std::cout << "vertex = " << v << "  Number of elements with this key v =" << vertex_to_face.count(v) << std::endl;
		mmm = vertex_to_face.equal_range(v);
		
		for (std::multimap<unsigned, unsigned>::iterator it2 = mmm.first; it2 != mmm.second; ++it2){
			std::cout << "  [" << (*it2).first << ", " << (*it2).second << "]" << std::endl;
		}
		std::cout << "-----------------------------------------------------------------------\n" << std::endl;
		
	}
	 
	
	
    std::cout << scaleVector.x << " " << scaleVector.y << " " << scaleVector.z << " <----- scaleVector" << std::endl;
    std::cout << translationVector.x << " " << translationVector.y << " " << translationVector.z << " <----- translationVector" << std::endl;
    std::cout << rotationAxis.x << " " << rotationAxis.y << " " << rotationAxis.z << " <----- rotationAxis" << std::endl;
    std::cout << rotationAngleDegrees << " <----- rotationAngleDegrees" << std::endl;
    
    for (int i=0; i<this->materials.size(); i++)
    {  
        std::cout << materials[i] << "---------------------------------------------" << std::endl;   
        std::cout << Ns[i] << " <- Ns specular exponent " << std::endl;
        std::cout << ka[i].x << " " << ka[i].y << " " << ka[i].z << " <----- Ambient" << std::endl;
        std::cout << kd[i].x << " " << kd[i].y << " " << kd[i].z << " <----- Diffuse" << std::endl;
        std::cout << ks[i].x << " " << ks[i].y << " " << ks[i].z << " <----- Specular" << std::endl;
        std::cout << kr.x << " " << kr.y << " " << kr.z << " <----- Attenuation Multiplier" << std::endl;
        std::cout << ke[i].x << " " << ke[i].y << " " << ke[i].z << " <----- Emission Coefficients" << std::endl;
        std::cout << Ni[i] << " <- Ni " << std::endl;
        std::cout << rho[i] << " <- rho " << std::endl;
    }
    
    std::cout << "End ----------------------------------------------------------------" << '\n' << std::endl;

}

Model::~Model(){
	for(int i = 0; i < n_axes; i++)
        {
            delete [] vertices[i];
        }
	delete [] vertices;

} // object destructor

#endif // MODEL_H_INCLUDED
