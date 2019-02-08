#ifndef RAY_H_INCLUDED
#define RAY_H_INCLUDED

#include <iostream>
#include <limits>   // to get infinity!

#include "vector3d.h"
#include "sphere.h"
#include "model.h"

#include <map>
#include <cassert>

class Ray{
public:

    Vector3D<double> eye; // s  // A2_2017:  // same as the l in matrices in Cramers method
    Vector3D<double> d;  // s  // A2_2017: unit ray though the pixel?  // same as d from M matrices in Cramers method
    Vector3D<double> pixel_xyz; // A2_2017: xyz vector representing the xyz coordinates of the pixel from which this ray passes

    unsigned objectType;

    // A2_2017: distance from eye to intersection?
    double bestT_sphere = INFINITY;  // s
    double bestT_model = INFINITY;
    int bestFaceIndex; // = -1;
	Vector3D<double> bestNormal;

    Sphere* ptr_bestSph = NULL;  // s
    Model* ptr_bestModel = NULL;
    Vector3D<double> best_intersection_sphere; // by default a vector is initialized to zero  // s
    Vector3D<double> best_intersection_model; // by default a vector is initialized to zero
    //best_intersection_sphere.x = best_intersection_sphere.y = best_intersection_sphere.z = 0;

    // to store the 0-index of the material of the best face
    // unsigned bestFaceMaterialIndex; // cannot be retrieved since it need access to facelist_material_index from model.h

    Ray();
    Ray(double, double, double, double, double, double);

    void calculateDeterminantMAndM3(const Vector3D<double> &a, const Vector3D<double> &b, const Vector3D<double> &c, const Vector3D<double> &d, const Vector3D<double> &l, double &detM, double &detM3);
    void calculateDeterminantM1AndM2(const Vector3D<double> &a, const Vector3D<double> &b, const Vector3D<double> &c, const Vector3D<double> &d, const Vector3D<double> &l, double &detM1, double &detM2);
    void initializeArrayElements(double arr[], int arrSize, double initValue);
    double getMinPositiveElement(const double arr[], const int arrSize, int& index);

    int raySphereIntersectionTest(Sphere* ptr_sph);
    int rayModelIntersectionTest(Model* ptr_model, bool debugFlag);

    // return 1 => sphere; return 2 => .obj object; return 0 => no intersection
    unsigned findBestIntersectionObject();

    void printRayAttributes();
    ~Ray();
};

Ray::Ray(){
    eye.x=0, eye.y=0, eye.z=0, d.x=0, d.y=0, d.z=0;
}

Ray::Ray(double eyeX, double eyeY, double eyeZ, double dX, double dY, double dZ){
    eye.x = eyeX;
    eye.y = eyeY;
    eye.z = eyeZ;
    d.x = dX;
    d.y = dY;
    d.z = dZ;
}

int Ray::raySphereIntersectionTest(Sphere* ptr_sph){

    double scalarC, scalarV, scalarD, scalarTEye;
    scalarC = (ptr_sph->center - this->eye).getMagnitude();
    //std::cout << " inside sphere test ........................ " << std::endl;
    //std::cout << scalarC << " scalarC" << std::endl;
    //ptr_sph->center.printVector3D();
    //this->eye.printVector3D();

    //this->d.printVector3D();
    scalarV = (ptr_sph->center - this->eye).dotProduct(this->d);
    //std::cout << scalarV << " scalarV" << std::endl;
    //std::cin.get();
    //******************************************** check if center is infront of the image plane ******************************************
    //std::cout << ptr_sph->radius << " ptr_sph->radius" << std::endl;
    scalarD = (ptr_sph->radius)*(ptr_sph->radius) - (scalarC*scalarC-scalarV*scalarV); // d^2
	scalarTEye = (scalarV)-(sqrt(scalarD)); // = (v-d) // Q-E = (scalarV-scalarD)*d

    if(scalarD < 0 || scalarTEye < 0.000001){
        /////std::cout << scalarD << " scalarD is less than 0" << std::endl;
        //this->ptr_bestSph = NULL;
        return 0;
    }

    else {
        /////std::cout << scalarD << " scalarD is > 0" << std::endl;
        // distance from intersection to eye

        // Following change removed ghost reflection that were appearing at the corners!!!!!!!!!!!!!!!!!!!!!!
        //~ scalarTEye = fabs(scalarV)-fabs(sqrt(scalarD)); // = (v-d) // Q-E = (scalarV-scalarD)*d
        //scalarTEye = (scalarV)-(sqrt(scalarD)); // = (v-d) // Q-E = (scalarV-scalarD)*d

        // Update best values for the ray
        if (scalarTEye < this->bestT_sphere && scalarTEye > 0.000001){
            this->bestT_sphere = scalarTEye;
            this->ptr_bestSph = ptr_sph;
            this->best_intersection_sphere = this->eye + (this->d)*scalarTEye;
            ///std::cout << "Ray's best intersection ------------->" << std::endl;
            ///this->best_intersection_sphere.printVector3D();
        } //if
        return 1;
    } // else

} // raySphereIntersectionTest


inline void Ray::calculateDeterminantMAndM3(const Vector3D<double> &a, const Vector3D<double> &b, const Vector3D<double> &c, const Vector3D<double> &d, const Vector3D<double> &l, double &detM, double &detM3){
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

inline void Ray::calculateDeterminantM1AndM2(const Vector3D<double> &a, const Vector3D<double> &b, const Vector3D<double> &c, const Vector3D<double> &d, const Vector3D<double> &l, double &detM1, double &detM2){
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

void Ray::initializeArrayElements(double arr[], int arrSize, double initValue){
    //std::cout << "----------- printing the array of t for a single pixel -----------" << std::endl;
    for(int i=0; i<arrSize; i++){
        arr[i] = initValue;
    }
}

//void printArrayElements(double arr[], int arrSize){
//    std::cout << "----------- printing the array of t for a single pixel -----------" << std::endl;
//    for(int i=0; i<arrSize; i++){
//        std::cout << arr[i] << ", ";
//    }
//    std::cout << '\n' << std::endl;
//}

inline double Ray::getMinPositiveElement(const double arr[], const int arrSize, int& index){
    double tMin = INFINITY; // wrong-> needs to be a negative number!!!!! because for negative t, we write 239, 239, 239 in the final step
    for(int i=0; i<arrSize; i++){
        //if(tMin>arr[i] && arr[i]>0){tMin = arr[i];} // find minimum positive array element
        if(arr[i]>0.00001 && arr[i]<tMin){
			tMin = arr[i];
			index = i;
		}
    }
    return tMin;
}

int Ray::rayModelIntersectionTest(Model* ptr_model, bool debugFlag){
    // For every face::::: Find beta, gamma, and tDepth values
    int face = 0;
    unsigned v1, v2, v3;
    unsigned v1Min, v2Min, v3Min;
    unsigned n1Min, n2Min, n3Min;

    Vector3D<double> vn_1; // store vn for each pixel
    Vector3D<double> vn_2; // store vn for each pixel
    Vector3D<double> vn_3; // store vn for each pixel

    Vector3D<double> v_a; // store v for that face
    Vector3D<double> v_b; // store v for that face
    Vector3D<double> v_c; // store v for that face
    Vector3D<double> v_E1; // for edge 1
    Vector3D<double> v_E2; // for edge 2

    Vector3D<double> vn_true; // true vn for each vertex (actually it is using cross product method, so it is for the face, but the same vn is then applied to vertices); added for a4 2018
	Vector3D<double> vNA; // vn after smoothing
	Vector3D<double> vNB; // vn after smoothing
	Vector3D<double> vNC; // vn after smoothing
    Vector3D<double> vn_current; // store vn current
	unsigned v1MinCurrent, v2MinCurrent, v3MinCurrent;


	
    Vector3D<double> vn_final; // store vn for each pixel

    Vector3D<double> a, b, c; // store temporary a, b, c, d

    double detM=0, detM1=0, detM2=0, detM3=0;

    double eps = 0.0000001;
    //std::cout << eps << " <<<----------------------------- epsilon value for [zero] comparison " << std::endl;

    double beta, gamma, tDepth;

    double tOfAPixel[ptr_model->n_faces]; //stores the t values for a single pixel after scanning all the faces; finally we pick minimum t>0
    double betaOfAPixel[ptr_model->n_faces]; //stores the beta values for a single pixel after scanning all the faces; finally we pick the one with index of best t (minimum t>0)
    double gammaOfAPixel[ptr_model->n_faces]; //stores the gamma values for a single pixel after scanning all the faces; finally we pick the one with index of best t (minimum t>0)
	double bestBeta, bestGamma; 
		
    initializeArrayElements(tOfAPixel, ptr_model->n_faces, -1.1-10); // set it to negative and not zero since object could be touching the image plane

    double tMin = -2.1111;
    int index = -1;  // to store the index where minimum positive t occurs for that pixel


    for(face=0; face<ptr_model->n_faces; face++){
        //std::cout << " >>>>>>>>>>>>>>>>>>> for face number " << face << " <<<<<<<<<<<<<<<<< " << std::endl;
        ptr_model->ptr_faceMat->getElementVertices(face, v1, v2, v3); // this fills in the v1 v2 and v3 with 1 indexing
        //std::cout << v1 << " " << v2 << " " << v3 << " <---------------------------------- printing vi values " << std::endl;

        // following function takes care of 1 indexing in .obj format
        ptr_model->ptr_verticesMat->getElementXYZTriangleA4(v1, v2, v3, a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z);
        //std::cout << " ----------------- printing a vector after setting --------------------------- " << std::endl;
        //a.printVector3D();

//////////                ptr_model->ptr_verticesMat->getElementXYZ(v2,0, b.x, b.y, b.z);
//                ptr_model->ptr_verticesMat->getElement1(v2,1, b.y);
//                ptr_model->ptr_verticesMat->getElement1(v2,2, b.z);
////                std::cout << " ----------------- printing b vector after setting --------------------------- " << std::endl;
////                b.printVector3D();

//////////                ptr_model->ptr_verticesMat->getElementXYZ(v3,0, c.x, c.y, c.z);
//                ptr_model->ptr_verticesMat->getElement1(v3,1, c.y);
//                ptr_model->ptr_verticesMat->getElement1(v3,2, c.z);
////                std::cout << " ----------------- printing c vector after setting --------------------------- " << std::endl;
////                c.printVector3D();

////                std::cout << " <----------------- printing determinant values -------------------------- " << std::endl;
                calculateDeterminantMAndM3(a, b, c, this->d, this->eye, detM, detM3);
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

                calculateDeterminantM1AndM2(a, b, c, this->d, this->eye, detM1, detM2);
                //detM1 = calculateDeterminantM(a, pixeL, c, d);                //NOTE: b is replaced by pixeL
                beta = detM1/detM;
                //beta = (float)beta;  // 2017 edge pixel issue?
////                std::cout << detM1/detM << " <------------ beta = detM1/detM " << std::endl;
                if(beta<0){
                    tOfAPixel[face] = -3; // Here, tDepth is positive already. But beta is negative => No intersection! So we gotta neglect this positive tDepth by setting it to -3.
////                    std::cout << " set t value to -3 (beta negative) and continue to  the next face! " << std::endl;
                    continue;
                } // if beta is negative

                //detM2 = calculateDeterminantM(a, b, pixeL, d);                //NOTE: c is replaced by pixeL
                gamma = detM2/detM;
                //gamma = (float)gamma;
////                std::cout << detM2/detM << " <------------ gamma = detM2/detM " << std::endl;
                if(gamma<0){
                    tOfAPixel[face] = -5; // Here, tDepth is positive already. But gamma is negative => No intersection! So we gotta neglect this positive tDepth by setting it to -5.
////                    std::cout << " set t value to -5 (gamma negative) and continue to  the next face! " << std::endl;
                    continue;
                } // if gamma is negative

                //Final step of sieving (EARLY EXIT/RETURN strategy), check if (beta+gamma) <= 1
                //////////////////////// CHANGE this for square shaped faces!!!!!!!!!!
                if((beta+gamma)>1.00001){
                    tOfAPixel[face] = -7; // Here, tDepth is positive already. But (beta+gamma) > 1) => No intersection! So we gotta neglect this positive tDepth by setting it to -7.
////                    std::cout << " set t value to -7 (gamma negative) and continue to  the next face! " << std::endl;
                    continue;
                } // if (beta+gamma) > 1

                // When everything from pdf notes (ray-triangle intersection) is satisfied then assign tDepth which is positive to tOfAPixel[face] array element
                if(tDepth>0){
                    tOfAPixel[face] = tDepth;
                    betaOfAPixel[face] = beta;
                    gammaOfAPixel[face] = gamma;
                    //std::cout << tDepth << " this is final step tDepth value; it has to be positive! " << std::endl;
                }

        } // for loop face
	
	
            // also fills up index (pass by value); face index is 0 indexed
            tMin = getMinPositiveElement(tOfAPixel, ptr_model->n_faces, index);  // tMin for that pixel?
	
	if (debugFlag){
		std::cout << tMin << " tMin debug mode " << std::endl;
	}
	
////            std::cout << tMin << " <------------ tMin for pixel number -> " << row << ", " << column << std::endl;
	
            /////ptr_cm->ptr_depthMat->setElement(row, column, tMin);

            if (tMin < this->bestT_model && tMin > 0.000001){
				if (index != -1) {
					//~ std::cout << index << " <------- index " << tMin << " <------------ tMin /////END " << '\n' << std::endl; // << row << ", " << column << std::endl;
					this->bestFaceIndex = index;
                    // this->bestFaceMaterialIndex =  // cannot be retrieved since it need access to facelist_material_index from model.h

          // "index" is zero indexed
					// ptr_model->ptr_faceMat->getElementVertices(index, v1Min, v2Min, v3Min); // the v1Mins returned here are 1-indexed //old: note it's just index and not index+1
					// ptr_model->ptr_faceMat->getElementNormals(index, n1Min, n2Min, n3Min); // removed 2017 // returns 1 indexed normals //old: note it's just index and not index+1
          //~ std::cout << v1Min << " " << v2Min << " " << v3Min << std::endl;
					//~ std::cout << n1Min << " " << n2Min << " " << n3Min << std::endl;

					// from this index get bestBeta and bestGamma
					bestBeta = betaOfAPixel[index];
					bestGamma = gammaOfAPixel[index];
					
					// added 2017  // "index" is zero indexed
					ptr_model->ptr_faceMat->getElementVertices(index, v1Min, v2Min, v3Min); // the v1Mins returned here are 1-indexed //old: note it's just index and not index+1
					// takes care of 1-indexing in v1Min
					ptr_model->ptr_verticesMat->getElementXYZTriangleA4(v1Min, v2Min, v3Min, v_a.x, v_a.y, v_a.z, v_b.x, v_b.y, v_b.z, v_c.x, v_c.y, v_c.z);

					v_E1 = v_a-v_b;
					//v_E2 = v_c-v_a;
					v_E2 = v_a-v_c;
          
					
					if (ptr_model->smoothingType == 0){  // no smoothing
						vn_final = v_E1.crossProduct(v_E2);
           				vn_final.normalize();
					}
					
					// Now, "index" (which is 0-indexed) is the index of face in face mat. The v1Mins (which are 1-indexed) are the vertex forming that face. 
					// Use these vMins as keys to get values (face index that are 0-indexed) from vertex_to_face matrix. For each v1Min iterate over all faces and find running average of vns of those faces if they are close to true normal
					// In vertex_to_face: key (vertex index) is 1-indexed, values (face index) are 0-indexed
					
					else if (ptr_model->smoothingType == 1){  // smoothing a4 2018
						vn_true = v_E1.crossProduct(v_E2);
           				vn_true.normalize();
						
						//vNA = vn_true;
						//vNB = vn_true;
						//vNC = vn_true;
						
						// for vNA ============
						ptr_model->mmm1 = ptr_model->vertex_to_face.equal_range(v1Min);
						
						for (std::multimap<unsigned, unsigned>::iterator it2 = ptr_model->mmm1.first; it2 != ptr_model->mmm1.second; ++it2){
							//std::cout << "  [" << (*it2).first << ", " << (*it2).second << "]" << std::endl;
							// (*it2).second  // this is the index of face (0-indexed), get vn for that face, add it to running average only if its below 22 degree limit
							// since "index" was also 0-indexed on the faces, we can reuse (from above) the function below
							
							// added 2018 a4  // "(unsigned)(*it2).second" is zero indexed
							ptr_model->ptr_faceMat->getElementVertices((unsigned)(*it2).second, v1MinCurrent, v2MinCurrent, v3MinCurrent); // the v1Mins returned here are 1-indexed
           					// takes care of 1-indexing in v1Min
           					ptr_model->ptr_verticesMat->getElementXYZTriangleA4(v1MinCurrent, v2MinCurrent, v3MinCurrent, v_a.x, v_a.y, v_a.z, v_b.x, v_b.y, v_b.z, v_c.x, v_c.y, v_c.z);
							
							v_E1 = v_a-v_b;
							//v_E2 = v_c-v_a;
							v_E2 = v_a-v_c;
							
							vn_current = v_E1.crossProduct(v_E2);
           					vn_current.normalize();
							
							if (vn_current.dotProduct(vn_true) > ptr_model->smoothingAngleCosineLimit){   // then add to running average
								vNA = (vNA + vn_current)*0.5;  // average of two
								vNA.normalize();
							}
							
						}  // for
						// check if final vNX is close to v_true (within 22 degrees)
						assert(vNA.dotProduct(vn_true) > ptr_model->smoothingAngleCosineLimit);
						
						
						// for vNB ============
						ptr_model->mmm2 = ptr_model->vertex_to_face.equal_range(v2Min);
						
						for (std::multimap<unsigned, unsigned>::iterator it2 = ptr_model->mmm2.first; it2 != ptr_model->mmm2.second; ++it2){
							//std::cout << "  [" << (*it2).first << ", " << (*it2).second << "]" << std::endl;
							// (*it2).second  // this is the index of face (0-indexed), get vn for that face, add it to running average only if its below 22 degree limit
							// since "index" was also 0-indexed on the faces, we can reuse (from above) the function below
							
							// added 2018 a4  // "(unsigned)(*it2).second" is zero indexed
							ptr_model->ptr_faceMat->getElementVertices((unsigned)(*it2).second, v1MinCurrent, v2MinCurrent, v3MinCurrent); // the v1Mins returned here are 1-indexed
           					// takes care of 1-indexing in v1Min
           					ptr_model->ptr_verticesMat->getElementXYZTriangleA4(v1MinCurrent, v2MinCurrent, v3MinCurrent, v_a.x, v_a.y, v_a.z, v_b.x, v_b.y, v_b.z, v_c.x, v_c.y, v_c.z);
							
							v_E1 = v_a-v_b;
							//v_E2 = v_c-v_a;
							v_E2 = v_a-v_c;
							
							vn_current = v_E1.crossProduct(v_E2);
           					vn_current.normalize();
							
							if (vn_current.dotProduct(vn_true) > ptr_model->smoothingAngleCosineLimit){   // then add to running average
								vNB = (vNB + vn_current)*0.5;  // average of two
								vNB.normalize();
							}
							
							
						}  // for
						// check if final vNX is close to v_true (within 22 degrees)
						assert(vNB.dotProduct(vn_true) > ptr_model->smoothingAngleCosineLimit);
						
						
						// for vNC ============
						ptr_model->mmm3 = ptr_model->vertex_to_face.equal_range(v3Min);
						
						for (std::multimap<unsigned, unsigned>::iterator it2 = ptr_model->mmm3.first; it2 != ptr_model->mmm3.second; ++it2){
							//std::cout << "  [" << (*it2).first << ", " << (*it2).second << "]" << std::endl;
							// (*it2).second  // this is the index of face (0-indexed), get vn for that face, add it to running average only if its below 22 degree limit
							// since "index" was also 0-indexed on the faces, we can reuse (from above) the function below
							
							// added 2018 a4  // "(unsigned)(*it2).second" is zero indexed
							ptr_model->ptr_faceMat->getElementVertices((unsigned)(*it2).second, v1MinCurrent, v2MinCurrent, v3MinCurrent); // the v1Mins returned here are 1-indexed
           					// takes care of 1-indexing in v1Min
           					ptr_model->ptr_verticesMat->getElementXYZTriangleA4(v1MinCurrent, v2MinCurrent, v3MinCurrent, v_a.x, v_a.y, v_a.z, v_b.x, v_b.y, v_b.z, v_c.x, v_c.y, v_c.z);
							
							v_E1 = v_a-v_b;
							//v_E2 = v_c-v_a;
							v_E2 = v_a-v_c;
							
							vn_current = v_E1.crossProduct(v_E2);
           					vn_current.normalize();
							
							if (vn_current.dotProduct(vn_true) > ptr_model->smoothingAngleCosineLimit){   // then add to running average
								vNC = (vNC + vn_current)*0.5;  // average of two
								vNC.normalize();
							}
							
							
						}  // for
						// check if final vNX is close to v_true (within 22 degrees)
						assert(vNC.dotProduct(vn_true) > ptr_model->smoothingAngleCosineLimit);
						
						// final vn
						vn_final = (vNA*(1 - bestBeta - bestGamma)) + (vNB*bestBeta) + (vNC*bestGamma);
           				vn_final.normalize();	
								
						/*
						std::cout << "----------------------------------------------------------------------- vns" << std::endl;
						std::cout << bestBeta << " " << bestGamma << " " << tMin << " <----- best Beta Gamma t" << std::endl;
						std::cout << vn_true.x << " " << vn_true.y << " " << vn_true.z << " <----- vn_true before smoothing" << std::endl;
						std::cout << vn_final.x << " " << vn_final.y << " " << vn_final.z << " <----- vn_final with smoothing" << std::endl;
						std::cout << vNA.x << " " << vNA.y << " " << vNA.z << " <----- vNA" << std::endl;
						std::cout << vNB.x << " " << vNB.y << " " << vNB.z << " <----- vNB" << std::endl;
						std::cout << vNC.x << " " << vNC.y << " " << vNC.z << " <----- vNC" << std::endl;
						*/
						
						// IMP set vNX to zero
						Vector3D<double> vNEmpty;
						vNA = vNEmpty;
						vNB = vNEmpty;
						vNC = vNEmpty;
												
					}  // else
					
           
					

          /*// 2017 working second method and it uses non-transformed normals
          ptr_model->ptr_faceMat->getElementNormals(index, n1Min, n2Min, n3Min); // removed 2017 // returns 1 indexed normals //old: note it's just index and not index+1
          // But following function takes care of 1 indexing of n1Min
					ptr_model->ptr_vnsMat->getElementXYZOfVn(n1Min, vn_1.x, vn_1.y, vn_1.z);
          ptr_model->ptr_vnsMat->getElementXYZOfVn(n2Min, vn_2.x, vn_2.y, vn_2.z);
          ptr_model->ptr_vnsMat->getElementXYZOfVn(n3Min, vn_3.x, vn_3.y, vn_3.z);

          vn_final = vn_1 + (vn_2 + vn_3);
          vn_final/3;
          vn_final.normalize();*/
					// ~ std::cout << vn.x << " " << vn.y << " " << vn.z << std::endl;

					// Already normalized. Unnecessary! Can be removed to save time. ----------------->
					//~ vn.normalize();
					//~ std::cout << vn.x << " " << vn.y << " " << vn.z << std::endl;
					this->bestNormal = vn_final;
					//~ std::cout << this->bestNormal.x << " " << this->bestNormal.y << " " << this->bestNormal.z << std::endl;
				}

				// same as sphere's --------------------------
				this->bestT_model = tMin;
				this->ptr_bestModel = ptr_model;
				this->best_intersection_model = this->eye + (this->d)*tMin;
				//~ std::cout << "Ray's best intersection ------------->" << std::endl;
				//~ this->best_intersection_model.printVector3D();

				return 1;
            } //if

            else {return 0;}

/*
    //Vector3D<double> d, intersection; // d is the unit ray
    double scalarC, scalarV, scalarD, scalarTEye;

    scalarC = (ptr_sph->center - this->eye).getMagnitude();
    //std::cout << scalarC << " scalarC" << std::endl;
    //ptr_sph->center.printVector3D();
    //this->eye.printVector3D();

    //this->d.printVector3D();
    scalarV = (ptr_sph->center - this->eye).dotProduct(this->d);
    //std::cout << scalarV << " scalarV" << std::endl;
    //std::cin.get();
    //******************************************** check if center is infront of the image plane ******************************************
    //std::cout << ptr_sph->radius << " ptr_sph->radius" << std::endl;
    scalarD = (ptr_sph->radius)*(ptr_sph->radius) - (scalarC*scalarC-scalarV*scalarV); // d^2

    if(scalarD<0){
        /////std::cout << scalarD << " scalarD is less than 0" << std::endl;
        //this->ptr_bestSph = NULL;
        return 0;
    }

    else {
        /////std::cout << scalarD << " scalarD is > 0" << std::endl;
        // distance from intersection to eye
        scalarTEye = fabs(scalarV)-fabs(sqrt(scalarD)); // = (v-d) // Q-E = (scalarV-scalarD)*d
        // Update best values for the ray
        if (scalarTEye < this->bestT_sphere && scalarTEye > 0.0){
            this->bestT_sphere = scalarTEye;
            this->ptr_bestSph = ptr_sph;
            this->best_intersection_sphere = this->eye + (this->d)*scalarTEye;
            ///std::cout << "Ray's best intersection ------------->" << std::endl;
            ///this->best_intersection_sphere.printVector3D();
        } //if
        return 1;
    } // else
*/
} // rayModelIntersectionTest

// return 1 => sphere; return 2 => .obj object; return 0 => no intersection
unsigned Ray::findBestIntersectionObject(){
	if (this->ptr_bestSph && this->bestT_sphere < this->bestT_model && this->bestT_sphere >= 0) {
		this->objectType = 1;
		return 1;
	}

	else if (this->ptr_bestModel && this->bestT_model >= 0 && this->bestT_model < this->bestT_sphere){
		this->objectType = 2;
		return 2;
	}

	else {
		this->objectType = 0;
		return 0;
	} // this means ray did not intersect with both sphere and object.
}

void Ray::printRayAttributes(){
    std::cout << "Ray is ------------->" << std::endl;
    std::cout << eye.x << " " << eye.y << " " << eye.z << " <----- eye (i.e., origin)" << std::endl;
    std::cout << d.x << " " << d.y << " " << d.z << " <----- direction" << std::endl;
    std::cout << bestT_sphere << " bestT_sphere " << bestT_model << " bestT_model " << ptr_bestSph << " (" << best_intersection_sphere.x << ", " << best_intersection_sphere.y << ", " << best_intersection_sphere.z << ")" << " <----- best T, sph pointer, intersection" << std::endl;
} // constructor for printing ray attributes

Ray::~Ray(){
} // Ray object destructor

#endif // RAY_H_INCLUDED
