#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "sceneFileParser.h"

// find the total light sources in the scene file
void getNumberOfLightSourcesSpheresModels(std::string inputFileName, unsigned& numberOfLightSources, unsigned& numberOfSpheres, unsigned& numberOfModels){
    std::string line;
    std::string searchKeywordLight("light "); // IMPORTANT: SPACE is necessary after the keyword
    std::string searchKeywordSphere ("sphere "); // IMPORTANT: SPACE is necessary after the keyword
    std::string searchKeywordModel ("model "); // IMPORTANT: SPACE is necessary after the keyword
    std::ifstream sceneFile(inputFileName.c_str());

    while (1)
        {
            getline(sceneFile, line);

            if (line.compare(0, searchKeywordLight.length(), searchKeywordLight) == 0){
                //std::cout << line << " light is present in the line " << std::endl;
                numberOfLightSources++;
            }

            if (line.compare(0, searchKeywordSphere.length(), searchKeywordSphere) == 0){
                //std::cout << line << " sphere is present in the line " << std::endl;
                numberOfSpheres++;
            }

            if (line.compare(0, searchKeywordModel.length(), searchKeywordModel) == 0){
                //std::cout << line << " model is present in the line " << std::endl;
                numberOfModels++;
            }

            //std::cin.get();
            if(sceneFile.eof()){break;} // break if end of file
        } // while1
    //std::cout << numberOfLightSources << " numberOfLightSources " << '\n';
    sceneFile.close();
    //return numberOfLightSources;
}

void getAmbientRGBVector(std::string inputFileName, Vector3D<double> &ambient){

    std::string line;
    std::string searchKeyword ("ambient "); // IMPORTANT: SPACE is necessary after the keyword
    std::string term1; // to store "ambient" word
    std::ifstream sceneFile(inputFileName.c_str());

    while (1)
        {
            getline(sceneFile, line);
            if (line.compare(0, searchKeyword.length(), searchKeyword) == 0){
                std::cout << line << " FOUND here " << '\n';
                std::stringstream ss(line);
                ss >> term1 >> ambient.x >> ambient.y >> ambient.z;
            }
            if(sceneFile.eof()){break;} // not redundant
        } // while1

    // Clipping Ambient RGB values to 1 ----------------------------------------------------------------------
    if (ambient.x > 1){
        std::cout << " WARNING::: AMBIENT LIGHT: R value is greater than 1. Clipping R to 1" << std::endl;
        ambient.x = 1;
    }
    if (ambient.y > 1){
        std::cout << " WARNING::: AMBIENT LIGHT: G value is greater than 1. Clipping G to 1" << std::endl;
        ambient.y = 1;
    }
    if (ambient.z > 1){
        std::cout << " WARNING::: AMBIENT LIGHT: B value is greater than 1. Clipping B to 1" << std::endl;
        ambient.z = 1;
    }
    // ------------------------------------------------------------------------------------------------------

    sceneFile.close();
}

void getRecursionLevel(std::string inputFileName, unsigned &recursionLevel){

    std::string line;
    std::string searchKeyword ("recursionLevel "); // IMPORTANT: SPACE is necessary after the keyword
    std::string term1; // to store "recursionLevel" word
    std::ifstream sceneFile(inputFileName.c_str());

    while (1)
        {
            getline(sceneFile, line);
            if (line.compare(0, searchKeyword.length(), searchKeyword) == 0){
                std::cout << line << " FOUND here " << '\n';
                std::stringstream ss(line);
                ss >> term1 >> recursionLevel;
            }
            if(sceneFile.eof()){break;} // not redundant
        } // while1
    sceneFile.close();
}

Light* getLightAttributes(std::string inputFileName, unsigned& lightNumber){

    Light* ptrLight;
    //std::cout << lightNumber << " inside light" << std::endl;
    ptrLight = new Light(0.2, 2, 3, 5, 7, 9, 10); // Initialize to arbitrary values

    std::string line;
    std::string searchKeyword ("light "); // IMPORTANT: SPACE is necessary after the keyword
    std::string term1; // to store "light" word
    unsigned counter = 0;
    std::ifstream sceneFile(inputFileName.c_str());

    while (counter<=lightNumber)
        {
            getline(sceneFile, line);
            if (line.compare(0, searchKeyword.length(), searchKeyword) == 0){
                //std::cout << "light is present in the line " << std::endl;
                counter++;
            }
            if(sceneFile.eof()){break;} // redundant
        } // while1

    std::cout << line << " FOUND here " << '\n';
    std::stringstream ss(line);
    ss >> term1 >> ptrLight->position.x >> ptrLight->position.y >> ptrLight->position.z >> ptrLight->w >> ptrLight->rgb.x >> ptrLight->rgb.y >> ptrLight->rgb.z;

    // Handling w=0 case where light source resides at infinity in that direction  -------------------
    if (ptrLight->w == 0){
		std::cout << " WARNING::: LIGHT: Since w=0 => source resides at infinity!" << std::endl;
		ptrLight->position = ptrLight->position + ptrLight->position*1000000000000000;  // 10000000 is assumed to be infinty in given scene context
	}
    // -----------------------------------------------------------------------------------------------

    // Clipping RGB values to 1 ----------------------------------------------------------------------
    if (ptrLight->rgb.x > 1){
        std::cout << " WARNING::: LIGHT: R value is greater than 1. Clipping R to 1" << std::endl;
        ptrLight->rgb.x = 1;
    }
    if (ptrLight->rgb.y > 1){
        std::cout << " WARNING::: LIGHT: G value is greater than 1. Clipping G to 1" << std::endl;
        ptrLight->rgb.y = 1;
    }
    if (ptrLight->rgb.z > 1){
        std::cout << " WARNING::: LIGHT: B value is greater than 1. Clipping B to 1" << std::endl;
        ptrLight->rgb.z = 1;
    }
    // -----------------------------------------------------------------------------------------------

    sceneFile.close();
    return ptrLight;
}

Sphere* getSphereAttributes(std::string inputFileName, unsigned& sphereNumber){

    Sphere* ptrSphere;
    ptrSphere = new Sphere(0.2, 2, 3, 5, 7, 9, 10); // Initialize to arbitrary values

    double rhoterm = 0.9;
    std::string line;
    std::string searchKeyword ("sphere "); // IMPORTANT: SPACE is necessary after the keyword
    std::string term1; // to store "sphere" word
    unsigned counter = 0;
    std::ifstream sceneFile(inputFileName.c_str());

    while (counter<=sphereNumber)
        {
            getline(sceneFile, line);
            if (line.compare(0, searchKeyword.length(), searchKeyword) == 0){
                counter++;
            }
            if(sceneFile.eof()){break;} // redundant
        } // while1

    std::cout << line << " FOUND here " << '\n';
    std::stringstream ss(line);
    // GET ka, kd, ks, kr ----------------->
    ss >> term1 >> ptrSphere->center.x >> ptrSphere->center.y >> ptrSphere->center.z >> ptrSphere->radius >> ptrSphere->ka.x >> ptrSphere->ka.y >> ptrSphere->ka.z
    >> ptrSphere->kd.x >> ptrSphere->kd.y >> ptrSphere->kd.z >> ptrSphere->ks.x >> ptrSphere->ks.y >> ptrSphere->ks.z >> ptrSphere->kr.x >> ptrSphere->kr.y >> ptrSphere->kr.z;

    ptrSphere->rho = rhoterm;

    //~ // Clipping RGB values to 1 ----------------------------------------------------------------------
    //~ if (ptrSphere->rgb.x > 1){
        //~ std::cout << " WARNING::: SPHERE: R value is greater than 1. Clipping R to 1" << std::endl;
        //~ ptrSphere->rgb.x = 1;
    //~ }
    //~ if (ptrSphere->rgb.y > 1){
        //~ std::cout << " WARNING::: SPHERE: G value is greater than 1. Clipping G to 1" << std::endl;
        //~ ptrSphere->rgb.y = 1;
    //~ }
    //~ if (ptrSphere->rgb.z > 1){
        //~ std::cout << " WARNING::: SPHERE: B value is greater than 1. Clipping B to 1" << std::endl;
        //~ ptrSphere->rgb.z = 1;
    //~ }
    //~ // -----------------------------------------------------------------------------------------------

    //~ // CHANGE THIS LATER ---------------------------------
    //~ // ka is a vector of 3 elements
    //~ ptrSphere->ka.x = 0.5;
    //~ ptrSphere->ka.y = 0.5;
    //~ ptrSphere->ka.z = 0.5;

    //~ // kd is a vector of 3 elements
    //~ // ks is a SCALAR but we took advantage of the vector stucture to make it vector as (ks, ks, ks)
    //~ ptrSphere->kd.x = ptrSphere->ks.x = ptrSphere->rgb.x;
    //~ ptrSphere->kd.y = ptrSphere->ks.y = ptrSphere->rgb.y;
    //~ ptrSphere->kd.z = ptrSphere->ks.z = ptrSphere->rgb.z;

    //~ // kr is a vector of 3 elements
    //~ ptrSphere->kr.x = 1.0;
    //~ ptrSphere->kr.y = 1.0;
    //~ ptrSphere->kr.z = 1.0;
    //~ // ---------------------------------------------------

    sceneFile.close();
    return ptrSphere;
}

Model* parseModelLine(std::string inputFileName, unsigned& modelNumber, std::string& modelName){

    Model* ptrModel;
    ptrModel = new Model(0,0,0,0);

    double s, tx, ty, tz, rx, ry, rz, angle;
    std::string line;
    std::string searchKeyword ("model "); // IMPORTANT: SPACE is necessary after the keyword
    std::string term1; // to store "model" word
    unsigned counter = 0;
    std::ifstream sceneFile(inputFileName.c_str());

    while (counter<=modelNumber)
        {
            getline(sceneFile, line);
            if (line.compare(0, searchKeyword.length(), searchKeyword) == 0){
                counter++;
                // std::cout << counter << " model counter " << std::endl;
            }
            if(sceneFile.eof()){break;} // redundant
        } // while1

    if (line.compare(0, searchKeyword.length(), searchKeyword) == 0){
        std::cout << line << " FOUND here " << '\n';
        std::stringstream ss(line);
        ss >> term1 >> rx >> ry >> rz >> angle >> s >> tx >> ty >> tz >> modelName;
    }

    ptrModel->scaleVector.x = s;
    ptrModel->scaleVector.y = s;
    ptrModel->scaleVector.z = s;
    std::cout << "Scaling vector (x, y, z) " << std::endl;
    ptrModel->scaleVector.printVector3D("");

    ptrModel->translationVector.x = tx;
    ptrModel->translationVector.y = ty;
    ptrModel->translationVector.z = tz;
    std::cout << "Translation vector (x, y, z) " << std::endl;
    ptrModel->translationVector.printVector3D("");

    ptrModel->rotationAxis.x = rx;
    ptrModel->rotationAxis.y = ry;
    ptrModel->rotationAxis.z = rz;
    std::cout << "Rotation axis (x, y, z) " << std::endl;
    ptrModel->rotationAxis.printVector3D("");

    ptrModel->rotationAngleDegrees = angle;
    std::cout << ptrModel->rotationAngleDegrees << " <- ptrModel->rotationAngleDegrees " << std::endl;

    ptrModel->modelName = modelName;
    std::cout << modelName << " <- modelName " << std::endl;
    sceneFile.close();
    return ptrModel;
}
