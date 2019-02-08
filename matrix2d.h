//===================================================================================
// Name        : matrix2d.h
// Created on  : 19 October 2016
// Description : Class for 2D matrix operations
//===================================================================================
#ifndef MATRIX2D_H_INCLUDED
#define MATRIX2D_H_INCLUDED

#include <iostream>
#include <vector>
#include <iomanip>  // std::setprecision(2)

#include "vector3d.h"
//#include "matrix2d.cpp"

template <class T>
class Matrix2D{
private:
    // created matrix as vector of vectors
    std::vector<std::vector<T> > matrix; // IMP: space between > > is necessary !!!!!

public:

    unsigned rows;
    unsigned columns;
    T determinant;  // available only for 3x3 matrix

    Matrix2D(unsigned n_rows, unsigned n_columns, const T& value);
    virtual ~Matrix2D();

    // Get the number of rows or columns ----------------------------
    unsigned getN_Rows() const;
    unsigned getN_Columns() const;

    // Get min and max element along row or column
//    void getMinElementArray(const unsigned& rOrC, const unsigned& index) const;

    // Access Individual Matrix Element1 -----------------------------
    const T& operator()(const unsigned &rowN, const unsigned &colN) const;
    // Access Individual Matrix Element2 -----------------------------
    T& operator()(const unsigned &rowN, const unsigned &colN);

    // Get or Set the elements of matrix using indexing (rows, columns)
    void setElement(const unsigned& row, const unsigned& col, const double value);
    T& getElement(const unsigned& row, const unsigned& col);
    void getElementXYZ(const unsigned& face, T& vertex0, T& vertex1, T& vertex2); //saves 9 seconds on 145 seconds
    void getElementXYZOfVn(const unsigned& index, T& vn0, T& vn1, T& vn2);
    void getElementVertices(const unsigned& face, T& vertex0, T& vertex1, T& vertex2); //saves 9 seconds on 145 seconds
    void getElementNormals(const unsigned& face, T& vertex0, T& vertex1, T& vertex2); //saves 9 seconds on 145 seconds
    void getElementXYZTriangle(const unsigned& vertex0, const unsigned& vertex1, const unsigned& vertex2, T& value00, T& value01, T& value02, T& value10, T& value11, T& value12, T& value20, T& value21, T& value22);
    void getElementXYZTriangleA4(const unsigned& vertex0, const unsigned& vertex1, const unsigned& vertex2, T& value00, T& value01, T& value02, T& value10, T& value11, T& value12, T& value20, T& value21, T& value22);
    //std::cout << this->matrix[i][j] << " <------- this->matrix[i][j] " << std::endl;

    // Assignment Operation -----------------------------------------
    Matrix2D<T>& operator=(const Matrix2D<T>& M2);

    // Scalar Operations on Matrix ----------------------------------
    Matrix2D<T> operator*(const T& M2);
    Matrix2D<T> operator/(const T& M2);

    // Arithmatic Operations ----------------------------------------
    Matrix2D<T> operator+(const Matrix2D<T>& M2);
    Matrix2D<T> operator-(const Matrix2D<T>& M2);
    Matrix2D<T> operator*(const Matrix2D<T>& M2);
	Vector3D<T> operator*(const Vector3D<T>& V2);

    // Print the entire matrix --------------------------------------
    void printMatrix();
};

// Matrix2D Constructor
template<class T>
Matrix2D<T>::Matrix2D(unsigned n_rows, unsigned n_columns, const T& value){
	unsigned i;
	//std::cout << n_rows << " " << n_columns << " " << value  << " <------- dimensions " << std::endl;
    matrix.resize(n_rows);
    for (i=0; i<matrix.size(); i++){
        matrix[i].resize(n_columns, value);
    }
    rows = n_rows;
    columns = n_columns;
} // constructor

// Set Matrix Element
template<class T>
void Matrix2D<T>::setElement(const unsigned& row, const unsigned& col, const double value){ //IMP removed () after setElements ***** (error: function returning a function
    this->matrix[row][col] = value;
    //return this->matrix[row][col];
}

// Get Matrix Element
template<class T>
T& Matrix2D<T>::getElement(const unsigned& row, const unsigned& col){
    //std::cout << this->matrix[i][j] << " <------- this->matrix[i][j] " << std::endl;
    return this->matrix[row][col];
}

template<class T>
void Matrix2D<T>::getElementXYZ(const unsigned& face, T& vertex0, T& vertex1, T& vertex2){
    //std::cout << this->matrix[i][j] << " <------- this->matrix[i][j] " << std::endl;
    vertex0 = this->matrix[face][1]; // start from 1 because column 0 is "the number of vertices for that face"
    vertex1 = this->matrix[face][2];
    vertex2 = this->matrix[face][3];
}

template<class T>
void Matrix2D<T>::getElementXYZOfVn(const unsigned& index, T& vn0, T& vn1, T& vn2){
    //std::cout << this->matrix[i][j] << " <------- this->matrix[i][j] " << std::endl;
    vn0 = this->matrix[index-1][0];
    vn1 = this->matrix[index-1][1];
    vn2 = this->matrix[index-1][2];
}

template<class T>
void Matrix2D<T>::getElementVertices(const unsigned& face, T& vertex0, T& vertex1, T& vertex2){
    // 0 indexing since facemat stores with zero indexing
    // vertex0 => type = unsigned
    //std::cout << this->matrix[i][j] << " <------- this->matrix[i][j] " << std::endl;
    vertex0 = this->matrix[face][0]; //
    vertex1 = this->matrix[face][3];
    vertex2 = this->matrix[face][6];
}

template<class T>
void Matrix2D<T>::getElementNormals(const unsigned& face, T& normal0, T& normal1, T& normal2){
    // ALL THREE VALUES ARE SAME
    //std::cout << this->matrix[i][j] << " <------- this->matrix[i][j] " << std::endl;
    normal0 = this->matrix[face][2]; // 
    normal1 = this->matrix[face][5];
    normal2 = this->matrix[face][8];
}

// for normal depth tracing => assignment 3 || for ply format vertices start from zeroth index
template<class T>
void Matrix2D<T>::getElementXYZTriangle(const unsigned& vertex0, const unsigned& vertex1, const unsigned& vertex2, T& value00, T& value01, T& value02, T& value10, T& value11, T& value12, T& value20, T& value21, T& value22){
    //std::cout << this->matrix[i][j] << " <------- this->matrix[i][j] " << std::endl;
    value00 = this->matrix[vertex0][0];
    value01 = this->matrix[vertex0][1];
    value02 = this->matrix[vertex0][2];

    value10 = this->matrix[vertex1][0];
    value11 = this->matrix[vertex1][1];
    value12 = this->matrix[vertex1][2];

    value20 = this->matrix[vertex2][0];
    value21 = this->matrix[vertex2][1];
    value22 = this->matrix[vertex2][2];
}

// for depth tracing with lights => assignment 4 || for .obj format, vertices start from first index
template<class T>
void Matrix2D<T>::getElementXYZTriangleA4(const unsigned& vertex0, const unsigned& vertex1, const unsigned& vertex2, T& value00, T& value01, T& value02, T& value10, T& value11, T& value12, T& value20, T& value21, T& value22){
    //std::cout << this->matrix[i][j] << " <------- this->matrix[i][j] " << std::endl;
    value00 = this->matrix[vertex0-1][0];
    value01 = this->matrix[vertex0-1][1];
    value02 = this->matrix[vertex0-1][2];

    value10 = this->matrix[vertex1-1][0];
    value11 = this->matrix[vertex1-1][1];
    value12 = this->matrix[vertex1-1][2];

    value20 = this->matrix[vertex2-1][0];
    value21 = this->matrix[vertex2-1][1];
    value22 = this->matrix[vertex2-1][2];
}

//void getMinElementArray(const unsigned& rOrC, const unsigned& index) const{
//
//
//};


// Access Individual Element
template<class T>
T& Matrix2D<T>::operator()(const unsigned &rowN, const unsigned &colN){ // const{
    //std::cout << this->matrix[i][j] << " <------- this->matrix[i][j] " << std::endl;
    //return this->matrix[rowN][colN];
    return this->matrix[rowN][colN];
}

// Access Individual Element with a Const output type
template<class T>
const T& Matrix2D<T>::operator()(const unsigned &rowN, const unsigned &colN) const{ // {
    //std::cout << this->matrix[i][j] << " <------- this->matrix[i][j] " << std::endl;
    return this->matrix[rowN][colN];
}

// Get number of rows
template<class T>
unsigned Matrix2D<T>::getN_Rows() const{
    return this->rows;
}

// Get number of columns
template<class T>
unsigned Matrix2D<T>::getN_Columns() const{
    return this->columns;
}

// Assignment Operation
template<class T>
Matrix2D<T>& Matrix2D<T>::operator=(const Matrix2D<T>& M2){
	unsigned  i, j;
    unsigned nRows;
    unsigned nCols;
    nRows = M2.getN_Rows();
    nCols = M2.getN_Columns();

    unsigned matSize = matrix.size();
    //std::cout << matSize << " <------- matSize " << std::endl;

    // First resize matrix to match the dimensions of M2
    matrix.resize(nRows);
    //matrix.resize(nCols);
    for (i=0; i<matSize; i++){
        //matrix[i].resize(nRows);
        matrix[i].resize(nCols);
    }
    //this->printMatrix();

    for (i=0; i<nRows; i++){
        for (j=0; j<nCols; j++){
            matrix[i][j] = M2(i,j);
        }
    }
    //std::cout << nRows << " " << nCols << " <------- dimensions " << std::endl;
    rows = nRows;
    columns = nCols;

    return *this;
}

// Scalar Operations ------------------------------------------------------------
// Scalar Multiplication
template<class T>
Matrix2D<T> Matrix2D<T>::operator*(const T& M2){
	unsigned  i, j;
    Matrix2D output(rows, columns, 0); // for storing output matrix
	//std::cout << rows << " " << cols << " <------- dimensions " << std::endl;

    for (i=0; i<rows; i++){
        for (j=0; j<columns; j++){
            output(i,j) = this->matrix[i][j] * M2;
        }
    }
    //output.printMatrix();
    return output;
}

// Scalar Division
template<class T>
Matrix2D<T> Matrix2D<T>::operator/(const T& M2){
	unsigned  i, j;
    Matrix2D output(rows, columns, 0);
    //std::cout << rows << " " << cols << " <------- dimensions " << std::endl;

    for (i=0; i<rows; i++){
        for (j=0; j<columns; j++){
            output(i,j) = this->matrix[i][j]/M2;
        }
    }
    //output.printMatrix();
    return output;
}

// Arithmatic Operations ------------------------------------------------------------
// Matrix Addition
template<class T>
Matrix2D<T> Matrix2D<T>::operator+(const Matrix2D<T>& M2){
	unsigned i, j;
	//std::cout << rows << " " << cols << " <------- dimensions " << std::endl;
    Matrix2D output(rows, columns, 0);

    for (i=0; i<rows; i++){
        for (j=0; j<columns; j++){
            output(i,j) = this->matrix[i][j] + M2(i,j);
            //std::cout << this->matrix[i][j] << " <------- this->matrix[i][j] " << std::endl;
        }
    }
    //this->printMatrix();
    //output.printMatrix();
    return output;
}

// Matrix Subtraction
template<class T>
Matrix2D<T> Matrix2D<T>::operator-(const Matrix2D<T>& M2){
	unsigned  i, j;
    unsigned nRows;
    unsigned nColumns;
    nRows = M2.getN_Rows();
    nColumns = M2.getN_Columns();
    //std::cout << nRows << " " << nColumns << " <------- dimensions " << std::endl;
    Matrix2D output(nRows, nColumns, 0);

    for (i=0; i<nRows; i++){
        for (j=0; j<nColumns; j++){
            output(i,j) = this->matrix[i][j] - M2(i,j);
        }
    }
    //output.printMatrix();
    return output;
}

// Matrix Multiplication: output = matrix*M2
template<class T>
Matrix2D<T> Matrix2D<T>::operator*(const Matrix2D<T>& M2){
	unsigned  i, j, k;
    unsigned nRows;
    unsigned nColumns;
    nRows = M2.getN_Rows();
    nColumns = M2.getN_Columns();
    //std::cout << nRows << " " << nColumns << " <------- dimensions " << std::endl;

    // check if number of columns of this matrix is equal to number of rows of M2
    if (this->columns != nRows){
        std::cout << "Dimension Mismatch! Number of columns of first matrix must be equal to number of rows on the second matrix! " << std::endl;
    }

    Matrix2D output(this->rows, nColumns, 0); // Very important! Note this step! Check the dimensions!

    for (i=0; i<this->rows; i++){
        for (j=0; j<nColumns; j++){
            for (k=0; k<nRows; k++){
                output(i,j) += this->matrix[i][k] * M2(k,j);
                //std::cout << output(i,j) << " -- i -> " << i << " -- j --> " << j  << " -- k --> " << k << std::endl;
            } // for i
        } // for j
    } // for k
    //output.printMatrix();
    return output;
}

template<class T>
Vector3D<T> Matrix2D<T>::operator*(const Vector3D<T>& V2){
	unsigned i, j;
	Vector3D<T> output(0,0,0);
	//output.printVector3D();

	output.x = this->matrix[0][0] * V2.x + this->matrix[0][1] * V2.y + this->matrix[0][2] * V2.z;
	output.y = this->matrix[1][0] * V2.x + this->matrix[1][1] * V2.y + this->matrix[1][2] * V2.z;
	output.z = this->matrix[2][0] * V2.x + this->matrix[2][1] * V2.y + this->matrix[2][2] * V2.z;

	return output;
}

// Print Entire Matrix
template<class T>
void Matrix2D<T>::printMatrix(){
	unsigned i, j;
    std::cout << "Matrix is ------------------> " << std::endl;
    for (i=0; i<this->rows; i++){
        for (j=0; j<this->columns; j++){
            std::cout << std::setprecision(3) << std::fixed << this->matrix[i][j] << ", ";
        }
        std::cout << std::endl;
    }
}

// Matrix2D Destructor
template<class T>
Matrix2D<T>::~Matrix2D(){
} // destructor

#endif // MATRIX2D_H_INCLUDED


//========================================================================================
// Referred Library  : http://www.stroustrup.com/ Programming/Matrix/Matrix.h
//             : https://www.quantstart.com/ articles/ Matrix-Classes-in-C-The-Header-File
//========================================================================================
