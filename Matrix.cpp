#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_linalg.h"

#include "Matrix.h"

using namespace std;

Vector::Vector(){
    this->size = 0;
    this->value = NULL;
}

Vector::Vector(int size){
    assert(size>=0);
    this->size = size;
    this->value = NULL;
    if(size>0)
	this->value = new double[size];
}

Vector::Vector(const Vector& b){
    this->size = b.size;
    this->value = NULL;

    if(size>0){
	this->value = new double[size];

	for(int i=0; i<size; i++)
	    value[i] = b.value[i];
    }
}

Vector::~Vector(){
    if(value!=NULL)
	delete[] value;
}



int Vector::len() const{
    return this->size;
}

void Vector::print() const{
    printf("[");    
    for(int i=0; i<size; i++)
	printf("%lf ", value[i]);
    printf("]\n");
}

Vector Vector::zeros(int size){
    assert(size>=0);
    Vector result(size);
    for(int i=0; i<size; i++)
	result.value[i]=0;

    return result;
}


void Vector::getValues(double* array) const{
    for(int i=0; i<size; i++)
	array[i] = value[i];
}

void Vector::setValues(const double* array){
    for(int i=0; i<size; i++)
	value[i] = array[i];
}

void Vector::cumSum(double* array) const{
	array[0] = value[0];
	for(int i=1; i<size; i++)
	array[i] = array[i-1] + value[i];
}


Matrix::Matrix(){
    this->nrow = 0;
    this->ncol = 0;
    /*
    this->inverted = false;
    this->invert = NULL;
    this->deted = false;
    this->determinant = 0.0;
    */

    this->value = NULL;
}

Matrix::Matrix(int nrow, int ncol){
    assert(nrow>=0 && ncol>=0);
    this->nrow = nrow;
    this->ncol = ncol;
    this->value = NULL;
    /*
    this->inverted = false;
    this->invert = NULL;
    this->deted = false;
    this->determinant = 0.0;
    */

    if(nrow*ncol>0)
	this->value = new double[nrow*ncol];
}

Matrix::Matrix(const Matrix& b){
    assert(b.nrow>=0 && b.ncol>=0);
    nrow = b.nrow;
    ncol = b.ncol;

    this->value = NULL;
    if(nrow*ncol>0){
	this->value = new double[nrow*ncol];
	for(int i=0; i<nrow*ncol; i++)
	    value[i] = b.value[i];
    }
    
}

Matrix::~Matrix(){
    if(value!=NULL)
	delete[] value;
}

Matrix Matrix::zeros(int nrow, int ncol){
    assert(nrow>=0 && ncol>=0);
    Matrix result(nrow, ncol);
    for(int i=0; i<nrow*ncol; i++)
	result.value[i] = 0.0;

    return result;
}



Matrix Matrix::inv() const{    
    Matrix invert(nrow, ncol);
    gsl_matrix* mat = gsl_matrix_alloc(nrow, ncol);
    gsl_matrix* inv_mat = gsl_matrix_alloc(nrow, ncol);
    gsl_permutation* perm = gsl_permutation_alloc(ncol);
    int signum;

    for(int i=0; i<nrow*ncol; i++)
	mat->data[i] = value[i];

    // LU decomposition
    gsl_linalg_LU_decomp(mat, perm, &signum);
	
    // inverse from LU
    gsl_linalg_LU_invert(mat, perm, inv_mat);

    // set return value
    invert.setValues(inv_mat->data);

    // deallocate all gsl object
    gsl_matrix_free(mat);
    gsl_matrix_free(inv_mat);
    gsl_permutation_free(perm);
    
    return invert;
}

double Matrix::det() const{
    assert(ncol==nrow);
    gsl_matrix* mat = gsl_matrix_alloc(nrow, ncol);
    gsl_permutation* perm = gsl_permutation_alloc(ncol);
    int signum;

    for(int i=0; i<nrow*ncol; i++)
	mat->data[i] = value[i];

    // LU decomposition
    gsl_linalg_LU_decomp(mat, perm, &signum);
	
    // inverse from LU
    double determinant = gsl_linalg_LU_det(mat, signum);

    // deallocate all gsl object
    gsl_matrix_free(mat);
    gsl_permutation_free(perm);
    
    return determinant;
}


void Matrix::print() const{
    printf("[");
    for(int r=0; r<nrow; r++){
	printf("[");
	for(int c=0; c<ncol; c++)
	    printf("%lf ", value[r*ncol + c]);
	printf("]\n");
    }
    printf("]\n");
}

void Matrix::getValues(double* array) const{
    for(int i=0; i<nrow*ncol; i++)
	array[i] = value[i];
}

void Matrix::setValues(const double* array){
    for(int i=0; i<nrow*ncol; i++)
	value[i] = array[i];
}

double sum(const Vector& vec){
    double s = 0.0;
    for(int i=0; i<vec.size; i++)
	s+=vec.value[i];
    return s;
}

Vector max(const Vector& vec1, const Vector& vec2){
	assert(vec1.size==vec2.size);
	Vector ret(vec1);

	for(int i=0; i<vec1.size; i++)
		if(vec2.value[i]>vec1.value[i])
			ret.value[i] = vec2.value[i];

	return ret;
}

Vector operator*(const Vector& vec1, const Vector& vec2){
    assert(vec1.size==vec2.size);
    Vector ret(vec1);

    for(int i=0; i<vec1.size; i++)
	ret.value[i]*=vec2.value[i];
    
    return ret;
}

Vector operator*(const Vector& vec, double scale){
    Vector res(vec);

    for(int i=0; i<res.size; i++)
	res.value[i]*=scale;

    return res;
}

Vector operator*(double scale, const Vector& vec){
    return operator*(vec, scale);
}

Vector operator/(const Vector& vec, double scale){
    return operator*(1.0/scale, vec);
}

Vector operator+(const Vector& vec1, const Vector& vec2){
    assert(vec1.size==vec2.size);
    Vector ret(vec1);
    
    for(int i=0; i<vec1.size; i++)
	ret.value[i] += vec2.value[i];

    return ret;
}

Vector operator-(const Vector& vec1, const Vector& vec2){
    assert(vec1.size==vec2.size);
    Vector ret(vec1);
    
    for(int i=0; i<vec1.size; i++)
	ret.value[i] -= vec2.value[i];

    return ret;
}


Matrix operator*(const Matrix& mat, double scale){
    Matrix res(mat);

    for(int i=0; i<res.nrow*res.ncol; i++)
	res.value[i]*=scale;

    return res;
}

Matrix operator*(double scale, const Matrix& mat){
    return operator*(mat, scale);
}

Matrix operator/(const Matrix& mat, double scale){
    return operator*(mat, 1.0/scale);
}

Matrix operator+(const Matrix& mat1, const Matrix& mat2){
    assert(mat1.nrow==mat2.nrow && mat1.ncol==mat2.ncol);
    Matrix res(mat1);

    for(int i=0; i<res.nrow*res.ncol; i++)
	res.value[i] += mat2.value[i];

    return res;
}

istream& operator>>(istream& is, Vector& vec){
    int size;
    is >> size;
    if(!is.good())
	return is;

    Vector input(size);
    for(int i=0; i<size; i++)
	is >> input.value[i];
    vec = input;

    return is;
}

istream& operator>>(istream& is, Matrix& mat){
    int nrow, ncol;

    is >> nrow >> ncol;
    if(!is.good())
	return is;
    
    Matrix input(nrow, ncol);

    for(int i=0; i<nrow*ncol; i++)
	is >> input.value[i];
    mat = input;

    return is;
}

ostream& operator<<(ostream& os, const Vector& vec){
    os << vec.size << endl;
    os.setf(ios::scientific);    
    os.precision(16);
    for(int i=0; i<vec.size; i++)
	os <<  vec.value[i] << ' ';
    os << endl;
    return os;
}

ostream& operator<<(ostream& os, const Matrix& mat){
    os << mat.nrow << ' ' << mat.ncol << endl;
    os.setf(ios::scientific);    
    os.precision(16);
    for(int row=0; row<mat.nrow; row++){
	for(int col=0; col<mat.ncol; col++)
	    os << mat.value[row*mat.ncol + col] << ' ';
	os << endl;
    }
    return os;
}

