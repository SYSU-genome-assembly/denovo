#include <istream>
#include <ostream>
#include <cassert>

#ifndef __Matrix_h__
#define __Matrix_h__

class Vector;
class Matrix;

double sum(const Vector& vec);
Vector max(const Vector& vec1, const Vector& vec2);

Vector operator*(const Vector& vec1, const Vector& vec2);
Vector operator*(const Vector& vec, double scale);
Vector operator*(double scale, const Vector& vec);
Vector operator/(const Vector& vec, double scale);
Vector operator+(const Vector& vec1, const Vector& vec2);
Vector operator-(const Vector& vec1, const Vector& vec2);

Matrix operator*(const Matrix& mat, double scale);
Matrix operator*(double scale, const Matrix& mat);
Matrix operator/(const Matrix& mat, double scale);
Matrix operator+(const Matrix& mat1, const Matrix& mat2);

std::istream& operator>>(std::istream& is, Vector& vec);
std::istream& operator>>(std::istream& is, Matrix& mat);
std::ostream& operator<<(std::ostream& os, const Vector& vec);
std::ostream& operator<<(std::ostream& os, const Matrix& mat);

class Vector{
 private:
    int size;
    double* value;

    friend class Matrix;
    friend double sum(const Vector& vec);
	friend Vector max(const Vector& vec1, const Vector& vec2);
    friend Vector operator*(const Vector& vec1, const Vector& vec2);    
    friend Vector operator*(const Vector& vec, double scale);
    friend Vector operator*(double scale, const Vector& vec);
    friend Vector operator/(const Vector& vec, double scale);
    friend Vector operator+(const Vector& vec1, const Vector& vec2);
    friend Vector operator-(const Vector& vec1, const Vector& vec2);

    friend std::istream& operator>>(std::istream& is, Vector& vec);
    friend std::ostream& operator<<(std::ostream& os, const Vector& vec);

    
 public:
    Vector();
    Vector(int size);
    Vector(const Vector& b);
    ~Vector();

    static Vector zeros(int size);
    
    double& operator[](int i);
    double operator[](int i) const;
    Vector& operator=(const Vector& b);

    double dot(const Vector& b) const;
    Vector dot(const Matrix& b) const;
    Matrix outer(const Vector& b) const;
    
    void getValues(double* array) const;
    void setValues(const double* array);
    void cumSum(double* array) const;
    int len() const;

    void print() const;
};

class Matrix{
 private:
    int nrow, ncol;
    double* value;

    /*
    bool inverted;
    Matrix* invert;

    bool deted;
    double determinant;
    */

    friend class Vector;
    friend Matrix operator*(const Matrix& mat, double scale);
    friend Matrix operator*(double scale, const Matrix& mat);
    friend Matrix operator/(const Matrix& mat, double scale);
    friend Matrix operator+(const Matrix& mat1, const Matrix& mat2);
    

    friend std::istream& operator>>(std::istream& is, Matrix& mat);
    friend std::ostream& operator<<(std::ostream& os, const Matrix& mat);
    
 public:
    Matrix();
    Matrix(int nrow, int ncol);
    Matrix(const Matrix& b);
    ~Matrix();

    static Matrix zeros(int nrow, int ncol);

    Matrix dot(const Matrix& b) const;
    Vector dot(const Vector& b) const;
    Matrix inv() const;
    double det() const;
    
    void getValues(double* array) const;
    void setValues(const double* array);
    double* operator[](int i);
    const double* operator[](int i)const;
    Matrix& operator=(const Matrix& b);

    void print() const;
};

// inlining operator functions
using namespace std;

inline
double Vector::dot(const Vector& b) const{
    assert(size == b.size);
    
    double ret = 0.0;
    for(int i=0; i<size; i++)
	ret += value[i]*b.value[i];

    return ret;
}

inline
Vector Vector::dot(const Matrix& b) const{
    assert(size == b.nrow);
    Vector vec(b.ncol);

    for(int i=0; i<b.ncol; i++){
	vec.value[i] = 0.0;
	for(int j=0; j<size; j++)
	    vec.value[i] += value[j]*b.value[j*b.ncol+i];
    }

    return vec;
}

inline
Matrix Vector::outer(const Vector& b) const{
    Matrix mat(size, b.size);

    for(int i=0; i<size; i++)
	for(int j=0; j<b.size; j++)
	    mat[i][j] = value[i]*b.value[j];

    return mat;
}

inline
double& Vector::operator[](int i){
    assert(i>=0 && i<size);
    return value[i];
}

inline
double Vector::operator[](int i) const{
    assert(i>=0 && i<size);
    return value[i];    
}

inline
Vector& Vector::operator=(const Vector& b){
    if(this!=&b){
	if(size!=b.size){
	    if(size==0 && b.size!=0){
		size = b.size;
		value = new double[size];
	    }
	    else if(size!=0 && b.size==0){
		delete[] value;
		value = NULL;
		size = 0;
	    }
	    else{ // both are not zero
		delete[] value;
		size = b.size;
		value = new double[size];
	    }
	}
	
	for(int i=0; i<size; i++)
	    value[i] = b.value[i];
    }
    return *this;
}

inline
Matrix Matrix::dot(const Matrix& b) const{
    assert(ncol==b.nrow);
    Matrix result(nrow, b.ncol);

    for(int row=0; row<nrow; row++)
	for(int col=0; col<b.ncol; col++){
	    result[row][col] = 0.0;
	    for(int i=0; i<ncol; i++)
		result[row][col] += value[row*ncol+i]*b.value[i*b.ncol+col];
	}

    return result;
}

inline
Vector Matrix::dot(const Vector& b) const{
    assert(ncol==b.size);
    Vector result(nrow);

    for(int row=0; row<nrow; row++){
	result[row]=0.0;
	for(int i=0; i<ncol; i++)
	    result[row]+=value[row*ncol+i]*b.value[i];
    }
	
    return result;
}

    
inline    
double* Matrix::operator[](int i){
    assert(i>=0 && i<nrow);
    return &value[i*ncol];
}

inline
const double* Matrix::operator[](int i)const{
    assert(i>=0 && i<nrow);
    return &value[i*ncol];    
}

inline
Matrix& Matrix::operator=(const Matrix& b){
    if(this!=&b){
	if(nrow*ncol!=b.nrow*b.ncol){
	    if(nrow*ncol==0 && b.nrow*b.ncol!=0){
		value = new double[b.nrow*b.ncol];
	    } else if(nrow*ncol!=0 && b.nrow*b.ncol==0){
		delete[] value;
		value = NULL;
	    } else{ // both are not zero
		delete[] value;
		value = new double[b.nrow*b.ncol];
	    }
	}
	nrow = b.nrow;
	ncol = b.ncol;
	
	for(int i=0; i<nrow*ncol; i++)
	    value[i] = b.value[i];
    }
    return *this;
}


#endif
