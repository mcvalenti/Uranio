#include "DVector.h"
#include <string.h>
#include <stdexcept>

DVector::DVector()
{
	// empty Dvector
    this->vec = NULL;
    this->asize = 0;
    //ctor
}

double DVector::getValue(unsigned index){
	// DVector elements
    return this->vec[index];
}

DVector::DVector(double* vec, unsigned asize)
{
	// DVector from array
    this->vec = new double[asize];
    this->asize = asize;
    //ctor
    memcpy(this->vec, vec, asize*sizeof(vec[0]));
}

DVector::DVector(unsigned asize)
{
	// DVector initialize in zeros
    this->vec = new double[asize];
    this->asize = asize;
    for(unsigned i=0;i<asize;i++){
        this->vec[i] = 0.0;
    }
    //ctor

}

DVector::DVector(const DVector &other)
{
	// DVector from another DVector
    this->vec = new double[other.asize];
    this->asize = other.asize;
    //ctor
    memcpy(this->vec, other.vec, other.asize*sizeof(vec[0]));
}

const DVector DVector::operator+(const DVector &other)const{
	// Sum of DVector elements with another DVector elements
    DVector result(other.asize);
    for(unsigned i=0;i<this->asize;i++){
        result.vec[i] = this->vec[i] + other.vec[i];
    }
    return result;
}


DVector& DVector::operator=(const DVector &other){
	// Set a DVector equal to other
    delete [] this->vec;
    this->vec = new double[other.asize];
    this->asize = other.asize;

    memcpy(this->vec, other.vec, asize*sizeof(vec[0]));
    return *this;
}


DVector& DVector::operator+=(const DVector &other){

    if(other.asize!=this->asize)
        throw std::runtime_error("out of bounds");

    for(unsigned i=0;i<this->asize;i++){
        this->vec[i] += other.vec[i];
    }
    return *this;
}


const DVector DVector::operator*(double val)const{
	// Multiply DVector elements with double const
    DVector result(this->asize);
    for(unsigned i=0;i<this->asize;i++){
        result.vec[i] = this->vec[i] * val;
    }
    return result;
}

double& DVector::operator[](unsigned index){
	// DVector element
    if(index>asize-1){
        throw std::runtime_error("out of bounds");
    }
    return this->vec[index];
}

const double& DVector::operator[](unsigned index)const{
	// DVector element
    if(index>asize-1){
        throw std::runtime_error("out of bounds");
    }
    return this->vec[index];
}

ostream& operator<<(ostream& os, const DVector& vec){
	// Print DVector
    for(unsigned i=0;i<vec.asize;i++){
        os<<vec.vec[i]<<"\t";
    }
    return os;
}



DVector::~DVector()
{
    //dtor
    delete [] this->vec;
}
