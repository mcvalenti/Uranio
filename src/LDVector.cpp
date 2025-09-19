/*
 * LLDVector.cpp
 *
 *  Created on: 19 may. 2023
 *      Author: ceci
 */

#include "LDVector.h"
#include <string.h>
#include <stdexcept>
#include <iomanip>
#include <cstdlib>
using namespace std;



// Constructors
LDVector::LDVector()
{
	// empty LDVector
    this->vec = NULL;
    this->asize = 0;
    //ctor
}
LDVector::LDVector(long double* vec, unsigned asize)
{
	// LDVector from array
    this->vec = new long double[asize];
    this->asize = asize;
    //ctor
    memcpy(this->vec, vec, asize*sizeof(vec[0]));
}
LDVector::LDVector(unsigned asize)
{
	// LDVector initialize in zeros
    this->vec = new long double[asize];
    this->asize = asize;
    for(unsigned i=0;i<asize;i++){
        this->vec[i] = 0.0;
    }
}
LDVector::LDVector(const LDVector &other)
{
	// LDVector from another LDVector
    this->vec = new long double[other.asize];
    this->asize = other.asize;
    //ctor
    memcpy(this->vec, other.vec, other.asize*sizeof(vec[0]));
}


// Operators
const LDVector LDVector::operator+(const LDVector &other)const{
	// Sum of LDVector elements with another LDVector elements
    LDVector result(other.asize);
    for(unsigned i=0;i<this->asize;i++){
        result.vec[i] = this->vec[i] + other.vec[i];
    }
    return result;
}
const LDVector LDVector::operator-(const LDVector &other)const{
	// Sum of LDVector elements with another LDVector elements
    LDVector result(other.asize);
    for(unsigned i=0;i<this->asize;i++){
        result.vec[i] = this->vec[i] - other.vec[i];
    }
    return result;
}
LDVector& LDVector::operator=(const LDVector &other){
	// Set a LDVector equal to other
    delete [] this->vec;
    this->vec = new long double[other.asize];
    this->asize = other.asize;

    memcpy(this->vec, other.vec, asize*sizeof(vec[0]));
    return *this;
}
LDVector& LDVector::operator+=(const LDVector &other){

    if(other.asize!=this->asize)
        throw std::runtime_error("out of bounds");

    for(unsigned i=0;i<this->asize;i++){
        this->vec[i] += other.vec[i];
    }
    return *this;
}
const LDVector LDVector::operator*(long double val)const{
	// Multiply LDVector elements with long double const
    LDVector result(this->asize);
    for(unsigned i=0;i<this->asize;i++){
        result.vec[i] = this->vec[i] * val;
    }
    return result;
}
long double& LDVector::operator[](unsigned index){
	// LDVector element
    if(index>asize-1){
        throw std::runtime_error("out of bounds");
    }
    return this->vec[index];
}
const long double& LDVector::operator[](unsigned index)const{
	// LDVector element
    if(index>asize-1){
        throw std::runtime_error("out of bounds");
    }
    return this->vec[index];
}
ostream& operator<<(ostream& os, const LDVector& vec){
	// Print LDVector
	os << std::fixed;
    os << std::setprecision(10);
    for(unsigned i=0;i<vec.asize;i++){
        os<<vec.vec[i]<<"\t";
    }
    return os;
}


// Methods
long double LDVector::getValue(unsigned index){
	// LDVector elements
    return this->vec[index];
}
/*bool LDVector::areEqual(const LDVector &other, double tolerance){
	for(unsigned i=0;i<this->asize;i++){
		if (abs(vec[i]-other[i])>tolerance){
			return false;
		}
	}
	return true;
}*/
/*double LDVector::get_max_absolute(const LDVector &other){
	unsigned aux=abs(vec[0]);
	for (unsigned i=0;i<this->asize;i++){
		if (abs(vec[i])>aux){
			aux=vec[i];
		}
	}
	return aux;
}*/

LDVector::~LDVector()
{
    //dtor
    delete [] this->vec;
}





