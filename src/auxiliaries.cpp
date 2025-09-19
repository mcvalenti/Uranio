/*
 * auxiliaries.cpp
 *
 *  Created on: 19 may. 2023
 *      Author: ceci
 */

#ifndef AUXILIARIES_CPP_
#define AUXILIARIES_CPP_

#include "LDVector.h"
#include <cmath>

double period(long double& semimajorAxis){
    /*
    -Computes orbital period in seconds-
    semimajorAxis: float [km]
     */
	double period;
	double mu=398600.448;
	if (semimajorAxis>0){
		period=2*M_PI*sqrt(semimajorAxis*semimajorAxis*semimajorAxis/mu);
	}
	return period;
};

LDVector cross(LDVector& a, LDVector& b){
	LDVector cross_prod(3);
	cross_prod[0]=a[1]*b[2]-b[1]*a[2];
	cross_prod[1]=-(a[0]*b[2]-b[0]*a[2]);
	cross_prod[2]=a[0]*b[1]-b[0]*a[1];
	return cross_prod;	
}


long double dot3(LDVector& a,LDVector& b){
	// Just for LDVectors of dim=3.
	// TO DO!! Check a and b dimensions.
	long double dot_prod;
	dot_prod=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	return dot_prod;
}



#endif /* AUXILIARIES_CPP_ */
