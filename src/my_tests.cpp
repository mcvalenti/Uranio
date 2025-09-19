/*
 * my_tests.cpp
 * =======================================================================	
 * This files will check some functions used as middle steps.
 * =======================================================================	
 *  Created on: 11 ago. 2023
 *      Author: ceci
 */

#include <iostream>
#include "LDVector.h"
#include "propagators.h"

void anomaly_propagation(){
	/* 
	 * Tests: "sv_from_true_anomaly" and including "lagrange_coeff_from_true_anomaly"
	 */
	long double sv_vector[6]={8182.4,-6865.9,0,0.47572,8.8116,0};
	long double sv_result[6]={1454.9864223053, 8251.4609450803,0.0000000000,-8.1323880267,5.6785406852,-0.0000000000};
	//double tolerance=1e-10;
	LDVector sv_0(sv_vector,6);
	LDVector result(sv_result,6);
	LDVector sv_1(6);
	double delta_nu=120.0;
	sv_1=sv_from_true_anomaly(sv_0,delta_nu);
/*	if (sv_1.areEqual(result,tolerance)){
		cout << "TEST OK\n" << endl;
	}else{
		cout << "TEST NOK\n" << endl;
	}*/
	
}

void run_tests(){
	cout << "---ANOMALY PROPAGATION---";
	anomaly_propagation();
}


/*
void RK45(double* y, double t, double h){


	double tolerance = 1e-5;

	// RK4(5) derivatives

		for (size_t i = 0; i < y.size(); ++i) {
		y_next_rk4[i] += 25 * k1[i] / 216 + 1408 * k3[i] / 2565 + 2197 * k4[i] / 4104 - k5[i] / 5;
		y_next_rk5[i] += 16 * k1[i] / 135 + 6656 * k3[i] / 12825 + 28561 * k4[i] / 56430 - 9 * k5[i] / 50 + 2 * k6[i] / 55;
        }


double computeError(const std::vector<double>& y_next_rk4, const std::vector<double>& y_next_rk5, double tolerance) {
    double max_error = 0.0;
    for (size_t i = 0; i < y_next_rk4.size(); ++i) {
        double error = std::abs(y_next_rk5[i] - y_next_rk4[i]);
        max_error = std::max(max_error, error / (tolerance * std::abs(y_next_rk5[i])));
    }
    return max_error;
}


}
*/
