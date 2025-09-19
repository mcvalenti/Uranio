/*
 * orbit_dynamics.h
 *
 *  Basics formulas for the different conics
 *
 *
 *  Created on: Jul 23, 2024
 *      Author: macec
 */

#include "orbit_dynamics.h"
#include <cmath>

double ellipse_period_from_radios(double mu_center, double R1, double R2){
	// Computes the period of an ellipse from radios of periapsis and apoapsis.
	double tc;
	tc=(2*M_PI/sqrt(mu_center))*(sqrt((R1+R2)*(R1+R2)*(R1+R2))/sqrt(8));
	return tc;
}

double elliptic_velocity(double semimajor_axis, double r_distance, double mu_central){
	double vel_elliptic;
	vel_elliptic=sqrt(mu_central*((2/r_distance)-(1/semimajor_axis)));
	return vel_elliptic;
}

double escape_vel(double mu_center, double r_distance){
	/* Escape Velocity -
	 *  At a given distance r from mu, the escape velocity is
	 *  v_esc=sqrt((2mu/r))
	 */
	double v_esc;
	return v_esc=sqrt(2*mu_center/r_distance);
}

double hyperbolic_velocity(double semimajor_axis, double r_distance, double mu_central){
	double hyp_vel;
	hyp_vel=sqrt(mu_central*((2/r_distance)+(1/fabs(semimajor_axis))));
	return hyp_vel;
}

double hyperbolic_excess_velocity(double mu_center, double semimajor_axis){
	// The speed at withc a body on a hyperbolic paht arrives at infinite
	double v_inf;
	v_inf=sqrt(mu_center/semimajor_axis);
	return v_inf;
}

double hyperbolic_escape_velocity(double v_escape, double v_infinity){
	/* The hyperbolic excess speed v_inf represents the excess kinetic energy
	 * over that which is required to simply escape from the center of attraction.
	 */
	double hyper_escape_vel;
	hyper_escape_vel=sqrt(v_escape*v_escape+v_infinity*v_infinity);
	return  hyper_escape_vel;
}

double hyperbolic_eccentricity_from_v_inf(double rp, double v_inf, double mu_center){
	double hyp_ecc;
	hyp_ecc=1+rp*v_inf*v_inf/mu_center;
	return hyp_ecc;
}

double hyperolic_semimajor_axis_from_v_inf(double mu_center, double v_inf){
	double hyp_a;
	hyp_a=mu_center/(v_inf*v_inf);
	return hyp_a;
}


double C3(double v_inf){
	/* C3 Characteristic energy
	 * The energy required for an interplanetary mission.
	 * Also the maximum energy a launch vehivle can impart
	 * to a SC of a given mass. C3)launch > C3)SC
	 */
	double C3=v_inf*v_inf;
	return C3;
}




