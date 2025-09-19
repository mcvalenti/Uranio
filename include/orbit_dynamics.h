/*
 * orbit_dynamics.h
 *
 *  Basics formulas for the different conics
 *  Description in the .cpp file
 *
 *  Created on: Jul 23, 2024
 *      Author: macec
 */

#ifndef ORBIT_DYNAMICS_H_
#define ORBIT_DYNAMICS_H_

double elliptic_velocity(double semimajor_axis, double r_distance, double mu_central);
double ellipse_period_from_radios(double mu_center, double R1, double R2);
double escape_vel(double mu_center, double r_distance);
double hyperbolic_eccentricity_from_v_inf(double rp, double v_inf, double mu_center);
double hyperolic_semimajor_axis_from_v_inf(double mu_center, double v_inf);
double hyperbolic_velocity(double semimajor_axis, double r_distance, double mu_central);
double hyperbolic_excess_speed(double mu_center, double semimajor_axis);
double hyperbolic_escape_velocity(double v_escape, double v_infinity);
double C3(double v_inf);




#endif /* ORBIT_DYNAMICS_H_ */
