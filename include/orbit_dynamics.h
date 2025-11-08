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
double eccentric_from_true_anomaly(double eccentricity, double true_anomaly);
double mean_from_eccentric(double eccentric, double eccentricity);
double kepler(double M, double e, double tol = 1e-10, int max_iter = 50);
double mean_motion(double mu_center, double semimajor_axis);
double delta_time_from_mean_anomaly(double delta_mean_anomaly, double mean_motion, double semimajor_axis);
double time_between_two_true_anomalies(double t_anomaly1, double t_anomaly2, double ecc, double sma, double mu_center);
double escape_vel(double mu_center, double r_distance);
double hyperbolic_eccentricity_from_v_inf(double rp, double v_inf, double mu_center);
double hyperolic_semimajor_axis_from_v_inf(double mu_center, double v_inf);
double hyperbolic_velocity(double semimajor_axis, double r_distance, double mu_central);
double hyperbolic_excess_speed(double mu_center, double semimajor_axis);
double hyperbolic_escape_velocity(double v_escape, double v_infinity);
double ellipse_period_from_semimajor_axis(double mu, double a);
double C3(double v_inf);





#endif /* ORBIT_DYNAMICS_H_ */
