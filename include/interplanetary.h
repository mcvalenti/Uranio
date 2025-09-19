/*
 * interplanetary.h
 *
 *  Created on: Jul 13, 2024
 *      Author: macec
 */

#ifndef INTERPLANETARY_H_
#define INTERPLANETARY_H_

void run_planetary_departure();
void run_rendezvous_opportunities();
void run_venus_rendezvous();
double delta_V_hohmann_heliocentric(double mu_sun, double R1, double R2);
double delta_v_to_hyperbola(double parking_v, double delta_v_inf);
double parking_v(double mu_center, double r_departure);
double beta_angle(double ecc);
double synodic_period_from_velocities(double n1, double n2);
double synodic_period_from_periods(double t1, double t2);
void planetary_rendezvous();
double time_wait(double phase_end, double n1, double n2);
double return_trip(double mu_center, double R_arrival, double R_departure);





#endif /* INTERPLANETARY_H_ */
