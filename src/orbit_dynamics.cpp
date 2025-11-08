/*
 * orbit_dynamics.h
 *
 *  Basics formulas for the different conics
 *
 *
 *  Created on: Jul 23, 2024
 *      Author: macec
 */
#include <iostream>
#include "orbit_dynamics.h"
#include <cmath>
#include <stdexcept>

using namespace std;

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

double eccentric_from_true_anomaly(double eccentricity, double true_anomaly){
	/* Computes eccentric anomaly from true anomaly */
	double eccentric;
	double eccentric_parabola_o2;
	// Ellipse
	// (Curtis, 4th ed., p. 151–152)
	if (eccentricity < 1){
		double denom = 1.0 + eccentricity *cos(true_anomaly);
		if (denom != 0){
			double cosE = (eccentricity + cos(true_anomaly)) / denom;
			double sinE = sqrt( 1.0 - eccentricity*eccentricity)*sin(true_anomaly) / denom;
			eccentric = atan2(sinE, cosE);
		}
	}
	// Hyperbola
	// (Curtis, 4th ed., p. 166)
	if (eccentricity > 1){
		double denom = 1.0 + eccentricity*cos(true_anomaly);
		if (denom != 0){
			double coshF = (eccentricity + cos(true_anomaly)) / denom;
			double sinhF = sqrt(eccentricity*eccentricity - 1.0)*sin(true_anomaly) / denom;
			if (fabs(sinhF) < 1e-10){
				eccentric = 0;
			}
			else {
				eccentric = asinh(sinhF);
			}
		}	
	}
	// parabola
	if (eccentricity == 1){
		eccentric_parabola_o2  = tan(true_anomaly/2.0);
		eccentric = 2*eccentric_parabola_o2;
	}
	return eccentric;
}

double mean_from_eccentric(double eccentric, double eccentricity){
	double mean_anomaly;
	// Ellipse
	if (eccentricity < 1){
		mean_anomaly = eccentric - eccentricity * sin(eccentric);
	}
	// Hyperbola
	if (eccentricity > 1){
		mean_anomaly = eccentricity* sinh(eccentric) - eccentric;
	}
	return mean_anomaly;
}

double kepler(double M, double e, double tol, int max_iter) {
    double E = M;
    for (int i = 0; i < max_iter; ++i) {
        double f = E - e * sin(E) - M;
        double fp = 1.0 - e * cos(E);
        double dE = -f / fp;
        E += dE;
        if (fabs(dE) < tol) break;
    }
    return E;
}

double mean_motion(double mu_center, double semimajor_axis){
	double mean_motion;
	double sma;
	// Ellipse
	if (semimajor_axis > 0){
		mean_motion = sqrt(mu_center/(semimajor_axis*semimajor_axis*semimajor_axis));
	}
	// Hyperbola
	if (semimajor_axis < 0){
		sma = abs(semimajor_axis);
		mean_motion = sqrt(mu_center/(sma*sma*sma));
	}
	return mean_motion;

}

double delta_time_from_mean_anomaly(double delta_mean_anomaly, double mean_motion, double semimajor_axis){
	double delta_time;
	// Ellipse
	if (semimajor_axis > 0){
		// TO DO check multiples of 2 pi.
		delta_time = delta_mean_anomaly/mean_motion;
	}
	// Hyperbola
	if (semimajor_axis < 0){
		delta_time = delta_mean_anomaly/mean_motion;
	}
	// Parabola . Barker formula
	// t - T = 1/2 sqrt((p*p*p)/mu_center) (D+ 1/3 D*D*D);	
	
	return delta_time;
}

double time_between_two_true_anomalies(double t_anomaly1, double t_anomaly2, double ecc, double sma, double mu_center){
	double E1, E2;
	double M1, M2, deltaM;
	double n;
	double delta_t;
	E1 = eccentric_from_true_anomaly(ecc, t_anomaly1);
	E2 = eccentric_from_true_anomaly(ecc, t_anomaly2);
	M1 = mean_from_eccentric(E1, ecc);
	M2 = mean_from_eccentric(E2, ecc);
	n = mean_motion(mu_center, sma);
	deltaM = M2-M1;
	delta_t = delta_time_from_mean_anomaly(deltaM, n, sma);
	return delta_t;
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

double ellipse_period_from_semimajor_axis(double mu, double a){
	// Computes period from Kepler's 3rd law. 
	double period;
	period = 2*M_PI*sqrt(a*a*a/mu);

 // Validation
	if (std::isnan(period) || std::isinf(period)) {
		throw std::runtime_error(" period Not valid");
	}
	return period;
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




