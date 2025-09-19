/*
 * interplanetary.cpp
 *  [Ref 1] Curtis, H. D. (2020). Orbital mechanics for engineering students:
 *	 		Revised Reprint. Butterworth-Heinemann.
 *  Plantery departure: concept of Hyperbolic Excess Velocity = V_inf
 *  Planetary Rendezvous
 *
 *  Created on: Jul 13, 2024
 *      Author: macec
 */
#include <iostream>
#include <cmath>
#include "interplanetary.h"
#include "orbit_dynamics.h"

using namespace std;

// TO DO: compute phase_end according to inter_to_outer // functions

void run_planetary_departure(){
	// variables
	double v_inf, v_elliptic, delta_v_final, v_park;
	double vel_escape, hyper_escape, v_eclip_norm;
	double dv_total_departure, C3_orbital, ecc;
	double beta, hyp_a;
	// Constants
	double pi = M_PI;
	//double AU = 150e6; //  km 149597870.7;
	double mu_earth = 398604; // km3/s2
	double r_earth=6378;
	double h_parking=1000; // km
	double r_parking=h_parking+r_earth; // km
	double mu_sun=1.33e11; // km3/s2

	// Delta V Heliocentric Hohmann Transfer
	double R1=148.7e06; // km
	double R2=39.784e06; // km venus
	double inc=22.5; // deg inclination wrt Ecliptic plane

	std::cout<<"================================================="<<std::endl;
	v_inf=delta_V_hohmann_heliocentric(mu_sun, R1, R2);
	std::cout<<"Delta V Hohmann Heliocentric - V_inf: "<<v_inf<<" km/s"<<std::endl;
	std::cout<<"================================================="<<std::endl;
	v_elliptic=elliptic_velocity((R1+R2)/2, R1, mu_sun);
	std::cout<<"S/C velocity in Elliptic Heliocentric - DV1: "<<v_elliptic<<std::endl;
	std::cout<<"================================================="<<std::endl;
	// Projection In the inclination plane
	v_eclip_norm=v_elliptic*tan(inc*pi/180.0);
	std::cout<<" V_eclip_norm : "<<v_eclip_norm<<std::endl;
	std::cout<<"================================================="<<std::endl;
	// Total DV departure (module)
	dv_total_departure=sqrt(v_inf*v_inf+v_eclip_norm*v_eclip_norm);
	std::cout<<"Total DV departure in inclined plane : "<<dv_total_departure<<std::endl;
	C3_orbital=dv_total_departure*dv_total_departure;
	std::cout<<" C3 (orbital): "<<C3_orbital<<std::endl;
	std::cout<<"================================================="<<std::endl;
	v_park=parking_v(mu_earth, r_parking);
	delta_v_final=delta_v_to_hyperbola(v_park, v_inf);
	std::cout<<"V Parking : "<<v_park<<std::endl;
	std::cout<<"From Parking to Hyperbolic excess - Total DV: "<<delta_v_final<<std::endl;
	std::cout<<"================================================="<<std::endl;
	vel_escape= escape_vel(mu_earth, r_parking);
	std::cout<<"V escape from parking orbit - v_escape: "<<vel_escape<<std::endl;
	hyper_escape=hyperbolic_escape_velocity(vel_escape, v_inf);
	std::cout<<"Total departure velocity - sqrt(v_inf2+v_esc2)- Total V: "<<hyper_escape<<std::endl;
	std::cout<<"================================================="<<std::endl;
	ecc= hyperbolic_eccentricity_from_v_inf(r_parking, v_inf, mu_earth);
	std::cout<<"Hyperbolic ecc: "<<ecc<<std::endl;
	std::cout<<"================================================="<<std::endl;
	beta= beta_angle(ecc);
	std::cout<<"Beta : "<<beta*180.0/M_PI<<std::endl;
	std::cout<<"================================================="<<std::endl;
	std::cout<<"R parking : "<<r_parking<<std::endl;
	std::cout<<"================================================="<<std::endl;
	hyp_a=hyperolic_semimajor_axis_from_v_inf(mu_earth, v_inf);
	std::cout<<"Hyperbola semimajor axis : "<<hyp_a<<std::endl;
}

void run_rendezvous_opportunities(){
	/* Chapter 8 - Section 3 (Curtis)
	 * Computes the minimum wait time to return from Mars to Earth
	 * Example 8.2 - [pag. 435]
	*/

	//Constants
	double R_departure = 227.9e06; // Km Mars Distance to Sun
	double R_arrival = 149.6e06;   // km Earth's Distance to Sun
	double mu_sun = 132.71e09;     // Sun mu
	double n_departure = 0.01720;    // rad/day (2*pi/365.26)
	double n_arrival = 0.0091327; // rad/day (2*pi/687.99)

	// Variables
	double t_12, t_12_days, phase_end, t_wait;
	t_12=ellipse_period_from_radios(mu_sun, R_departure, R_arrival); // sec
	t_12_days = (t_12/2)/86400; // days
	phase_end=M_PI-n_departure*t_12_days; // rad
	t_wait=time_wait(phase_end, n_departure, n_arrival); // days
	std::cout<<"================================================="<<std::endl;
	std::cout<<"Time to wait t_wait: "<<t_wait<<std::endl;
	std::cout<<"================================================="<<std::endl;
}

void run_venus_rendezvous(){
	/*
	 *  Venus rendezvous
	 */
	// Variables (TO DO list and describe)
	double t_12, phi_end;
	double v_inf, hyp_a, v_hyp, t_soi;
	double h_hyp, theta_hyp, F, ecc, M_hyp;

	// Constants
	//double n_earth= 0.01720; // [rad/day]
	double n_venus=0.0279;	// [rad/day]
	double R_earth = 149.6e06;   // km Earth's Distance to Sun
	double R_venus = 107.5e06;   // km Venus's Distance to Sun
	double mu_sun = 132.71e09;     // Sun mu
	double mu_earth = 398604; // km3/s2
	double r_earth=6378;
	double h_parking=300; // km
	double r_parking=h_parking+r_earth; // km
	double r_soi= 9.24e05; // km

	t_12=ellipse_period_from_radios(mu_sun, R_earth, R_venus)/(2*86400.0); // sec
	phi_end=n_venus*t_12-M_PI;
	std::cout<<"================================================="<<std::endl;
	std::cout<<"Time to Venus: "<<t_12<<std::endl;
	std::cout<<"================================================="<<std::endl;
	std::cout<<"Phase with Earth: "<<phi_end*180.0/M_PI<<std::endl;
	std::cout<<"================================================="<<std::endl;
	v_inf=delta_V_hohmann_heliocentric(mu_sun, R_earth, R_venus);
	std::cout<<"================================================="<<std::endl;
	ecc= hyperbolic_eccentricity_from_v_inf(r_parking, v_inf, mu_earth);
	std::cout<<"Hyperbolic ecc: "<<ecc<<std::endl;
	std::cout<<"================================================="<<std::endl;
	hyp_a=hyperolic_semimajor_axis_from_v_inf(mu_earth, v_inf);
	std::cout<<"Hyperbola semimajor axis : "<<hyp_a<<std::endl;


	// Computation sequence
	// 1- Compute angular momentum - h=rp*vp
	// 2- Compute true anomaly at SOI
	// 3- Compute Eccentric anomaly for hyperbola F
	// 4- Compute Mean anomaly from F
	// 5- Compute Time from periapsis when arriving soi

	std::cout<<"================================================="<<std::endl;
	v_hyp = hyperbolic_velocity(hyp_a, r_parking, mu_earth);
	std::cout<<"Hyperbolic velocity : "<<v_hyp<<std::endl;
	std::cout<<"================================================="<<std::endl;
	h_hyp = r_parking*v_hyp;
	std::cout<<"Angular Momentum of the Hyperbola : "<<h_hyp<<std::endl;
	std::cout<<"================================================="<<std::endl;
	theta_hyp = acos((((h_hyp*h_hyp)/(r_soi*mu_earth))-1)/ecc);
	std::cout<<"Theta : "<<theta_hyp*180.0/M_PI<<std::endl;
	std::cout<<"================================================="<<std::endl;
	F=2*atanh(sqrt((ecc-1)/(ecc+1))*tan(theta_hyp/2));
	std::cout<<"F : "<<F*180.0/M_PI<<std::endl;
	std::cout<<"================================================="<<std::endl;
	M_hyp=ecc*sinh(F)-F;
	std::cout<<"M_hyp : "<<M_hyp*180.0/M_PI<<std::endl;
	t_soi=h_hyp*h_hyp*h_hyp*M_hyp/(mu_earth*mu_earth*(sqrt((ecc*ecc-1)*(ecc*ecc-1)*(ecc*ecc-1))));
	std::cout<<"================================================="<<std::endl;
	std::cout<<"Time to SOI : "<<t_soi/86400<<std::endl;
}

double synodic_period_from_velocities(double n1, double n2){
	/*
	 * Considering two orbiting bodies, the time required for the phase
	 * angle to return to its initial value value is called the Synodic Period (t_syn).
	 * If n1>n1 is the time required for phi to change  from phi0 to phi0-2*PI.
	 * If n2>n1 is the time required for phi to change from phi0 to phi0+s*PI.
	 * t_syn:[days] Earth days 24hs
	  */
	double t_syn;
	if (n1>n2){
		t_syn=2*M_PI/(n1-n2);
	} else {
		t_syn=2*M_PI/(n2-n1);
	}
	return t_syn;
}

double synodic_period_from_periods(double t1, double t2){
	/*
	 * Synodic Period (t_syn), is the orbital period of planet 2
	 * relative to planet 1.
	 * t1, t2, t_syn:[days] Earth days 24hs
	 */
	double t_syn;
	if (t1>t2){
		t_syn=t1*t2/(t1-t2);
	} else {
		t_syn=t1*t2/(t2-t1);
	}
	return t_syn;
}

void planetary_rendezvous(char inner_outer){
	// Planetary rendezvous & flyby
	//double v_inf, v_sc, v_planet;
	//seguir aca

}

double delta_V_hohmann_heliocentric(double mu_sun, double R1, double R2){
	/* Planetary departure
	* Heliocentric velocity for a SC departing on a Hohmann trajectory,
	* to a point farther from the Sun.
	* R1: Distance of the departure planet (planet 1) to the Sun
	* R2: Distance of the arriving point/planet (planet 2) to the Sun
	*/
	double delta_v;
	delta_v=sqrt(mu_sun/R1)*(sqrt(2*R2/(R1+R2))-1);

	return delta_v;
}


double parking_v(double mu_center, double r_departure){
	double parking_v;
	parking_v=sqrt(mu_center/r_departure);

	return parking_v;
}

double delta_v_to_hyperbola(double parking_v, double delta_v_inf){
	double delta_v_hyp;
	delta_v_hyp=parking_v*(sqrt(2+(delta_v_inf*delta_v_inf/(parking_v*parking_v)))-1);

	return delta_v_hyp;
}

double beta_angle(double ecc){
	// Orientation of the apse line of the hyperbola
	// to the planet's heliocentric velocity vector.
	double beta;
	beta=acos(1/ecc);
	return beta;
}

double time_wait(double phase_end, double n1, double n2){
	/*
	 * n1: Departure Planet Angular velocity
	 * n2: Arrival Planet Angular velocity
	 */
	double t_wait=-1;
	int N=0;
	if (n1>n2){
		while (t_wait < 0){
			t_wait=(-2*phase_end-2*M_PI*N)/(n2-n1);
			N=N+1;
		}
	}else{
		while (t_wait < 0){
			t_wait=(-2*phase_end+2*M_PI*N)/(n2-n1);
			N=N+1;
			}
		}
	return t_wait;
};

double return_trip(double mu_center, double R_arrival, double R_departure, double n_arrival, double n_departure){
	/*
	 * Minimum time to return from p1 to p2
	 *
	 */
	double t_wait, t_12, phase_end;

	t_12=ellipse_period_from_radios(mu_center, R_arrival, R_departure);
	phase_end= M_PI-n_arrival*t_12;
	t_wait=time_wait(phase_end, n_departure, n_arrival);

	return t_wait;
}



// Extra data
//double falconH=11.2; // km/s
//double distance_mars_sun = 1.524*AU;
//double distance_venus_sun = 1.077504e8; // 1.077504 x 10^8 km
//double distance_mercury_sun = 5.833821e7; // 5.833821 x 10^7 km
//double mass_sun = 1.989e30; // Masa del Sol en kg
//double G = 6.67430e-11;  // m^3 kg^-1 s^-2
// Flybys data
// Venus_SOI=6.14e05 // km
// Venus_SOI_bd=100
// Earth_SOI=9.24e05 // km
// Earth_SOI_bd=145
//te=365.26; // Earth's period days
//tm=687.99; // Mars's period days
//period t_earth=365.26;
//period t_venus=225;
//period t_mercury=88;
// Venus radio = 6051.8; /km
