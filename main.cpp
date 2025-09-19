#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <time.h>
#include "propagators.h"
#include "auxiliaries.h"
#include "LDVector.h"
#include "my_tests.h"
#include "mission.h"
#include "interplanetary.h"
#include "orbit_dynamics.h"

using namespace std;


int main() {

	unsigned size=7;
	float mass;
	long double sv_m[size]={6858,0,0,0,7.7102,0,2000};

	LDVector init_sv(sv_m, 7);
	mass=init_sv[6];
	SpaceVehicle my_sat(init_sv, mass);
	// TESTS
	/*	cout << "--------TESTs-------------\n";
	run_tests();
	cout << "-------- END OF TESTs-------------\n";*/

	//run_venus_rendezvous();
	//run_planetary_departure();
	double mu_central=1.33e11;
	double r_distance=39800000;
	double R1=39800000;
	double R2=155358642;
	double semimajor_axis=(R1+R2)/2;
	double semimajor_axis2=(R1+149000000)/2;
	double v1, v2, dv;
	v1=elliptic_velocity(semimajor_axis, r_distance, mu_central);
	v2=elliptic_velocity(semimajor_axis2, r_distance, mu_central);
	dv=v1-v2;
	cout << "-------- DV : "<<dv<<endl;
	cout << "-------- Semimajor a : "<< semimajor_axis2 <<endl;


	//double Rt=6378;
	//double hp=480; // Perigee altitude
	//double ha=800; // Apogee altitude
	//double rp=Rt+hp; // Perigee radius
	//double ra=Rt+ha; // Apogee radius
	//long double semimajorAxis=(rp+ra)/2;
	//double aux=get_max_absolute(init_sv); ????????????????????????????????
	//cout << aux << endl;




	//-----------------------------------------
	// Orbit 2 - With Continuous thrust
	//-----------------------------------------
	/*
	std::cout <<"---------------------"<<std::endl;
	std::cout <<"- CONTINUOUS THRUST "<<std::endl;
	std::cout <<"---------------------"<<std::endl;
	float thrust_time=261.1127;
	float thrust_step=0.1;
	string filename1="output_files/orbit2.csv";
	cBody_param cbody;
	cbody.mu=398600.448;
	propagator thrust_prop(init_sv,thrust_time,thrust_step);
	thrust_param tparam;
	tparam.isp = 300.0;
	tparam.thrust = 10000.0;
	thrust_prop.addPerturbation(&central_body, &cbody);
	thrust_prop.addPerturbation(&thrust, &tparam);
	thrust_prop.propagate(filename1);

	cout<<" End-of-burn state vector: "<<thrust_prop.last_sv<<endl;
	thrust_prop.sv2oe(thrust_prop.last_sv);
	std::cout <<"True anomaly nu: "<<thrust_prop.nu <<std::endl;
	double delta_nu=(M_PI-thrust_prop.nu)*180.0/M_PI; // computes delta true anomaly to apogee
	std::cout <<"Delta anomaly nu: "<<delta_nu <<std::endl;
	LDVector final = sv_from_true_anomaly(thrust_prop.last_sv,delta_nu); //to compute radius at final point
	cout << "At apogee: " << final << endl;*/



	/*

	//-----------------------------------------
	// Orbit 1 - Initial - Central Body
	//-----------------------------------------
	std::cout <<std::endl;
	std::cout <<std::endl;
	std::cout <<"---------------------"<<std::endl;
	std::cout <<"- CENTRAL BODY "<<std::endl;
	std::cout <<"---------------------"<<std::endl;
	float step=1.0;
	double period0=period(semimajorAxis);
	int total_time=ceil(period0);
	string filename0="output_files/orbit1.csv";

	// First propagation instance
	propagator cbody_prop(init_sv,total_time,step);
	cBody_param cbody;
	cbody.mu=398600.448;
	LDVector cbody_last_sv;

	// Add perturbation
	cbody_prop.addPerturbation(&central_body, &cbody); // args: function and structure
	const clock_t begin_time = clock();
	cbody_prop.propagate(filename0);
	cbody_last_sv=cbody_prop.last_sv;

	std::cout <<cbody_prop.last_sv <<std::endl;
	cbody_prop.sv2oe(cbody_prop.last_sv);
	std::cout <<"Semimajor Axis : "<<cbody_prop.a <<std::endl;


	//-----------------------------------------
	// Orbit 2 - With Continuous thrust
	//-----------------------------------------
	std::cout <<"---------------------"<<std::endl;
	std::cout <<"- CONTINUOUS THRUST "<<std::endl;
	std::cout <<"---------------------"<<std::endl;
	float thrust_time=261.1127;
	float thrust_step=0.1;
	string filename1="output_files/orbit2.csv";
	propagator thrust_prop(init_sv,thrust_time,thrust_step);
	thrust_param tparam;
	tparam.isp = 300.0;
	tparam.thrust = 10000.0;
	thrust_prop.addPerturbation(&central_body, &cbody);
	thrust_prop.addPerturbation(&thrust, &tparam);
	thrust_prop.propagate(filename1);


	cout<<" End-of-burn state vector: "<<thrust_prop.last_sv<<endl;
	thrust_prop.sv2oe(thrust_prop.last_sv);
	std::cout <<"True anomaly nu: "<<thrust_prop.nu <<std::endl;
	double delta_nu=(M_PI-thrust_prop.nu)*180.0/M_PI; // computes delta true anomaly to apogee
	std::cout <<"Delta anomaly nu: "<<delta_nu <<std::endl;
	LDVector final = sv_from_true_anomaly(thrust_prop.last_sv,delta_nu); //to compute radius at final point
	cout << "At apogee: " << final << endl;

	// if not ra=22378 --> increase thrust_time

	//Estimation
	long double semimajor_axis_transfer=14618.0;
	double period2=period(semimajor_axis_transfer);
	// transfer Free propagation
	propagator free_prop(thrust_prop.last_sv,period2-thrust_time,step);
	cBody_param free_body;
	free_body.mu=398600.448;
	LDVector free_body_last_sv;

	// Add perturbation
	free_prop.addPerturbation(&central_body, &free_body); // args: function and structure
	string filename3="output_files/orbit3.csv";
	free_prop.propagate(filename3);
	free_prop.sv2oe(free_prop.last_sv);
	std::cout <<"Semimajor Axis : "<<free_prop.a <<std::endl;

	// End of Propagation
	std::cout << "Propagation time: "<<float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;





	// End of Code
	std::cout << "Propagation time: "<<float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
	cout << endl;
	cout<<"End of processing!";
	*/

	return 0;
}













