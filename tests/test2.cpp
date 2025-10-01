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
	double Rt=6378; // Earth Radius
	double hp=480; // Perigee altitude
	double ha=800; // Apogee altitude
	double rp=Rt+hp; // Perigee radius
	double ra=Rt+ha; // Apogee radius
	long double semimajorAxis=(rp+ra)/2;
	long double sv_m[size]={6858,0,0,0,7.7102,0,2000};

	LDVector init_sv(sv_m, 7);
	//double aux=get_max_absolute(init_sv); ??????

	//-----------------------------------------
	// Orbit 2 - With Continuous thrust
	//-----------------------------------------
	
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
	cout << "At apogee: " << final << endl;

	cout<<"End of processing!"<<std::endl;
	
	system("pause"); // Wait until press a key. 
	return 0;
}













