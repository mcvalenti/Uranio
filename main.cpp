#include <filesystem>
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
//#include "interplanetary.h"
#include "orbit_dynamics.h"

using namespace std;

float distance_residual (const LDVector sv, float const distance){
	return sqrt(sv[0]*sv[0]+sv[1]*sv[1]+sv[2]*sv[2]);
}

int main() {

    // ---- OUTPUT FILES PATH -------------------    
    auto exec_path = std::filesystem::current_path();
    std::cout << "Running from: " << exec_path << std::endl;

    // Up 1 level if in build
    auto project_root = exec_path;
    if (exec_path.filename() == "build")
        project_root = exec_path.parent_path();
    // ---- ------------------- ------------------   

	// SET INIT CONDITIONS
	// ---- ------------------- ------------------   
	// Curtis Esample - page 323 and 369
	// ---- ------------------- ------------------   

	unsigned size=7;
	double Rt=6378; // Earth Radius
	double hp=480; // Perigee altitude
	double ha=800; // Apogee altitude
	double rp=Rt+hp; // Perigee radius
	double ra=Rt+ha; // Apogee radius
	double mass = 2000.0; // mass [kg]
	long double semimajorAxis=(rp+ra)/2;
	long double sv_m[size]={rp,0,0,0,7.7102,0,mass};

	LDVector init_sv(sv_m, 7);
	//double aux=get_max_absolute(init_sv); ??????


	//-----------------------------------------
	// Orbit 1 - Parking Orbit - Central Body
	//-----------------------------------------
	std::cout <<std::endl;
	std::cout <<std::endl;
	std::cout <<"---------------------"<<std::endl;
	std::cout <<"- PARKING ORBIT - "<<std::endl;
	std::cout <<"---------------------"<<std::endl;
	float step=1.0;
	double period0=period(semimajorAxis);
	int total_time=ceil(period0);
	string filename0=(project_root / "output_files/orbit1.csv").string();

	// First propagation instance
	propagator parking_orbit(init_sv,total_time,step);
	cBody_param cbody;
	cbody.mu=398600.448;
	LDVector cbody_last_sv;
	parking_orbit.addPerturbation(&central_body, &cbody); // args: function and structure
	parking_orbit.propagate(filename0);
	cbody_last_sv=parking_orbit.last_sv;
	parking_orbit.sv2oe(parking_orbit.last_sv);
	// TODO: create a method that print orbit keplerian elements
	std::cout <<"Semimajor Axis : "<<parking_orbit.a <<std::endl;
	std::cout <<"True Anomaly : "<<parking_orbit.nu <<std::endl;


	//-----------------------------------------
	// Orbit 2 - TRANSFER ORBIT
	// With Continuous thrust
	//-----------------------------------------
	std::cout <<"---------------------"<<std::endl;
	std::cout <<"- TRANSFER ORBIT "<<std::endl;
	std::cout <<"---------------------"<<std::endl;
	float thrust_time=261.0; //1127; // NOTE! Only represents the complete time if step=1.0
	float thrust_step=1.0;
	string filename1=(project_root / "output_files/orbit2.csv").string();
	propagator transfer_orbit(init_sv,thrust_time,thrust_step);
	thrust_param tparam;
	tparam.isp = 300.0;
	tparam.thrust = 10000.0;
	transfer_orbit.addPerturbation(&central_body, &cbody);
	transfer_orbit.addPerturbation(&thrust, &tparam);
	transfer_orbit.propagate(filename1);
	transfer_orbit.sv2oe(transfer_orbit.last_sv);
	std::cout <<"End-of-burn state vector : "<<transfer_orbit.last_sv <<std::endl;
	std::cout <<"Semimajor Axis : "<<transfer_orbit.a <<std::endl;
	std::cout <<"Eccentricity : "<<transfer_orbit.e <<std::endl;
	std::cout <<"True anomaly nu: "<<transfer_orbit.nu <<std::endl;
	double delta_nu=(M_PI-transfer_orbit.nu)*180.0/M_PI; // computes delta true anomaly to apogee
	std::cout <<"Delta anomaly nu: "<<delta_nu <<std::endl;
	// State Vector in Apogee
	LDVector apogee_sv = sv_from_true_anomaly(transfer_orbit.last_sv,delta_nu); //to compute radius at final point
	long double apogee_with_mass[7]={apogee_sv[0], apogee_sv[1],apogee_sv[2],apogee_sv[3],apogee_sv[4],apogee_sv[5],transfer_orbit.last_sv[6]} ;
	LDVector apogee_sv1(apogee_with_mass, 7);
	std::cout <<"Apogee state vector: "<<apogee_sv1 <<std::endl;

	std::cout <<"---------------------"<<std::endl;
	std::cout <<"- Propagation without thrusting "<<std::endl;
	std::cout <<"---------------------"<<std::endl;
	// TODO - PROPAGATE TO APOGEE how to get time to propagate ?????

	std::cout <<"---------------------"<<std::endl;
	std::cout <<"- Circularization Orbit "<<std::endl;
	std::cout <<"---------------------"<<std::endl;	
	float thrust_time_circ=118.88; 
	float thrust_step_circ=1.0;
	string filename2=(project_root / "output_files/orbit3.csv").string();
	propagator circularization_orbit(apogee_sv1,thrust_time_circ,thrust_step_circ);
	circularization_orbit.addPerturbation(&central_body, &cbody);
	circularization_orbit.addPerturbation(&thrust, &tparam);
	circularization_orbit.propagate(filename2);
	circularization_orbit.sv2oe(circularization_orbit.last_sv);
	std::cout <<"End-of-burn state vector : "<<circularization_orbit.last_sv <<std::endl;
	std::cout <<"Semimajor Axis : "<<circularization_orbit.a <<std::endl;

	std::cout <<"Eccentricity : "<<circularization_orbit.e <<std::endl;
	std::cout <<"True anomaly nu: "<<circularization_orbit.nu <<std::endl;

	std::cout <<"---------------------"<<std::endl;
	std::cout <<"- Final Orbit "<<std::endl;
	std::cout <<"---------------------"<<std::endl;
	float final_time=20000.88; 
	float final_step=1.0;
	string filename4=(project_root / "output_files/orbit4.csv").string();
	propagator final_orbit(circularization_orbit.last_sv, final_time, final_step);
	final_orbit.addPerturbation(&central_body, &cbody); // args: function and structure
	final_orbit.propagate(filename4);

	cout << endl;
	cout<<"End of processing!";
	
	return 0;

}
