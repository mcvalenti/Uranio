#include <filesystem>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <time.h>
#include "constants.h"
#include "propagators.h"
#include "auxiliaries.h"
#include "LDVector.h"
#include "my_tests.h"
#include "mission.h"
//#include "interplanetary.h"
#include "orbit_dynamics.h"
#include <optional>

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
	double Rt=6378; 		// Earth Radius [km]
	double hp=480; 			// Perigee altitude [km]
	double ha=800; 			// Apogee altitude [km]
	double rp=Rt+hp; 		// Perigee radius [km]
	double ra=Rt+ha; 		// Apogee radius [km]
	double mass = 2000.0; 	// mass [kg]
	long double semimajorAxis=(rp+ra)/2;
	long double sv_m[size]={rp*1000,0,0,0,7.7102*1000,0,mass};

	LDVector init_sv(sv_m, 7);
	//double aux=get_max_absolute(init_sv); ??????


	//-----------------------------------------
	// Orbit 0 - Parking Orbit - Central Body
	//-----------------------------------------
	std::cout <<std::endl;
	std::cout <<std::endl;
	std::cout <<"---------------------"	<<std::endl;
	std::cout <<"- PARKING ORBIT - "	<<std::endl;
	std::cout <<"---------------------"	<<std::endl;
	float step=1.0;
	double period0=period(semimajorAxis);
	int total_time=ceil(period0);
	string filename0=(project_root / "output_files/orbit0.csv").string();

	// First propagation instance
	propagator parking_orbit(init_sv,total_time,step);
	cBody_param cbody;
	cbody.mu=GM; // [m*3/s*2]
	LDVector cbody_last_sv;
	parking_orbit.addPerturbation(&central_body, &cbody); // args: function and structure
	parking_orbit.propagate();
	parking_orbit.ephemeris_tofile(filename0);
	cbody_last_sv=parking_orbit.last_sv;
	parking_orbit.sv2oe(parking_orbit.last_sv);
	// TODO: create a method that print orbit keplerian elements
	std::cout <<"Semimajor Axis : "<<parking_orbit.a <<std::endl;
	std::cout <<"True Anomaly : "<<parking_orbit.nu <<std::endl;

	//-----------------------------------------
	// Orbit 1 - TRANSFER ORBIT
	// With Continuous thrust
	//-----------------------------------------
	std::cout <<"---------------------"	<<std::endl;
	std::cout <<"- TRANSFER ORBIT "		<<std::endl;
	std::cout <<"---------------------"	<<std::endl;
	float thrust_time=255.0; 
	float thrust_step=0.5;
	thrust_param tparam;
	tparam.isp = 300.0;
	tparam.thrust = 10000.0; // [kN]
	string filename1=(project_root / "output_files/orbit1.csv").string();
	
	double delta_nu;
	LDVector apogee_sv1;
	
	// Loop variables init
	double res_ra_min = 1000000.0; // res_ra_min: min residual to target apogee radius  [m]
	double res_ra = 0.0;		 // res_ra: residual to target apogee radius 
	double tol = 50000.0;		 // tolerance [m]
	double thrust_time_min = 0;	 // thurst_time_min: minimum propulsion interval to get to target

	propagator transfer_orbit;
	//LDVector ld(transfer_orbit.current_sv);  ;

	int iter_num=1;
	string iter_file="output_files/orbit"+to_string(iter_num)+".csv";
	while ((res_ra_min > tol && transfer_orbit.current_sv[6]> 0)){
		string filename0=(project_root / iter_file).string();
		transfer_orbit.reset_init(init_sv,thrust_time,thrust_step);
		transfer_orbit.addPerturbation(&central_body, &cbody); 
		transfer_orbit.addPerturbation(&thrust, &tparam);
		transfer_orbit.propagate();
		transfer_orbit.sv2oe(transfer_orbit.last_sv);
		//std::cout <<"End-of-burn state vector : "<<transfer_orbit.last_sv <<std::endl;
		//std::cout <<"Semimajor Axis : "<<transfer_orbit.a <<"	Eccentricity : "<<transfer_orbit.e<<std::endl;
		//std::cout <<"Eccentricity : "<<transfer_orbit.e <<std::endl;
		// std::cout <<"True anomaly nu: "<<transfer_orbit.nu <<std::endl;
		delta_nu=(180.0-transfer_orbit.nu); // computes delta true anomaly to apogee [deg]
		//std::cout <<"Delta anomaly nu: "<<delta_nu <<std::endl;
		// State Vector in Apogee
		LDVector apogee_sv = sv_from_true_anomaly(transfer_orbit.last_sv,delta_nu); //to compute radius at final point
		long double apogee_with_mass[7]={apogee_sv[0], apogee_sv[1],apogee_sv[2],apogee_sv[3],apogee_sv[4],apogee_sv[5],transfer_orbit.last_sv[6]} ;
		LDVector apogee_sv1(apogee_with_mass, 7);
		//std::cout <<"Apogee state vector: "<<apogee_sv1 <<std::endl;
		double ra_apogee = sqrt(apogee_sv1[0]*apogee_sv1[0]+apogee_sv1[1]*apogee_sv1[1]+apogee_sv1[2]*apogee_sv1[2]);
		//std::cout <<" ra apogee computed: " << ra_apogee/1000.0 << std::endl;
		std::cout <<" Iteration final mass: " << transfer_orbit.last_sv[6] << std::endl;
		double ra_target = 22378000.0; // [m]
		double res_ra = abs(ra_apogee - ra_target);
		//std::cout <<" Semimajor axis post burn: " << transfer_orbit.a << std::endl;
		if (res_ra < res_ra_min) {
			res_ra_min = res_ra;
			thrust_time_min = thrust_time;
		}
		thrust_time = thrust_time + 1;
		iter_num=iter_num + 1;
	}

	transfer_orbit.ephemeris_tofile(filename1);
	std::cout <<" Thrust time min: " << thrust_time_min << std::endl;
	//-----------------------------------------
	// Orbit 2 - Propagation without thrusting
	//-----------------------------------------
	std::cout <<"--------------------------------"	<<std::endl;
	std::cout <<"- Propagation without thrusting "	<<std::endl;
	std::cout <<"--------------------------------"	<<std::endl;
	string filename2=(project_root / "output_files/orbit2.csv").string();
	double time_to_apogee = time_between_two_true_anomalies(transfer_orbit.nu, M_PI, transfer_orbit.e, transfer_orbit.a, cbody.mu);
	propagator onTransfer_orbit(transfer_orbit.last_sv,time_to_apogee,thrust_step);
	onTransfer_orbit.addPerturbation(&central_body, &cbody);
	onTransfer_orbit.propagate();
	onTransfer_orbit.ephemeris_tofile(filename2);

	//-----------------------------------------
	// Orbit 3 - Circularization Orbit
	//-----------------------------------------
	std::cout <<"------------------------"	<<std::endl;
	std::cout <<"- Circularization Orbit "	<<std::endl;
	std::cout <<"------------------------"	<<std::endl;	
	float thrust_time_circ=1.0; 
	float thrust_step_circ=1.0;
	string filename3=(project_root / "output_files/orbit_3.csv").string();
	propagator circularization_orbit(onTransfer_orbit.last_sv,thrust_time_circ,thrust_step_circ);
	circularization_orbit.addPerturbation(&central_body, &cbody);
	circularization_orbit.addPerturbation(&thrust, &tparam);
	circularization_orbit.propagate();
	circularization_orbit.ephemeris_tofile(filename3);
	circularization_orbit.sv2oe(circularization_orbit.last_sv);
	std::cout <<"End-of-burn state vector : "<<circularization_orbit.last_sv <<std::endl;
	std::cout <<"Semimajor Axis : "<<circularization_orbit.a <<std::endl;

	std::cout <<"Eccentricity : "<<circularization_orbit.e <<std::endl;
	std::cout <<"True anomaly nu: "<<circularization_orbit.nu <<std::endl;

	//-----------------------------------------
	// Orbit 4 - Final Orbit
	//-----------------------------------------

	std::cout <<"---------------------"	<<std::endl;
	std::cout <<"- Final Orbit "		<<std::endl;
	std::cout <<"---------------------"	<<std::endl;
	float final_time=1.0; 
	float final_step=1.0;
	string filename4=(project_root / "output_files/orbit_4.csv").string();
	propagator final_orbit(circularization_orbit.last_sv, final_time, final_step);
	final_orbit.addPerturbation(&central_body, &cbody); // args: function and structure
	final_orbit.propagate();
	final_orbit.ephemeris_tofile(filename4);

	cout << endl;
	cout<<"End of processing!";
	return 0;

}
