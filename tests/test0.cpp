#include <iostream>
#include <filesystem>
#include <string>
#include "propagators.h"
#include "constants.h"

int main() {
    // ---- OUTPUT FILES PATH -------------------    
    auto exec_path = std::filesystem::current_path();
    std::cout << "Running from: " << exec_path << std::endl;

    // Up 1 level if in build
    auto project_root = exec_path;
    if (exec_path.filename() == "build")
        project_root = exec_path.parent_path();
    // ---- ------------------- ------------------ 

    cBody_param cbody;
	cbody.mu=GM;
    float thrust_time=260.0;  // [s]
	float thrust_step=0.5;
	thrust_param tparam;
	tparam.isp = 300.0;
	tparam.thrust = 10000.0; // [N]
    string filename1=(project_root / "output_files/orbit0.csv").string();
    propagator leop_orbit;
    std::cout <<"Init state vector : "<< leop_orbit.current_sv<<std::endl; 
    LDVector ld(leop_orbit.current_sv);
    leop_orbit.reset_init(ld,thrust_time,thrust_step);
    // Propagation
    std::cout <<"Init state vector : "<< ld <<std::endl;
    leop_orbit.addPerturbation(&central_body, &cbody); 
    leop_orbit.addPerturbation(&thrust, &tparam);
	leop_orbit.propagate();
    leop_orbit.ephemeris_tofile(filename1);
    
    std::cout <<"End-of-burn state vector : "<< leop_orbit.last_sv<<std::endl;
    
    return 0;
}