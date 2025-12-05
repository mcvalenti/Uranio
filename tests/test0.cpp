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
    float propaga_time=90*60;  // [s] 90 minutes
	float propaga_step=0.5;
    string filename1=(project_root / "output_files/orbit0.csv").string();
    long double val[7]= {6.016759234991846e+06, -3.330707211868651e+06, 0.0, 3.651024513484383e+03,  6.766879765678871e+03, 0.0, 2000.0};
    LDVector validation_vec(val, 7); 
    propagator circular_orbit;
    LDVector ld(circular_orbit.current_sv);
    circular_orbit.reset_init(ld,propaga_time,propaga_step);
    // Propagation
    std::cout <<"Init state vector : "<< ld <<std::endl;
    circular_orbit.addPerturbation(&central_body, &cbody); 
	circular_orbit.propagate();
    circular_orbit.ephemeris_tofile(filename1);
    // Validation
    std::cout <<"End-of-burn state vector : "<< circular_orbit.last_sv<<std::endl;
    double tol=100;
    int verify_counter = 0;
    for (int i=0; i<=6; i++){
        if (abs(circular_orbit.last_sv[i]-validation_vec[i]) < tol){
            verify_counter = verify_counter+1;
        }
    }
    if (verify_counter == 7){
        std::cout <<"Test0 Successfull !!!: "<< circular_orbit.last_sv<<std::endl;
    }
    
    return 0;
}