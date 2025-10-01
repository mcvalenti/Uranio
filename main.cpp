#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// ----------------------------
// Struct to define the state
// ----------------------------
struct State {
    double x;      // posición
    double v;      // velocidad
    double lx;     // coestado asociado a x
    double lv;     // coestado asociado a v
};

// ----------------------------
// Dynamics of the Extended Sys
// ----------------------------
State dynamics(const State &s) {
    State ds;
    double u = -0.5 * s.lv;  // optimal control

    ds.x  = s.v;        // dx/dt
    ds.v  = u;          // dv/dt
    ds.lx = 0.0;        // dλx/dt
    ds.lv = -s.lx;      // dλv/dt

    return ds;
}

// ----------------------------
// Paso de integración RK4
// ----------------------------
State RK4_step(const State &s, double h) {
    State k1 = dynamics(s);

    State s2 = { s.x + 0.5*h*k1.x, s.v + 0.5*h*k1.v,
                 s.lx + 0.5*h*k1.lx, s.lv + 0.5*h*k1.lv };
    State k2 = dynamics(s2);

    State s3 = { s.x + 0.5*h*k2.x, s.v + 0.5*h*k2.v,
                 s.lx + 0.5*h*k2.lx, s.lv + 0.5*h*k2.lv };
    State k3 = dynamics(s3);

    State s4 = { s.x + h*k3.x, s.v + h*k3.v,
                 s.lx + h*k3.lx, s.lv + h*k3.lv };
    State k4 = dynamics(s4);

    State out;
    out.x  = s.x  + (h/6.0)*(k1.x  + 2*k2.x  + 2*k3.x  + k4.x);
    out.v  = s.v  + (h/6.0)*(k1.v  + 2*k2.v  + 2*k3.v  + k4.v);
    out.lx = s.lx + (h/6.0)*(k1.lx + 2*k2.lx + 2*k3.lx + k4.lx);
    out.lv = s.lv + (h/6.0)*(k1.lv + 2*k2.lv + 2*k3.lv + k4.lv);

    return out;
}

// ----------------------------
// Propagation until T
// ----------------------------
State integrate_system(const State &s0, double T, int Nsteps) {
    double h = T / Nsteps;
    State s = s0;
    for (int i=0; i<Nsteps; i++) {
        s = RK4_step(s, h);
    }
    return s;
}

// ----------------------------
// Shooting (Newton 2D)
// ----------------------------
int main() {
    
    double x0 = 0.0, v0 = 0.0; 	// Init conditions   
    double xf = 10.0, vf = 0.0; // End Restrictions or desired conditions
    double T = 5.0;

    // Suggested initial costates
    State s0 = {x0, v0, 1.0, 1.0};

    double tol = 1e-6;
    int max_iter = 20;

    for (int iter=0; iter<max_iter; iter++) {
        // Propagation with lambda guess
        State sf = integrate_system(s0, T, 500);

        // Error in final conditions
        double F1 = sf.x - xf;
        double F2 = sf.v - vf;

        cout << "Iter " << iter << "  Error = (" << F1 << ", " << F2 << ")\n";

        if (fabs(F1) < tol && fabs(F2) < tol) {
            cout << "Converged!\n";
            break;
        }

        // Newton implementation in 2D - Jacobian w.r.t  lx0, lv0
        double eps = 1e-5;
        // Vary lx0
        State s_test = {x0, v0, s0.lx + eps, s0.lv};
        State sf_test = integrate_system(s_test, T, 500);
        double dF1_dlx = (sf_test.x - xf - F1)/eps;
        double dF2_dlx = (sf_test.v - vf - F2)/eps;
        // Vary lv0
        s_test = {x0, v0, s0.lx, s0.lv + eps};
        sf_test = integrate_system(s_test, T, 500);
        double dF1_dlv = (sf_test.x - xf - F1)/eps;
        double dF2_dlv = (sf_test.v - vf - F2)/eps;

        // solve the linear system J*delta = -F
        double det = dF1_dlx*dF2_dlv - dF1_dlv*dF2_dlx;
        double d_lx = (-F1*dF2_dlv + F2*dF1_dlv) / det;
        double d_lv = (-dF1_dlx*F2 + F1*dF2_dlx) / det;

        // Update guess
        s0.lx += d_lx;
        s0.lv += d_lv;
    }

    cout << "Lx0 = " << s0.lx << " , Lv0 = " << s0.lv << endl;
	system("pause"); // Wait until press a key.
    return 0;
}
