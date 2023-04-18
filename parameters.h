#pragma once
#include <complex>
#include <cmath>
#include <vector>

typedef std::complex<double> dcomp;


struct Parameters
{
    double onsite_cor;         // onsite energy in the correlated region
	double onsite_l; //onsite in the left lead
	double onsite_r; //onsite in the right lead
    double hopping;      // the hopping 
	double gamma; //
	double bandwidth;    // the hopping 
    int num_ky_points; // this is the number of k in the y direction for the scattering region
    int num_kx_points; // This is the number of points in the x direction.
    double chemical_potential;
    double temperature;
    double e_upper_bound;       // this is the max energy value
    double e_lower_bound;       // this is the min energy value
    double hubbard_interaction; // this is the hubbard interaction
    double self_consistent_steps; // this is the number of self consistent steps my code needs
    double delta_gf; //delta in the scattering region
    double spin_up_occup;
    double spin_down_occup;
    double convergence;
    int interaction_order; 
    int steps; // number of energy points we take
    std::vector<double> energy;
    dcomp j1; // this is a complex number class defined within the complex library
    bool spin_polarised;
    int hubbard_steps;
    double delta_onsite;
    double delta_energy;
    int NIV_points;//number of IV points
	int NIV_start; //starting bias for the calculation. 0 for equilibrium
    double delta_v; //the voltage step between IV points
    std::vector<double> voltage_r;
    std::vector<double> voltage_l;
    static Parameters from_file();
};




double fermi_function(double energy, const Parameters &parameters);
void print_parameters(Parameters& parameters);
