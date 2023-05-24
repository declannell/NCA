#include "parameters.h"
#include <complex>  //this contains complex numbers and trig functions
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void simple_tokenizer(std::string s, std::string &variable, std::string &value)
{
    std::stringstream ss(s);
    std::string word;
	int count = 0;
    while (ss >> word) {
        //std::cout << word << std::endl;
		if (count == 0) {
			variable = word;
		} else if (count == 1) {
			value = word;
		}
		count++;
    }
}


Parameters Parameters::from_file()
{
	Parameters parameters;
	std::string line, variable, value;
	std::ifstream input_file;

	input_file.open("input_file");

	if(!input_file.is_open())
	{
		std::cout << "The input_file doesn't exist \n" << std::endl;
		std::exit(1);
	} 
    double onsite_cor;         // onsite energy in the correlated region
	double onsite_l; //onsite in the left lead
	double onsite_r; //onsite in the right lead
    double hopping;      // the hopping 
	double gamma; //
	double bandwidth; //this is the bandwidth of the leads
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



	while (getline(input_file, line)) {
			//std::cout << line << '\n';
			simple_tokenizer(line, variable, value);
			//std::cout << "The variable name is " << variable << " The value is  " << value << std::endl;
			if (variable == "onsite_cor") {
				parameters.onsite_cor = std::stod(value);
			} else if (variable == "onsite_l") {
				parameters.onsite_l = std::stod(value);
			} else if (variable == "onsite_r") {
				parameters.onsite_r = std::stod(value);
			} else if (variable == "hopping") {
				parameters.hopping = std::stod(value);
			} else if (variable == "bandwidth") {
				parameters.bandwidth = std::stod(value);
			} else if (variable == "gamma") {
				parameters.gamma = std::stod(value);
			} else if (variable == "num_ky_points") {
				parameters.num_ky_points = std::stoi(value);
			} else if (variable == "num_kx_points") {
				parameters.num_kx_points = std::stoi(value);
			} else if (variable == "chemical_potential") {
				parameters.chemical_potential = std::stod(value);
			} else if (variable == "temperature") {
				parameters.temperature = std::stod(value);
			} else if (variable == "e_upper_bound") {
				parameters.e_upper_bound = std::stod(value);
			} else if (variable == "e_lower_bound") {
				parameters.e_lower_bound = std::stod(value);
			} else if (variable == "hubbard_interaction") {
				parameters.hubbard_interaction = std::stod(value);
			} else if (variable == "self_consistent_steps") {
				parameters.self_consistent_steps = std::stod(value);
			} else if (variable == "delta_gf") {
				parameters.delta_gf = std::stod(value);
			} else if (variable == "spin_up_occup") {
				parameters.spin_up_occup = std::stod(value);
			} else if (variable == "spin_down_occup") {
				parameters.spin_down_occup = std::stod(value);
			} else if (variable == "convergence") {
				parameters.convergence = std::stod(value);
			} else if (variable == "hubbard_steps") {
				parameters.hubbard_steps = std::stod(value);
			} else if (variable == "delta_onsite") {
				parameters.delta_onsite = std::stod(value);
			} else if (variable == "steps") {
				parameters.steps = std::stoi(value);
			} else if (variable == "interaction_order") {
				parameters.interaction_order = std::stoi(value);
			} else if (variable == "spin_polarised") {
				std::istringstream(value) >> parameters.spin_polarised;
			} else if (variable == "NIV_points") {
				parameters.NIV_points = std::stoi(value);
			} else if (variable == "NIV_start") {
                parameters.NIV_start = std::stoi(value);
            } else if (variable == "delta_v") {
				parameters.delta_v = std::stod(value);
			} 
	}
	input_file.close();
	
	if (parameters.hubbard_interaction == 0.0) {
		parameters.interaction_order =
		    0;  // this is the order the green function will be calculated too in terms of interaction strength. this can be equal to 0 , 1 or 2//
	} 

	parameters.voltage_l.resize(parameters.NIV_points);
	parameters.voltage_r.resize(parameters.NIV_points);
	for (int i = 0; i < parameters.NIV_points; i++) {
		parameters.voltage_l.at(i) = parameters.delta_v * (double)(i);
		parameters.voltage_r.at(i) = -parameters.delta_v * (double)(i);
	}

	parameters.energy.resize(parameters.steps);

	parameters.j1 = -1;
	parameters.j1 = sqrt(parameters.j1);

	parameters.delta_energy =
	    (parameters.e_upper_bound - parameters.e_lower_bound) / (double)parameters.steps;

	for (int i = 0; i < parameters.steps; i++) {
		parameters.energy.at(i) = parameters.e_lower_bound + parameters.delta_energy * (double)i;
	}
	return parameters;


}

double fermi_function(double energy, const Parameters& parameters)
{
	if (parameters.temperature == 0) {
		if (energy < parameters.chemical_potential) {
			return 1.0;
		} else {
			return 0.0;
		}
	} else {
	}
	return 1.0 / (1.0 + exp((energy - parameters.chemical_potential) / parameters.temperature));
}
//Parameters params = Parameters::from_file();
void print_parameters(Parameters& parameters)
{
	std::cout << " .onsite_cor = " << parameters.onsite_cor << std::endl;
	std::cout << "hopping = " << parameters.hopping << std::endl;
	std::cout << "num_ky_points = " << parameters.num_ky_points << std::endl;
	std::cout << "num_kx_points = " << parameters.num_kx_points << std::endl;
	std::cout << "chemical_potential = " << parameters.chemical_potential << std::endl;
	std::cout << "temperature = " << parameters.temperature << std::endl;
	std::cout << "e_upper_bound = " << parameters.e_upper_bound << std::endl;
	std::cout << "e_lower_bound = " << parameters.e_lower_bound << std::endl;
	std::cout << "hubbard_interaction = " << parameters.hubbard_interaction << std::endl;
	std::cout << "self_consistent_steps = " << parameters.self_consistent_steps << std::endl;
	std::cout << "delta_gf = " << parameters.delta_gf << std::endl;
    std::cout << "parameters.interaction_order = " << parameters.interaction_order << std::endl;
	std::cout << "parameters.steps = " << parameters.steps << std::endl;
	std::cout << "parameters.spin_up_occup = " << parameters.spin_up_occup << std::endl;
	std::cout << "parameters.spin_down_occup = " << parameters.spin_down_occup << std::endl;
	std::cout << "parameters.spin_polarised = " << parameters.spin_polarised << std::endl;
}
