#include "parameters.h"
#include "impurity_solver.h"
#include "green_function.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <limits>
#include <iomanip>  

void get_momentum_vectors(std::vector<double> &kx, std::vector<double> &ky, Parameters &parameters) {
	for (int i = 0; i < parameters.num_kx_points; i++) {
		if (parameters.num_kx_points != 1) {
			kx.at(i) = 2 * M_PI * i / parameters.num_kx_points;
		} else if (parameters.num_kx_points == 1) {
			kx.at(i) = M_PI / 2.0;
		}
	}

	for (int i = 0; i < parameters.num_ky_points; i++) {
		if (parameters.num_ky_points != 1) {
			ky.at(i) = 2 * M_PI * i / parameters.num_ky_points;
		} else if (parameters.num_ky_points == 1) {
			ky.at(i) = M_PI / 2.0;
		}
	}
}

void print_gf_lesser_greater(const Parameters &parameters, const int &voltage_step, std::vector<dcomp> &gf_lesser, std::vector<dcomp> &gf_greater) {
	
	std::ostringstream ossgf;
	ossgf << voltage_step << ".greater_gf.dat";
	std::string var = ossgf.str();
	std::ofstream gf_greater_file;
	gf_greater_file.open(var);
	for (int r = 0; r < parameters.steps; r++) {
		gf_greater_file << parameters.energy.at(r) << "  " << gf_greater.at(r).real() << "   " << gf_greater.at(r).imag() << "  \n";
	}
	gf_greater_file.close();

    ossgf.str("");
    ossgf.clear();
	ossgf << voltage_step << ".lesser_gf.dat";
	var = ossgf.str();
	std::ofstream gf_lesser_file;
	gf_lesser_file.open(var);
	for (int r = 0; r < parameters.steps; r++) {
		gf_lesser_file << parameters.energy.at(r) << "  " << gf_lesser.at(r).real() << "   " << gf_lesser.at(r).imag() << "  \n";
	}
	gf_lesser_file.close();

}

void print_dos(const Parameters &parameters, const int &voltage_step, std::vector<dcomp> &gf_lesser, std::vector<dcomp> &gf_greater) {
	std::ostringstream ossgf;
	ossgf << voltage_step << ".dos.dat";
	double dos_integral = 0.0;
	std::string var = ossgf.str();
	std::ofstream dos;
	dos.open(var);
	for (int r = 0; r < parameters.steps; r++) {
		dos << parameters.energy.at(r) << "  " << 0.5 * (gf_lesser.at(r).imag() - gf_greater.at(r).imag()) << "  \n";
		dos_integral += 0.5 * (gf_lesser.at(r).imag() - gf_greater.at(r).imag());
	}
	std::cout << "The dos integrated for all energies is " << dos_integral * parameters.delta_energy << std::endl;
	dos.close();
}

void get_lesser_greater_gf(const Parameters &parameters, const Interacting_GF &boson, const Interacting_GF &fermion, 
     const double &z_prefactor, std::vector<dcomp> &gf_lesser, std::vector<dcomp> &gf_greater) {

	//std::cout << parameters.steps << std::endl;

	for (int r = 0; r < parameters.steps; r++) {
		for (int i = 0; i < parameters.steps; i++) {
			if (((i + r) >= (parameters.steps / 2)) && ((i + r) < 3 * (parameters.steps / 2))) {
				//std::cout << r << " " << i << "  " << parameters.energy.at(r) + parameters.energy.at(i)  << "  " << parameters.energy.at(r + i - (parameters.steps / 2)) << "\n ";
				gf_lesser.at(r) += boson.greater_gf.at(i) * fermion.lesser_gf.at(i + r - (parameters.steps / 2));
				gf_greater.at(r) += boson.lesser_gf.at(i) * fermion.greater_gf.at(i + r - (parameters.steps / 2));				
			}		
		}
		//std::cout << parameters.steps / 2 << std::endl;
		gf_lesser.at(r) = gf_lesser.at(r) * z_prefactor * parameters.j1 / (2.0 * M_PI) * parameters.delta_energy;
		gf_greater.at(r) = gf_greater.at(r) * z_prefactor * parameters.j1 / (2.0 * M_PI) * parameters.delta_energy;		
	}
}

void get_current(const Parameters &parameters, const std::vector<dcomp> &gf_lesser, const std::vector<dcomp> &gf_greater, const int voltage_step,
	 double &current_left, double &current_right) {
	for (int r = 0; r < parameters.steps; r++) {
		double hybridisation = (parameters.energy.at(r) * parameters.energy.at(r) + parameters.bandwidth * parameters.bandwidth);
		current_left += ((1.0 - fermi_function(parameters.energy.at(r) - parameters.voltage_l[voltage_step], parameters)) * gf_lesser.at(r).imag() 
			+ fermi_function(parameters.energy.at(r) - parameters.voltage_l[voltage_step], parameters) * gf_greater.at(r).imag()) / hybridisation;
		current_right += ((1.0 - fermi_function(parameters.energy.at(r) - parameters.voltage_r[voltage_step], parameters)) * gf_lesser.at(r).imag() 
			+ fermi_function(parameters.energy.at(r) - parameters.voltage_r[voltage_step], parameters) * gf_greater.at(r).imag()) / hybridisation;
	}
	current_left = current_left * parameters.delta_energy * parameters.gamma * parameters.bandwidth * parameters.bandwidth;
	current_right = current_right * parameters.delta_energy * parameters.gamma * parameters.bandwidth * parameters.bandwidth;	
}

void get_conductance(const Parameters &parameters, const int &voltage_step, std::vector<dcomp> &gf_lesser, std::vector<dcomp> &gf_greater) {
	double conductance = 0;

	for (int r = 0; r < parameters.steps; r++) {
		double dos = 0.5 * (gf_lesser.at(r).imag() - gf_greater.at(r).imag());
		double coupling = parameters.gamma * parameters.bandwidth * parameters.bandwidth / (parameters.energy.at(r) * parameters.energy.at(r) + parameters.bandwidth * parameters.bandwidth);
		double derivative =  exp((parameters.energy.at(r) - parameters.chemical_potential) / parameters.temperature) / (parameters.temperature * 
		(1.0 + exp((parameters.energy.at(r) - parameters.chemical_potential) / parameters.temperature)) *
		(1.0 + exp((parameters.energy.at(r) - parameters.chemical_potential) / parameters.temperature))); 
		if (std::isnan(derivative) == true) { //this is a big hack but i think it is just numerical problems that the derivative can be nan away from the chemical potential
			derivative = 0;
		}

		//std::cout << parameters.energy.at(r) << " " << dos << " " << coupling << " " << derivative << "\n";
		conductance += dos * coupling * derivative;
	}

	std::cout << "The conductance is " << conductance * parameters.delta_energy * 2.0 * M_PI << "\n";
}


void get_spin_occupation(const Parameters &parameters, const std::vector<dcomp> &gf_lesser, double &z_prefactor)
{
	double occupation = 0.0;

	for (int r = 0; r < parameters.steps; r++) {
		occupation += gf_lesser.at(r).imag();
	}
	
	occupation = occupation * parameters.delta_energy * z_prefactor / (2.0 * M_PI);
	std::cout << "the pseudo fermion occupation is " << occupation << std::endl;
}

int main(int argc, char **argv)
{
	Parameters parameters = Parameters::from_file();
	//std::vector<double> kx(parameters.num_kx_points, 0);
	//std::vector<double> ky(parameters.num_ky_points, 0);
    //get_momentum_vectors(kx, ky, parameters);
	std::vector<double> current_left(parameters.NIV_points, 0);
	std::vector<double> current_right(parameters.NIV_points, 0);

	for (int m = parameters.NIV_start; m < parameters.NIV_points; m++) {
		std::cout << std::setprecision(15) << "The voltage is " << parameters.voltage_l[m] * 2 << ". \n";
		Interacting_GF boson(parameters);
		Interacting_GF fermion_up(parameters);
		Interacting_GF fermion_down(parameters);	

		std::vector<dcomp> gf_lesser_up(parameters.steps, 0), gf_greater_up(parameters.steps, 0);
		std::vector<dcomp> gf_lesser_down(parameters.steps, 0), gf_greater_down(parameters.steps, 0);
		double z_prefactor = 0;

		double up_occupation = parameters.spin_up_occup, down_occupation = parameters.spin_down_occup;

		impurity_solver(parameters, boson, fermion_up, fermion_down, m, z_prefactor);


		std::cout << "The ratio of Z_0 / Z_1 is " << z_prefactor << std::endl;
		get_lesser_greater_gf(parameters, boson, fermion_up, z_prefactor, gf_lesser_up, gf_greater_up);
		get_lesser_greater_gf(parameters, boson, fermion_down, z_prefactor, gf_lesser_down, gf_greater_down);

		//boson.print_green_function(parameters, m, "boson");
		//fermion_up.print_green_function(parameters, m, "fermion_up");
		//fermion_down.print_green_function(parameters, m, "fermion_down");
	
		print_gf_lesser_greater(parameters, m, gf_lesser_up, gf_greater_up);
		print_dos(parameters, m, gf_lesser_up, gf_greater_up);
		get_spin_occupation(parameters, fermion_up.lesser_gf, z_prefactor);
		if (m == 0) {
			get_conductance(parameters, m, gf_lesser_up, gf_greater_up);
		}
		get_current(parameters, gf_lesser_up, gf_greater_up, m, current_left.at(m), current_right.at(m));
	}

	for (int m = 0; m < parameters.NIV_points; m++) {
		std::cout << "The voltage is " << parameters.voltage_l[m] * 2 << ". The left current is " << current_left.at(m) << ". The right current is " << current_right.at(m) << "\n";
	}
}