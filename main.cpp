#include "parameters.h"
#include "slave_boson.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <limits>

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

void get_interacting_gf(const Parameters &parameters, const dcomp hamiltonian, const double &occupation, std::vector<dcomp> &interacting_gf){
    dcomp inverse_gf;
    for(int r = 0; r < parameters.steps; r++){
        inverse_gf = parameters.energy.at(r) + parameters.j1 * parameters.delta_gf - hamiltonian - occupation * parameters.hubbard_interaction;
        interacting_gf.at(r) = 1.0 / inverse_gf;
    }
}


void get_hamiltonian(Parameters const &parameters, const double kx, const double ky, dcomp &hamiltonian){
    hamiltonian = parameters.onsite_cor + 2 * parameters.hopping * (cos(kx) + cos(ky)); 
}


void get_local_gf(const Parameters &parameters, double occupation, std::vector<dcomp> &gf_local, const std::vector<std::vector<dcomp>> &hamiltonian){

    int n_x = parameters.num_kx_points;
    int n_y = parameters.num_ky_points;

    double num_k_points = n_x * n_y;
    
	for(int r = 0; r < parameters.steps; r++){
        gf_local.at(r) = 0;
    }


    for(int kx_i = 0; kx_i < n_x; kx_i++) {
        for(int ky_i = 0; ky_i < n_y; ky_i++) {
            std::vector<dcomp> gf_k_dep(parameters.steps, 0);
            get_interacting_gf(parameters, hamiltonian.at(kx_i).at(ky_i), occupation, gf_k_dep);
            for(int r = 0; r < parameters.steps; r++){
                gf_local.at(r) +=  gf_k_dep.at(r) * ( 1.0 / num_k_points);
            }
        }
    }
}

void get_spin_occupation(const Parameters &parameters, const std::vector<dcomp> &gf_r_up,
                        const std::vector<dcomp> &gf_r_down, double *spin_up, double *spin_down)
{
	double delta_energy = (parameters.e_upper_bound - parameters.e_lower_bound) / (double)(parameters.steps);
	double result_up = 0.0, result_down = 0.0;

	for (int r = 0; r < parameters.steps; r++) {
		result_up -= (delta_energy) * 2.0 * fermi_function(parameters.energy.at(r), parameters) * gf_r_up.at(r).imag();
		result_down -= (delta_energy) * 2.0 * fermi_function(parameters.energy.at(r), parameters) * gf_r_down.at(r).imag();
	}
	
	result_up = 1.0 / (2.0 * M_PI) * result_up;
	result_down = 1.0 / (2.0 * M_PI) * result_down;

    *spin_up = result_up;
    *spin_down = result_down;
}

double absolute_value(double num1) {
	return std::sqrt((num1 ) * (num1));
}

void get_difference(const Parameters &parameters, std::vector<dcomp> &gf_local_up, std::vector<dcomp> &old_green_function,
                double &difference, int &index){
	difference = - std::numeric_limits<double>::infinity();
	double old_difference = 0;
	double real_difference = 0, imag_difference = 0;
	for (int r = 0; r < parameters.steps; r++) {
		real_difference = absolute_value(gf_local_up.at(r).real() - old_green_function.at(r).real());
		imag_difference = absolute_value(gf_local_up.at(r).imag() - old_green_function.at(r).imag());
		//std::cout << gf_local_up.at(r)(i, j).real() << " " << old_green_function.at(r)(i, j).real() << std::endl;
		//std::cout << real_difference << "  " << imag_difference << "  "  << difference << "\n";
		difference = std::max(difference, std::max(real_difference, imag_difference));
		old_green_function.at(r) = gf_local_up.at(r);
		
        if (difference > old_difference) {
			index = r;
		}
		old_difference = difference;
	}
	//std::cout << "The difference is " << difference <<  std::endl;
}

void dmft(const Parameters &parameters, double *up_occup, double *down_occup, std::vector<dcomp> &gf_local_up, std::vector<dcomp> &gf_local_down,
    const std::vector<std::vector<dcomp>> &hamiltonian)
{
	double difference = std::numeric_limits<double>::infinity();
	int index, count = 0;

	std::vector<dcomp> old_green_function(parameters.steps, 0);

	while (difference > parameters.convergence && count < parameters.self_consistent_steps) {

        get_local_gf(parameters, *down_occup, gf_local_up, hamiltonian);
        get_local_gf(parameters, *up_occup, gf_local_down, hamiltonian);

		get_difference(parameters, gf_local_up, old_green_function, difference, index);

		if (difference < parameters.convergence) {
			break;
		}


        get_spin_occupation(parameters, gf_local_up, gf_local_down, up_occup, down_occup);
        //std::cout << "The count is " << count << ". The spin up occupation is " << *up_occup << ". The spin down occupation is " << *down_occup << std::endl;
		count++;
    }
	std::cout << "The count is " << count << ". The difference is " << difference << "\n" 
		"hubbard U = " << parameters.hubbard_interaction << " spin up occupation = " << *up_occup << " spin down occupation = " << *down_occup << std::endl;
}

int main(int argc, char **argv)
{
	Parameters parameters = Parameters::from_file();
	std::vector<double> kx(parameters.num_kx_points, 0);
	std::vector<double> ky(parameters.num_ky_points, 0);
    get_momentum_vectors(kx, ky, parameters);

	Interacting_GF boson(parameters);
	Interacting_GF fermion_up(parameters);
	Interacting_GF fermion_down(parameters);	

	double up_occupation = parameters.spin_up_occup, down_occupation = parameters.spin_down_occup;
	for (int m = parameters.NIV_start; m < parameters.NIV_points; m++) {
		impurity_solver(parameters, boson, fermion_up, fermion_down, m);
	}






    for (int kx_i = 0; kx_i < parameters.num_kx_points; kx_i++) {
		for (int ky_i = 0; ky_i < parameters.num_ky_points; ky_i++) {
            get_hamiltonian(parameters, kx.at(kx_i), ky.at(ky_i), hamiltonian.at(kx_i).at(ky_i));
        }
    }
	std::cout << "The initial spin occupation is " << parameters.spin_up_occup << ". The initial spin down occupation is " 
		<< parameters.spin_down_occup << std::endl;

	for (int i = 14; i < parameters.hubbard_steps; i++) {
		parameters.hubbard_interaction = (double) i * 0.5;
		std::cout << parameters.hubbard_interaction << std::endl;
    	dmft(parameters, &up_occupation, &down_occupation, gf_local_up, gf_local_down, hamiltonian);
		bool accepted = false;
		while (accepted != true) {
			double total_occupation = up_occupation + down_occupation;
			if (total_occupation < 1.4) {
				parameters.onsite_cor -= parameters.delta_onsite;
				accepted = false;
			} else if (total_occupation > 1.45) {
				parameters.onsite_cor += parameters.delta_onsite;
				accepted = false;				
			} else {
				accepted = true;
			}
			std::cout << accepted << std::endl;
			if (accepted == false) {
				up_occupation = parameters.spin_up_occup, down_occupation = parameters.spin_down_occup;
				gf_local_up.resize(parameters.steps, 0);
				gf_local_down.resize(parameters.steps, 0);		
				for (int kx_i = 0; kx_i < parameters.num_kx_points; kx_i++) {
					for (int ky_i = 0; ky_i < parameters.num_ky_points; ky_i++) {
            			get_hamiltonian(parameters, kx.at(kx_i), ky.at(ky_i), hamiltonian.at(kx_i).at(ky_i));
        			}
    			}
				dmft(parameters, &up_occupation, &down_occupation, gf_local_up, gf_local_down, hamiltonian);			
			}	
			total_occupation = up_occupation + down_occupation;
			std::cout << total_occupation << " " << up_occupation << " " <<  down_occupation << " " << parameters.onsite_cor << "\n";	
		}

		std::ostringstream ossgf;
		ossgf << i <<  ".gf.txt";
		std::string var = ossgf.str();
		std::ofstream gf_local_file;
		gf_local_file.open(var);
		for (int r = 0; r < parameters.steps; r++) {
			gf_local_file << parameters.energy.at(r) << "  " << gf_local_up.at(r).real() << "   " << gf_local_up.at(r).imag() << " "
				<< gf_local_down.at(r).real() << "   " << gf_local_down.at(r).imag() << "  \n";
		}
		gf_local_file.close();
	}


	

    //std::ofstream fermi_function_right;
	//fermi_function_right.open("fermi_function_right.txt");
	//for (int r = 0; r < parameters.steps; r++) {
    //    fermi_function_right << parameters.energy.at(r) << "  " << fermi_function(parameters.energy.at(r) + 0.5, parameters) << "   " 
    //        << fermi_function(parameters.energy.at(r) - 0.5, parameters) << " \n";
    //}

            //std::cout << parameters.energy.at(r) << "  " << gf_local_up.at(r).real() << "   " << gf_local_up.at(r).imag() << " "
			//<< gf_local_down.at(r).real() << "   " << gf_local_down.at(r).imag() << "  \n";
	//fermi_function_right.close();	

}