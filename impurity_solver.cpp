#include "parameters.h"
#include "green_function.h"
#include "impurity_solver.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <limits>

void get_greater_lesser_se_fermion(const Parameters &parameters, const Interacting_GF &boson, Interacting_GF &fermion, int voltage_step) {
    
    fermion.lesser_se.clear();
    fermion.greater_se.clear();

    fermion.lesser_se.resize(parameters.steps, 0);
    fermion.greater_se.resize(parameters.steps, 0);

    
    //this uses the WBA 
    for (int r = 0; r < parameters.steps; r++) {
        //std::cout <<  boson.lesser_gf.at(r) << std::endl;
        for (int i = 0; i < parameters.steps; i++) {

            double coupling = parameters.gamma * parameters.bandwidth * parameters.bandwidth /
             ((parameters.energy.at(r) - parameters.energy.at(i)) * (parameters.energy.at(r) - parameters.energy.at(i)) + parameters.bandwidth * parameters.bandwidth);
            //if (r == 0) {
            //    std::cout << parameters.energy.at(r)  << " " << parameters.energy.at(i) << " " << fermi_function(parameters.energy.at(r) - parameters.energy.at(i) - parameters.voltage_l[voltage_step], parameters) << std::endl;
            //}
            fermion.greater_se.at(r) += 1.0 / (2 * M_PI) * coupling * boson.greater_gf.at(i) * (2.0 - 
                fermi_function(parameters.energy.at(r) - parameters.energy.at(i) - parameters.voltage_l[voltage_step], parameters) - 
                fermi_function(parameters.energy.at(r) - parameters.energy.at(i) - parameters.voltage_r[voltage_step], parameters)) * parameters.delta_energy;

            fermion.lesser_se.at(r) -= 1.0 / (2 * M_PI) * coupling * boson.lesser_gf.at(i) * ( 
                fermi_function(parameters.energy.at(r) - parameters.energy.at(i) - parameters.voltage_l[voltage_step], parameters) + 
                fermi_function(parameters.energy.at(r) - parameters.energy.at(i) - parameters.voltage_r[voltage_step], parameters)) * parameters.delta_energy;

            //fermion.lesser_se.at(r) += - 1.0 / (2 * M_PI) * parameters.gamma * boson.lesser_gf.at(i) *  
            //    (fermi_function(parameters.energy.at(r) - parameters.energy.at(i) - parameters.voltage_l[voltage_step], parameters) +
            //    fermi_function(parameters.energy.at(r) - parameters.energy.at(i) - parameters.voltage_r[voltage_step], parameters)) * parameters.delta_energy;
        }

    //std::cout << parameters.energy.at(r) << " " <<  fermion.lesser_se.at(r) << " " << fermion.greater_se.at(r) << "\n";
    }

}

void get_retarded_se(const Parameters &parameters, Interacting_GF &green_function) {
    
    green_function.retarded_se.clear();
    green_function.retarded_se.resize(parameters.steps, 0);  

    for (int r = 0; r < parameters.steps; r++) {
        green_function.retarded_se.at(r) = 0.5 * green_function.greater_se.at(r);
    }

    for (int r = 0; r < parameters.steps; r++) {
        double se_real = 0;
        for (int i = 0; i < parameters.steps; i++) {
            if (i != r) {
                se_real += green_function.retarded_se.at(i).imag() / (parameters.energy.at(i) - parameters.energy.at(r));
            }
        }
        green_function.retarded_se.at(r) += se_real * parameters.delta_energy / M_PI;
    }

    //for (int r = 0; r < parameters.steps; r++) {
    //    for (int i = 0; i < parameters.steps; i++) {
    //        dcomp term = parameters.energy.at(r) - parameters.energy.at(i) + parameters.delta_gf * parameters.j1;
    //        green_function.retarded_se.at(r) += green_function.greater_se.at(i) / term;
    //    }
    //    green_function.retarded_se.at(r) = green_function.retarded_se.at(r) * parameters.j1 * parameters.delta_energy / (2 * M_PI);
    //}
}

void get_retarded_gf_fermion(const Parameters &parameters, Interacting_GF &green_function) {   
    for( int r = 0; r < parameters.steps; r++) {
        green_function.retarded_gf.at(r) = 1.0 / (parameters.energy.at(r) - parameters.onsite_cor - green_function.retarded_se.at(r));
    }
}

void get_greater_lesser_gf(const Parameters &parameters, Interacting_GF &green_function) {
    for (int r = 0; r < parameters.steps; r++) {
        green_function.greater_gf.at(r) = green_function.retarded_gf.at(r) * green_function.greater_se.at(r) * std::conj(green_function.retarded_gf.at(r));
        green_function.lesser_gf.at(r) = green_function.retarded_gf.at(r) * green_function.lesser_se.at(r) * std::conj(green_function.retarded_gf.at(r));
    }
}

void get_greater_lesser_gf(const Parameters &parameters, Interacting_GF &green_function, int count) {
    if (count < 10) {
        for (int r = 0; r < parameters.steps; r++) {
            green_function.greater_gf.at(r) = 0.5 * green_function.greater_gf.at(r) + 0.5 * green_function.retarded_gf.at(r) * green_function.greater_se.at(r) * std::conj(green_function.retarded_gf.at(r));
            green_function.lesser_gf.at(r) = 0.5 * green_function.lesser_gf.at(r) + 0.5 * green_function.retarded_gf.at(r) * green_function.lesser_se.at(r) * std::conj(green_function.retarded_gf.at(r));
        }
    } else if (count > 10 && count < 0) {
        for (int r = 0; r < parameters.steps; r++) {
            green_function.greater_gf.at(r) = 0.75 * green_function.greater_gf.at(r) + 0.25 * green_function.retarded_gf.at(r) * green_function.greater_se.at(r) * std::conj(green_function.retarded_gf.at(r));
            green_function.lesser_gf.at(r) = 0.75 * green_function.lesser_gf.at(r) + 0.25 * green_function.retarded_gf.at(r) * green_function.lesser_se.at(r) * std::conj(green_function.retarded_gf.at(r));
        }    
    } else {
        for (int r = 0; r < parameters.steps; r++) {
            green_function.greater_gf.at(r) = 0.9 * green_function.greater_gf.at(r) + 0.1 * green_function.retarded_gf.at(r) * green_function.greater_se.at(r) * std::conj(green_function.retarded_gf.at(r));
            green_function.lesser_gf.at(r) = 0.9 * green_function.lesser_gf.at(r) + 0.1 * green_function.retarded_gf.at(r) * green_function.lesser_se.at(r) * std::conj(green_function.retarded_gf.at(r));
        }
    }
}

void get_greater_lesser_se_boson(const Parameters &parameters, Interacting_GF &boson, const Interacting_GF &fermion_up,
    const Interacting_GF &fermion_down, int voltage_step) {

    boson.lesser_se.clear();
    boson.greater_se.clear();

    boson.lesser_se.resize(parameters.steps, 0);
    boson.greater_se.resize(parameters.steps, 0);

    //std::cout << "make sure the spin summation is correct in the spin non-polarised case. \n";

    for (int r = 0; r < parameters.steps; r++) {
        //std::cout << boson.lesser_se.at(r) << std::endl;
        for (int i = 0; i < parameters.steps; i++) {

            double coupling = parameters.gamma * parameters.bandwidth * parameters.bandwidth /
             ((parameters.energy.at(i) - parameters.energy.at(r))* (parameters.energy.at(i) - parameters.energy.at(r)) + parameters.bandwidth * parameters.bandwidth);

            boson.greater_se.at(r) += (1.0 / (2.0 * M_PI)) * coupling * parameters.delta_energy * (fermion_up.greater_gf.at(i) + fermion_down.greater_gf.at(i))
                * (fermi_function(parameters.energy.at(i) - parameters.energy.at(r) - parameters.voltage_l[voltage_step], parameters) + 
                   fermi_function(parameters.energy.at(i) - parameters.energy.at(r) - parameters.voltage_r[voltage_step], parameters));

            boson.lesser_se.at(r) -= (1.0 / (2.0 * M_PI)) * parameters.delta_energy * coupling * (fermion_up.lesser_gf.at(i) + fermion_down.lesser_gf.at(i)) *
                (2.0 - fermi_function(parameters.energy.at(i) - parameters.energy.at(r) - parameters.voltage_l[voltage_step], parameters) 
                    -  fermi_function(parameters.energy.at(i) - parameters.energy.at(r) - parameters.voltage_r[voltage_step], parameters));
            //if (r == 1) {
            //    std::cout << fermion_up.lesser_gf.at(i) << " " << 2.0 - fermi_function(parameters.energy.at(r) - parameters.energy.at(i) - parameters.voltage_l[voltage_step], parameters) 
            //        - fermi_function(parameters.energy.at(r) - parameters.energy.at(i) - parameters.voltage_r[voltage_step], parameters) << std::endl;
            //}
        }
    }
}

void get_retarded_gf_boson(const Parameters &parameters, Interacting_GF &green_function) {   
    for( int r = 0; r < parameters.steps; r++) {
        green_function.retarded_gf.at(r) = 1.0 / (parameters.energy.at(r) - green_function.retarded_se.at(r));
    }
}

double absolute_value(double num1) {
	return std::sqrt((num1 ) * (num1));
}

void get_difference(const Parameters &parameters, std::vector<dcomp> &old_self_energy, std::vector<dcomp> &new_self_energy,
                double &difference){
    difference = 0;
	double old_difference = 0;
	double real_difference = 0, imag_difference = 0;
	for (int r = 0; r < parameters.steps; r++) {
		real_difference = absolute_value(old_self_energy.at(r).real() - new_self_energy.at(r).real());
		imag_difference = absolute_value(old_self_energy.at(r).imag() - new_self_energy.at(r).imag());
		difference = std::max(difference, std::max(real_difference, imag_difference));
		old_self_energy.at(r) = new_self_energy.at(r);
    }
}
void intialise_gf(const Parameters &parameters, Interacting_GF &boson, Interacting_GF &fermion_up, Interacting_GF &fermion_down) {
    int a = parameters.steps / 3;
    for (int r = a; r < a *2; r++) {
        //fermion_up.lesser_gf.at(r) = parameters.j1;
        //fermion_down.lesser_gf.at(r) = parameters.j1;
        boson.lesser_gf.at(r) = - 1.0 * parameters.j1;
        //fermion_up.greater_gf.at(r) = - parameters.j1;
        //fermion_down.greater_gf.at(r) = - parameters.j1;
        boson.greater_gf.at(r) =  - 1.0 * parameters.j1;
    }
}

double get_z_prefactor(const Parameters &parameters, Interacting_GF &boson, Interacting_GF &fermion_up, Interacting_GF &fermion_down) {
    double z_prefactor = 0;

    for (int r = 0; r < parameters.steps; r++) {
        z_prefactor += - boson.lesser_gf.at(r).imag() + fermion_up.lesser_gf.at(r).imag() + fermion_down.lesser_gf.at(r).imag();
    }

    return  parameters.delta_energy * z_prefactor / (2.0 * M_PI);
}

void test_retarded_gf(const Parameters &parameters, Interacting_GF &boson, Interacting_GF &fermion_up, Interacting_GF &fermion_down) {

    double test_boson = 0.0, test_fermion_up = 0.0, test_fermion_down = 0.0;

    for (int r = 0; r < parameters.steps; r++) {
        test_boson += boson.retarded_gf.at(r).imag();
        test_fermion_up += fermion_up.retarded_gf.at(r).imag();
        test_fermion_down += fermion_down.retarded_gf.at(r).imag();
    }

    std::cout << "The imaginary part of the boson retarded gf integrates to " <<  test_boson * -1.0 * parameters.delta_energy / (M_PI) << "\n";
    std::cout << "The imaginary part of the fermion up retarded gf integrates to " <<  test_fermion_up * -1.0 * parameters.delta_energy / (M_PI) << "\n";
    std::cout << "The imaginary part of the fermion down retarded gf integrates to " <<  test_fermion_down * -1.0 * parameters.delta_energy / (M_PI) << "\n";
}

void impurity_solver(const Parameters &parameters, Interacting_GF &boson, Interacting_GF &fermion_up, Interacting_GF &fermion_down, int voltage_step,
     double &z_prefactor) {

    std::vector<dcomp> old_self_energy(parameters.steps, 0);
    int count = 0;
    double difference = std::numeric_limits<double>::infinity();
    std::vector<double> coupling;

    intialise_gf(parameters, boson, fermion_up, fermion_down);

    while (difference > parameters.convergence && count < parameters.self_consistent_steps) {
        if (parameters.spin_polarised == 1) {
            
            get_greater_lesser_se_fermion(parameters, boson, fermion_up, voltage_step);
            get_greater_lesser_se_fermion(parameters, boson, fermion_down, voltage_step);

            get_retarded_se(parameters, fermion_up);
            get_retarded_se(parameters, fermion_down);

            get_retarded_gf_fermion(parameters, fermion_up);
            get_retarded_gf_fermion(parameters, fermion_down);

            get_greater_lesser_gf(parameters, fermion_up);
            get_greater_lesser_gf(parameters, fermion_down);

            get_greater_lesser_se_boson(parameters, boson, fermion_up, fermion_down, voltage_step);

            get_retarded_se(parameters, boson);
            
            get_retarded_gf_boson(parameters, boson);
            
            get_greater_lesser_gf(parameters, boson, count);
        }

        get_difference(parameters, old_self_energy, fermion_up.retarded_se, difference);
        count++;
        std::cout << "The count is " << count << ". The difference is " << difference << std::endl;
    }

    test_retarded_gf(parameters, boson, fermion_up, fermion_down);


    z_prefactor =  1.0 / get_z_prefactor(parameters, boson, fermion_up, fermion_down);
}

