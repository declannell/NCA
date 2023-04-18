#include "parameters.h"
#include "slave_boson.h"
#include "impurity_solver.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <limits>

void get_greater_lesser_se_fermion(Parameters &parameters, Interacting_GF &boson, Interacting_GF &fermion, int voltage_step) {
    
    dcomp prefactor = parameters.delta_energy * parameters.j1 * parameters.j1 * parameters.bandwidth * parameters.bandwidth * parameters.gamma / (2 * M_PI);
    for (int r = 0; r < parameters.steps; r++) {
        for (int i = 0; i < parameters.steps; i++) {
            double term = (parameters.energy.at(r) - parameters.energy.at(i)) * (parameters.energy.at(r) - parameters.energy.at(i))
                + parameters.bandwidth * parameters.bandwidth;

            fermion.greater_se.at(r) += (2 - fermi_function(parameters.energy.at(r) - parameters.energy.at(i) - parameters.voltage_l.at(i))
                                         - fermi_function(parameters.energy.at(r) - parameters.energy.at(i) - parameters.voltage_r.at(i)))
                                         * boson.greater_gf.at(i) / term;

            fermion.lesser_se.at(r) += (fermi_function(parameters.energy.at(r) - parameters.energy.at(i) - parameters.voltage_l.at(i))
                                         + fermi_function(parameters.energy.at(r) - parameters.energy.at(i) - parameters.voltage_r.at(i)))
                                         * boson.greater_gf.at(i) / term;
        }
        fermion.lesser_se.at(r) = fermion.lesser_se.at(r) * prefactor;
        fermion.greater_se.at(r) = fermion.greater_se.at(r) * prefactor;
    }
}



void impurity_solver(Parameters &parameters, Interacting_GF &boson, Interacting_GF &fermion_up, Interacting_GF &fermion_down, int voltage_step) {

    get_greater_lesser_se_fermion(parameters, boson, fermion_up, int voltage_step);
    get_greater_lesser_se_fermion(parameters, boson, fermion_down, int voltage_step);
}