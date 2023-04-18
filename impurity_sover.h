#include "parameters.h"
#include "slave_boson.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <limits>

void impurity_solver(Parameters &parameters, Interacting_GF &boson, Interacting_GF &fermion_up, Interacting_GF &fermion_down, int voltage_step);


void get_greater_lesser_se_fermion(Parameters &parameters, Interacting_GF &boson, Interacting_GF &fermion_up, int voltage_step);