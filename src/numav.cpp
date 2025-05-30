
#include "numav.hpp"

Simulation::Simulation(){};

Simulation::~Simulation(){};

void Simulation::set_value(const double& value) {
    _value = value;
};

void Simulation::run() {
    _value *= 2;
};

double Simulation::get_result() {
    return _value;
};
