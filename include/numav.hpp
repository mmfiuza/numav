#pragma once



class Simulation {
public:
    Simulation::Simulation();
    Simulation::~Simulation(){};
    void Simulation::set_value(const double& value);
    void Simulation::run();
    double Simulation::get_result();
private:
    double _value;
};
