#pragma once

class Simulation {
public:
    Simulation();
    ~Simulation();
    void set_value(const double& value);
    void run();
    double get_result();
private:
    double _value;
};
