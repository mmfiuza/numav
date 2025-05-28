
#include "jlcxx/jlcxx.hpp"



class Simulation {
public:
    Simulation(){};
    ~Simulation(){};
    void set_value(const double& value) {
        _value = value;
    };
    void run() {
        _value *= 2;
    };
    double get_result() {
        return _value;
    };
private:
    double _value;
};

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    mod.add_type<Simulation>("Simulation")
        .constructor<const double&>()
        .method("set_value", &Simulation::set_value)
        .method("run", %Simulation::run)
        .method("get_result", &Simulation::get_result);
}