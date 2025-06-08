#include "numav.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    mod.add_type<Simulation>("Simulation")
        .constructor<>()
        .method("set_value", &Simulation::set_value)
        .method("run", &Simulation::run)
        .method("get_result", &Simulation::get_result);
}