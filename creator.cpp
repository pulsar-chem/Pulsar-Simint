#include "SimintERI.hpp"

using pulsar::modulemanager::ModuleCreationFuncs;


extern "C" {

ModuleCreationFuncs insert_supermodule(void)
{
    ModuleCreationFuncs cf;
    cf.add_cpp_creator<SimintERI>("SimintERI");
    return cf;
}



}

