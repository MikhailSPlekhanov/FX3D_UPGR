#pragma once

#include "../defines.hpp"
#include "../lbm.hpp"
#include "../shapes.hpp"
#include "../stat_funcs.hpp"

#ifdef TEMPERATURE
void NN_MPI();

void set_geometry(LBM& lbm);
#endif