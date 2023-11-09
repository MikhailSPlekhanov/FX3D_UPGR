#pragma once

#include "../units.hpp"
#include "../defines.hpp"
#include "../lbm.hpp"
#include "../shapes.hpp"
#include "../stat_funcs.hpp"
#include <fstream>

void TVS_PZ();

void set_geometry(LBM& lbm, float zazor);

void Neumann(LBM& lbm, float drho);

void kill_vortices(LBM& lbm);

void export_png(LBM& lbm, std::string export_path);