#pragma once

#include <iostream>
#include <functional>
#include <string.h>
#include "omp.h"

#include "defines.hpp"
#include "utilities.hpp"
#include "lbm.hpp"

using namespace std;

class LBM_state {
public:
	std::vector<float3> u;
	std::vector<float> rho;
#ifdef TEMPERATURE
	std::vector<float> T;
#endif
	int Nf;

public:
	LBM_state(LBM& lbm, const string us, const string rhos, const string Ts) {
		this->Nf = 0u;
		for (int i = 0; i < lbm.get_N(); i++) {
			if (lbm.flags[i] != TYPE_S) {
				this->Nf += 1u;
			}
		}
		if (us == "u") this->u.reserve(lbm.get_N());
		if (rhos == "rho") this->rho.reserve(lbm.get_N());
#ifdef TEMPERATURE
		if (Ts == "T") this->T.reserve(lbm.get_N());
#endif
	};
	LBM_state() = default;
	~LBM_state() = default;
	void copy_u_data(LBM& lbm) {
		lbm.u.read_from_device();
#pragma omp parallel for
		for (int i = 0; i < lbm.get_N(); i++) {
			this->u[i] = float3(lbm.u.x[i], lbm.u.y[i], lbm.u.z[i]);
		}
	}
	void write_u_data(LBM& lbm) {
#pragma omp parallel for
		for (int i = 0; i < lbm.get_N(); i++) {
			lbm.u.x[i] = this->u[i].x;
			lbm.u.y[i] = this->u[i].y;
			lbm.u.z[i] = this->u[i].z;
		}
		lbm.u.write_to_device();
	}
	void copy_rho_data(LBM& lbm) {
		lbm.rho.read_from_device();
#pragma omp parallel for
		for (int i = 0; i < lbm.get_N(); i++) {
			this->rho[i] = lbm.rho[i];
		}
	}
	void write_rho_data(LBM& lbm) {
#pragma omp parallel for
		for (int i = 0; i < lbm.get_N(); i++) {
			lbm.rho[i] = this->rho[i];
		}
		lbm.rho.write_to_device();
	}
#ifdef TEMPERATURE 
	void copy_T_data(LBM& lbm) {
		lbm.T.read_from_device();
#pragma omp parallel for
		for (int i = 0; i < lbm.get_N(); i++) {
			this->T[i] = lbm.T[i];
		}
	}
	void write_T_data(LBM& lbm) {
#pragma omp parallel for
		for (int i = 0; i < lbm.get_N(); i++) {
			lbm.T[i] = this->T[i];
		}
		lbm.T.write_to_device();
	}
#endif
};

float average_u_difference(LBM& new_lbm, LBM_state& old_lbm); //both should have the same dimensions

float find_u_max_2D(LBM& lbm);

float find_u_max_2D_omp(LBM& lbm, bool update = true);

float average_rho_plane(LBM& lbm, char plane, uint cnst, bool update = false);

float average_u_plane(LBM& lbm, char plane, uint cnst, bool update = false);

float average_rhou_plane_omp(LBM& lbm, char plane, uint cnst, bool update = false);

float average_rhou_plane(LBM& lbm, char plane, uint cnst, bool update = false);

float u_from_rhou_plane(LBM& lbm, char plane, uint cnst, bool update = false);

float integr_rhou_plane(LBM& lbm, char plane, uint cnst, bool update = false);

uint fluid_cells_in_plane(LBM& lbm, char plane, uint cnst, bool update = false);

void set_eq_u_plane(LBM& lbm, char plane, uint cnst, const float u);

#ifdef TEMPERATURE
std::vector<std::vector<std::vector<float>>>* temperature_matrix_3D(LBM& lbm, uint xdim, uint ydim, uint zdim);

void T_sensors_2D_to_vector(LBM& lbm, std::vector<float>& matr, int xdim, int ydim, uint z, bool update = true);

void T_sensors_2D_to_vector_omp(LBM& lbm, std::vector<float>& matr, int xdim, int ydim, uint z, bool update = true);

float Rayleigh_Benard_Nusselt(LBM& lbm, float down_av = 1.0f, bool update = true);

float Rayleigh_Benard_Nusselt_noomp(LBM& lbm);

float Rayleigh_Benard_Nusselt_2D(LBM& lbm, float T_down_av = 1.0f, float T_top = 1.0f, bool update = true);

float Rayleigh_Benard_Nusselt_2D_noomp(LBM& lbm, float T_down_av = 1.0f, float T_top = 1.0f, bool update = true);

void set_complex_plane_temp(LBM& lbm, std::vector<std::vector<float>>& temp, int cnst, char plane);

void set_complex_zplane_temp_alt(LBM& lbm, std::vector<float>& temp, int nx, int ny, int z);

void set_complex_zplane_temp_noomp(LBM& lbm, std::vector<std::vector<float>>& temp, int z);

void set_complex_line_temp(LBM& lbm, std::vector<float>& temp, int cnst1, int cnst2, char var);
#endif

void print_matr_from_pointer(float* begin, int Num, int rowlen);

class Index {
private:
	uint Nx, Ny, Nz;
public:
	Index(uint nx, uint ny, uint nz)
		:Nx{ nx }, Ny{ ny }, Nz{ nz } {}
	ulong operator ()(int x, int y, int z) {
		return (ulong)x + ((ulong)y + (ulong)z * (ulong)(this->Ny)) * (ulong)(this->Nx);
	}
};

inline ulong zindex(Index& index, int i, int j, int cnst) {
	return index(i, j, cnst);
}

inline ulong yindex(Index& index, int i, int j, int cnst) {
	return index(i, cnst, j);
}

inline ulong xindex(Index& index, int i, int j, int cnst) {
	return index(cnst, i, j);
}

inline ulong zlindex(LBM& lbm, int i, int j, int cnst) {
	return lbm.iindex(i, j, cnst);
}

inline ulong ylindex(LBM& lbm, int i, int j, int cnst) {
	return lbm.iindex(i, cnst, j);
}

inline ulong xlindex(LBM& lbm, int i, int j, int cnst) {
	return lbm.iindex(cnst, i, j);
}

inline float uz(LBM& lbm, ulong n) {
	return lbm.u.z[n];
}

inline float uy(LBM& lbm, ulong n) {
	return lbm.u.y[n];
}

inline float ux(LBM& lbm, ulong n) {
	return lbm.u.x[n];
}

