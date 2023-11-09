#include "stat_funcs.hpp"

float average_u_difference(LBM& new_lbm, LBM_state& old_lbm) { // both should have the same dimensions
	new_lbm.update_fields();
	new_lbm.u.read_from_device();
	float delta = 0.0f;
#pragma omp parallel for reduction (+:delta)
	for (int i = 0; i < new_lbm.get_N(); i++) {
		if (new_lbm.flags[i] != TYPE_S) {
#ifndef D2Q9
			delta += length(float3(new_lbm.u.x[i], new_lbm.u.y[i], new_lbm.u.z[i]) - old_lbm.u[i]);
#else
			delta += length(float3(new_lbm.u.x[i], new_lbm.u.y[i], 0.0f) - old_lbm.u[i]);
#endif //D2Q9
		}
	}

	delta /= (float)old_lbm.Nf;
	return delta;
}

float find_u_max_2D(LBM& lbm)
{
	lbm.update_fields();
	lbm.u.read_from_device();
	float umax = 0.0f, ucur;
	for (uint i = 0u; i < lbm.get_N(); i++) {
		ucur = sqrt(sq(lbm.u.x[i]) + sq(lbm.u.y[i]));
		if (ucur > umax) {
			umax = ucur;
		}
	}
	return umax;
}

float find_u_max_2D_omp(LBM& lbm, bool update)
{
	if (update) {
		lbm.update_fields();
		lbm.u.read_from_device();
	}
	int MAX_THREADS = omp_get_max_threads(); omp_set_num_threads(MAX_THREADS);
	std::vector<float> maxs(MAX_THREADS);
	float umax, ucur;
#pragma omp parallel shared (maxs) private (ucur)
	{
		int id = omp_get_thread_num();
		maxs[id] = 0.0f;
#pragma omp for
		for (int i = 0u; i < lbm.get_N(); i++) {
			ucur = sqrt(sq(lbm.u.x[i]) + sq(lbm.u.y[i]));
			if (ucur > maxs[id]) {
				maxs[id] = ucur;
			}
		}

#pragma omp flush(maxs)
#pragma omp master
		{
			umax = maxs[0];
			for (int i = 1; i < maxs.size(); i++) {
				if (maxs[i] > umax) {
					umax = maxs[i];
				}
			}
		}
	}
	return umax;
}

float average_rho_plane(LBM& lbm, char plane, uint cnst, bool update)
{
	if (update) {
		lbm.update_fields(); lbm.rho.read_from_device();
	}
	float rhoav = 0.0f; int Nfluid = 0;
	uint3 dim = lbm.get_Dim(); uint n; 
	uint idim = min(min(dim.x, dim.y), dim.z), jdim = idim;
	std::function<ulong(LBM&, int, int, int)> index;
	switch (plane) {
	case 'z': { index = zlindex; idim = dim.x; jdim = dim.y; break; }
	case 'y': { index = ylindex; idim = dim.x; jdim = dim.z; break; }
	case 'x': { index = xlindex; idim = dim.y; jdim = dim.z; break; }
	}
	for (int i = 0; i < idim; i++) {
		for (int j = 0; j < jdim; j++) {
			n = index(lbm, i, j, cnst);
			if (lbm.flags[n] != TYPE_S) {
				rhoav += lbm.rho[n];
				Nfluid += 1;
			}
		}
	}
	rhoav /= Nfluid;
	return rhoav;
}

float average_u_plane(LBM& lbm, char plane, uint cnst, bool update)
{
	if (update) {
		lbm.update_fields();
		lbm.u.read_from_device();
	}
	float uav = 0.0f; int Nfluid = 0;
	uint3 dim = lbm.get_Dim(); uint n;
	uint idim = min(min(dim.x, dim.y), dim.z), jdim = idim;
	std::function<ulong(LBM&, int, int, int)> index;
	std::function<float(LBM&, ulong)> uperp;
	switch (plane) {
	case 'z': { index = zlindex; uperp = uz; idim = dim.x; jdim = dim.y; break; }
	case 'y': { index = ylindex; uperp = uy; idim = dim.x; jdim = dim.z; break; }
	case 'x': { index = xlindex; uperp = ux; idim = dim.y; jdim = dim.z; break; }
	}
	for (int i = 0; i < idim; i++) {
		for (int j = 0; j < jdim; j++) {
			n = index(lbm, i, j, cnst);
			if (lbm.flags[n] != TYPE_S) {
				uav += uperp(lbm, n);
				Nfluid += 1;
			}
		}
	}
	uav /= Nfluid;
	return uav;
}

float average_rhou_plane_omp(LBM& lbm, char plane, uint cnst, bool update)
{
	if (update) {
		lbm.update_fields(); lbm.rho.read_from_device(); lbm.u.read_from_device();
	}
	float rhouav = 0.0f; int Nfluid = 0;
	uint3 dim = lbm.get_Dim(); uint n; Index ind(dim.x, dim.y, dim.z);
	uint idim = min(min(dim.x, dim.y), dim.z), jdim = idim;
	std::function<ulong(Index&, int, int, int)> index;
	std::function<float(LBM&, ulong)> uperp;
	switch (plane) {
	case 'z': { index = zindex; uperp = uz; idim = dim.x; jdim = dim.y; break; }
	case 'y': { index = yindex; uperp = uy; idim = dim.x; jdim = dim.z; break; }
	case 'x': { index = xindex; uperp = ux; idim = dim.y; jdim = dim.z; break; }
	}
	int MAX_THREADS = omp_get_max_threads(); omp_set_num_threads(MAX_THREADS);
	std::vector<float> rhous(MAX_THREADS); 
	std::vector<int> NFs(MAX_THREADS);
#pragma omp parallel shared (rhous, NFs) firstprivate (idim, jdim, cnst, ind) private (n)
	{
		int id = omp_get_thread_num();
		rhous[id] = 0.0f; NFs[id] = 0;
#pragma omp for
		for (int i = 0; i < idim; i++) {
			for (int j = 0; j < jdim; j++) {
				n = index(ind, j, i, cnst); // выгоднее в случае x,y-плоскости, т к в памяти рядом лежат строчки по x 
				if (lbm.flags[n] != TYPE_S) {
					rhous[id] += lbm.rho[n] * uperp(lbm, n);
					NFs[id] += 1;
				}
			}
		}
#pragma omp flush(rhous, NFs)
#pragma omp master
		{
			for (int i = 1; i < rhous.size(); i++) {
				rhouav += rhous[i];
				Nfluid += NFs[i];
			}
		}
	}
	rhouav /= Nfluid;
	return rhouav;
}

float average_rhou_plane(LBM& lbm, char plane, uint cnst, bool update)
{
	if (update) {
		lbm.update_fields(); lbm.rho.read_from_device(); lbm.u.read_from_device();
	}
	float rhouav = 0.0f; int Nfluid = 0;
	uint3 dim = lbm.get_Dim(); uint n; Index ind(dim.x, dim.y, dim.z);
	uint idim = min(min(dim.x, dim.y), dim.z), jdim = idim;
	std::function<ulong(Index&, int, int, int)> index;
	std::function<float(LBM&, ulong)> uperp;
	switch (plane) {
	case 'z': { index = zindex; uperp = uz; idim = dim.x; jdim = dim.y; break; }
	case 'y': { index = yindex; uperp = uy; idim = dim.x; jdim = dim.z; break; }
	case 'x': { index = xindex; uperp = ux; idim = dim.y; jdim = dim.z; break; }
	}
	for (int i = 0; i < idim; i++) {
		for (int j = 0; j < jdim; j++) {
			n = index(ind, i, j, cnst); // выгоднее в случае x,y-плоскости, т к в памяти рядом лежат строчки по x 
			if (lbm.flags[n] != TYPE_S) {
				rhouav += lbm.rho[n] * uperp(lbm, n);
				Nfluid += 1;
			}
		}
	}
	rhouav /= Nfluid;
	return rhouav;
}

float u_from_rhou_plane(LBM& lbm, char plane, uint cnst, bool update)
{
	if (update) { lbm.update_fields(); lbm.rho.read_from_device(); lbm.u.read_from_device(); }
	float rhou = 0.0f; float rho = 0.0f;
	uint3 dim = lbm.get_Dim(); uint n; 
	Index ind(dim.x, dim.y, dim.z);
	uint idim = min(min(dim.x, dim.y), dim.z), jdim = min(min(dim.x, dim.y), dim.z);
	std::function<ulong(Index&, int, int, int)> index;
	std::function<float(LBM&, ulong)> uperp;
	switch (plane) {
	case 'z': { index = zindex; uperp = uz; idim = dim.x; jdim = dim.y; break; }
	case 'y': { index = yindex; uperp = uy; idim = dim.x; jdim = dim.z; break; }
	case 'x': { index = xindex; uperp = ux; idim = dim.y; jdim = dim.z; break; }
	}
//#pragma omp parallel for reduction (+ : rhoav, Nfluid)
	for (int i = 0; i < idim; i++) {
		for (int j = 0; j < jdim; j++) {
//#pragma omp atomic read
			n = index(ind, i, j, cnst);
			if (lbm.flags[n] != TYPE_S) {
				rhou += lbm.rho[n] * uperp(lbm, n);
				rho += lbm.rho[n];
			}
		}
	}
	return (rhou / rho);
}

float integr_rhou_plane(LBM& lbm, char plane, uint cnst, bool update)
{
	if (update) { lbm.update_fields(); lbm.rho.read_from_device(); lbm.u.read_from_device(); }
	float rhou = 0.0f; 
	uint3 dim = lbm.get_Dim(); uint n; Index ind(dim.x, dim.y, dim.z);
	uint idim = min(min(dim.x, dim.y), dim.z), jdim = min(min(dim.x, dim.y), dim.z);
	std::function<ulong(Index&, int, int, int)> index;
	std::function<float(LBM&, ulong)> uperp;
	switch (plane) {
	case 'z': { index = zindex; uperp = uz; idim = dim.x; jdim = dim.y; break; }
	case 'y': { index = yindex; uperp = uy; idim = dim.x; jdim = dim.z; break; }
	case 'x': { index = xindex; uperp = ux; idim = dim.y; jdim = dim.z; break; }
	}
	//#pragma omp parallel for reduction (+ : rhoav, Nfluid)
	for (int i = 0; i < idim; i++) {
		for (int j = 0; j < jdim; j++) {
			//#pragma omp atomic read
			n = index(ind, i, j, cnst);
			if (lbm.flags[n] != TYPE_S) {
				rhou += lbm.rho[n] * uperp(lbm, n);
			}
		}
	}
	return (rhou);
}

uint fluid_cells_in_plane(LBM& lbm, char plane, uint cnst, bool update)
{
	if (update) {
		lbm.update_fields(); lbm.rho.read_from_device();lbm.u.read_from_device();
	}
	uint3 dim = lbm.get_Dim(); uint n;
	uint idim = min(min(dim.x, dim.y), dim.z), jdim = min(min(dim.x, dim.y), dim.z);
	std::function<ulong(LBM&, int, int, int)> index;
	switch (plane) {
	case 'z': { index = zlindex; idim = dim.x; jdim = dim.y; break; }
	case 'y': { index = ylindex; idim = dim.x; jdim = dim.z; break; }
	case 'x': { index = xlindex; idim = dim.y; jdim = dim.z; break; }
	}
	uint count = 0u; 
	int MAX_THREADS = omp_get_max_threads(); omp_set_num_threads(MAX_THREADS);
	std::vector<float> cnts(MAX_THREADS);
#pragma omp parallel shared (cnts) firstprivate(idim, jdim, cnst) private (n)
	{
		int id = omp_get_thread_num();
		cnts[id] = 0.0f;
#pragma omp for
		for (int i = 0; i < idim; i++) {
			for (int j = 0; j < jdim; j++) {
				n = index(lbm, i, j, cnst);
				if (lbm.flags[n] != TYPE_S) {
					cnts[id]++;
				}
			}
		}
#pragma omp flush(cnts)
#pragma omp master
		{
			for (int i = 1; i < cnts.size(); i++) {
				count += cnts[i];
			}
		}
	}
	return count;
}

void set_eq_u_plane(LBM& lbm, char plane, uint cnst, const float u)
{
	uint3 dim = lbm.get_Dim(); int N = lbm.get_N();
	switch (plane) {
	case 'z': {
		//#pragma omp parallel for reduction (+ : rhoav, Nfluid)
		for (int i = 0; i < dim.x; i++) {
			for (uint j = 0; j < dim.y; j++) {
				//#pragma omp atomic read
				uint n = lbm.index((uint)i, j, cnst);
				if (lbm.flags[n] == TYPE_E) {
					lbm.u.z[n] = u;
				}
			}
		}
		break;
	}
	case 'y': {
		//#pragma omp parallel for reduction (+ : rhoav, Nfluid)
		for (int i = 0; i < dim.x; i++) {
			for (uint j = 0; j < dim.z; j++) {
				//#pragma omp atomic read
				uint n = lbm.index((uint)i, cnst, j);
				if (lbm.flags[n] == TYPE_E) {
					lbm.u.y[n] = u;

				}
			}
		}
		break;
	}
	case 'x': {
		//#pragma omp parallel for reduction (+ : rhoav, Nfluid)
		for (int i = 0; i < dim.y; i++) {
			for (uint j = 0; j < dim.z; j++) {
				//#pragma omp atomic read
				uint n = lbm.index(cnst, (uint)i, j);
				if (lbm.flags[n] == TYPE_E) {
					lbm.u.x[n] = u;
				}
			}
		}
		break;
	}
	}
	lbm.u.write_to_device();
}

#ifdef TEMPERATURE

std::vector<std::vector<std::vector<float>>>* temperature_matrix_3D(LBM& lbm, uint xdim, uint ydim, uint zdim) { //хреновая реализация
	lbm.update_fields();
	lbm.T.read_from_device();
	uint Nx = lbm.get_Nx(); uint Ny = lbm.get_Ny(); uint Nz = lbm.get_Nz(); uint N = lbm.get_N();
	std::vector<std::vector<std::vector<float>>> res(Nx);
	for (uint n = 0; n < Nx; n++) res[n].reserve(Ny);
	for (uint n = 0; n < Ny; n++) {
		for (uint m = 0; m < Ny; m++) {
			res[n][m].reserve(Nz);
		}
	}
	int xstep = Nx / xdim; int ystep = Ny / ydim; int zstep = Nz / zdim;
	int n;
	int i = 0, j = 0, k = 0;
	//#pragma omp parallel for
	for (int x = xstep / 2u; x < Nx; x += xstep) {
		for (int y = ystep / 2u; y < Ny; y += ystep) {
			for (int z = zstep / 2u; z < Nz; z += zstep) {
				n = lbm.index(x, y, z);
				if (lbm.flags[n] != TYPE_S) { res[i][j][k] = lbm.T[n]; }
				k++;
			}
			j++;
		}
		i++;
	}
	return &res;
}

void T_sensors_2D_to_vector(LBM& lbm, std::vector<float>& matr, int xdim, int ydim, uint z, bool update)
{
	if (update) { lbm.update_fields(); lbm.T.read_from_device();	}
	uint3 dim = lbm.get_Dim(); ulong N = lbm.get_N();
	uint xstep = dim.x / xdim; uint ystep = dim.y / ydim;
	for (int i = 0; i < xdim; i++) {
		for (int j = 0; j < ydim; j++) {
			//std::cout << "reading cell x = " << (2 * i + 1) * xstep / 2u << ", y = " << (2 * j + 1) * ystep / 2u << ", temp = " << lbm.T[lbm.index((2 * i + 1) * xstep / 2u, (2 * j + 1) * ystep / 2u, z)] << std::endl;
			matr[i + j * ydim] = lbm.T[lbm.index((2 * i + 1) * xstep / 2u, (2 * j + 1) * ystep / 2u, z)];
		}
	}
}

void T_sensors_2D_to_vector_omp(LBM& lbm, std::vector<float>& matr, int xdim, int ydim, uint z, bool update)
{
	if (update) {
		lbm.update_fields(); lbm.T.read_from_device();
	}
	int MAX = omp_get_max_threads(); omp_set_num_threads(MAX);
	uint3 dim = lbm.get_Dim(); ulong N = lbm.get_N();
	uint xstep = dim.x / xdim; uint ystep = dim.y / ydim;
#pragma omp parallel for firstprivate(xstep, ystep, xdim, ydim, z)
	for (int i = 0; i < xdim; i++) {
		for (int j = 0; j < ydim; j++) {
			//std::cout << "i = " << i << ", j = " << j << " reading cell x = " << (2 * i + 1) * xstep / 2u << ", y = " << (2 * j + 1) * ystep / 2u << ", temp = " << lbm.T[lbm.index((2 * i + 1) * xstep / 2u, (2 * j + 1) * ystep / 2u, z)] << std::endl;
			matr[i + j * ydim] = lbm.T[lbm.index((2 * i + 1) * xstep / 2u, (2 * j + 1) * ystep / 2u, z)];
		}
	}
}

float Rayleigh_Benard_Nusselt(LBM& lbm, float down_av, bool update) // можно еще работать над ускорением
{
	if (update) {
		lbm.update_fields();
		lbm.u.read_from_device();
		lbm.T.read_from_device();
	}
	float nus = 0.0f; float plane_var = 0.0f; 
	int Nx = lbm.get_Nx(); int Ny = lbm.get_Ny(); int Nz = lbm.get_Nz(); int N = lbm.get_N();
#pragma omp parallel for reduction (+ : nus)
	for (int z = 2; z < Nz - 2; z++) {
		for (int x = 0; x < Nx; x++) {
			for (int y = 0; y < Ny; y++) {
				/*int n, nd, nu;
				n = lbm.index(x, y, z);
				nd = lbm.index(x, y, z - 1);
				nu = lbm.index(x, y, z + 1);*/
				nus += lbm.u.z[lbm.index(x, y, z)] * lbm.T[lbm.index(x, y, z)] - (lbm.T[lbm.index(x, y, z + 1)] - lbm.T[lbm.index(x, y, z - 1)]) / 2; //вроде так лучше, если не создавать переменных индексов
			}
		}
	}
	return nus / (float)((Nz - 2) * Nx * Ny);
}

float Rayleigh_Benard_Nusselt_noomp(LBM& lbm) // можно еще работать над ускорением
{
	lbm.update_fields();
	lbm.u.read_from_device();
	lbm.T.read_from_device();
	float nus = 0.0f; float plane_var = 0.0f; int n, nd, nu;
	int Nx = lbm.get_Nx(); int Ny = lbm.get_Ny(); int Nz = lbm.get_Nz(); int N = lbm.get_N();
	for (int z = 2; z < Nz - 2; z++) {
		for (int x = 0; x < Nx; x++) {
			for (int y = 0; y < Ny; y++) {
				n = lbm.index(x, y, z);
				nd = lbm.index(x, y, z - 1);
				nu = lbm.index(x, y, z + 1);
				nus += lbm.u.z[n] * lbm.T[n] - (lbm.T[nu] - lbm.T[nd]) / 2;
			}
		}
	}
	return nus / (float)((Nz - 4) * Nx * Ny);
}

float Rayleigh_Benard_Nusselt_2D(LBM& lbm, float T_down_av, float T_top, bool update)
{
	if (update) {
		lbm.update_fields(); lbm.u.read_from_device(); lbm.T.read_from_device();
	}
	float nus = 0.0f;
	uint3 dim = lbm.get_Dim(); ulong N = lbm.get_N(); float delta0 = T_down_av - T_top;
	float delta_T = (delta0 < -0.1f || delta0 > 0.1f) ? delta0 : (delta0 > 0) ? 0.1f : -0.1f;
	uint n = 0, nu = 0, nd = 0; float alpha;

	int MAX_THREADS = omp_get_max_threads(); omp_set_num_threads(MAX_THREADS);
	std::vector<float> maxs(MAX_THREADS);

#pragma omp parallel shared (maxs) private (n, nu, nd, alpha)
	{
		alpha = lbm.get_alpha();
		int id = omp_get_thread_num();
		maxs[id] = 0.0f;
#pragma omp for
		for (int x = 0; x < dim.x; x++) {
			for (int y = 1; y < dim.y - 1u; y++) {
				nu = lbm.index(x, y + 1, 0);
				nd = lbm.index(x, y - 1, 0);
				n = lbm.index(x, y, 0);
				maxs[id] += lbm.u.y[n] * (lbm.T[n] - 1.0f) - alpha * (lbm.T[nu] - lbm.T[nd]) / 2.0f;
				//maxs[id] += lbm.u.y[n] * lbm.T[n];
			}
		}

#pragma omp flush(maxs)
#pragma omp master
		{
			for (int i = 0; i < maxs.size(); i++) {
				nus += maxs[id];
			}
		}
	}
	return (nus / lbm.get_alpha() / delta_T / (float)dim.x);
	//return (1.0f +  nus / lbm.get_alpha() / delta_T / (float)dim.x );
}

float Rayleigh_Benard_Nusselt_2D_noomp(LBM& lbm, float T_down_av, float T_top, bool update)
{
	if (update) {
		lbm.update_fields(); lbm.u.read_from_device(); lbm.T.read_from_device();
	}
	float nus = 0.0f; 
	uint3 dim = lbm.get_Dim(); ulong N = lbm.get_N();
	float delta_T = loc::clamp(T_down_av - T_top, -0.1, 0.1);
	std::cout << "T_down_av = T_top = " << (T_down_av - T_top) << ", delta_T = " << delta_T << std::endl;
	uint n = 0;
	for (int x = 0; x < dim.x; x++) {
		for (int y = 1; y < dim.y - 1u; y++) {
			n = lbm.index(x, y, 0);
			nus += lbm.u.y[n] * lbm.T[n];
		}
	}
	return (1 + nus / (float)(dim.x) / lbm.get_alpha() / delta_T);
}

void set_complex_plane_temp(LBM& lbm, std::vector<std::vector<float>>& temp, int cnst, char plane)
{
	//lbm.update_fields();
	//lbm.T.read_from_device();
	uint3 dim = lbm.get_Dim(); int N = lbm.get_N();
	int i = 0; int j = 0; uint n = 0u;
	switch (plane) {
	case'z':{
		int xsize = temp.size(); int ysize = temp[0].size();
		int xchunk = dim.x / xsize; int ychunk = dim.y / ysize;
		for (int x = 0; x < dim.x; x++) {
			for (int y = 0; y < dim.y; y++) {
				i = x / xchunk; j = y / ychunk; n = lbm.index(x, y, cnst);
				if (lbm.flags[n] == TYPE_T)
					lbm.T[n] = temp[i][j];
			}
		}
		break;
	}
	case'y': {
		int xsize = temp.size(); int zsize = temp[0].size();
		int xchunk = dim.x / xsize; int zchunk = dim.z / zsize;
		for (int x = 0; x < dim.x; x++) {
			for (int z = 0; z < dim.z; z++) {
				i = x / xchunk; j = z / zchunk; n = lbm.index(x, cnst, z);
				if (lbm.flags[n] == TYPE_T)
					lbm.T[n] = temp[i][j];
			}
		}
		break;
	}
	case'x': {
		int ysize = temp.size(); int zsize = temp[0].size();
		int ychunk = dim.y / ysize; int zchunk = dim.z / zsize;
		for (int y = 0; y < dim.y; y++) {
			for (int z = 0; z < dim.z; z++) {
				i = y / ychunk; j = z / zchunk; n = lbm.index(cnst, y, z);
				if (lbm.flags[n] == TYPE_T)
					lbm.T[n] = temp[i][j];
			}
		}
		break;
	}
	}

	lbm.T.write_to_device();
}

void set_complex_zplane_temp_alt(LBM& lbm, std::vector<float>& temp, int nx, int ny, int z)
{
//	lbm.update_fields();
//	lbm.T.read_from_device();
//	int Nx = lbm.get_Nx(); int Ny = lbm.get_Ny(); int Nz = lbm.get_Nz(); int N = lbm.get_N();
//	int xsize = temp.size(); int ysize = temp[0].size();
//	int xchunk = Nx / xsize; int ychunk = Ny / ysize;
//#pragma omp parallel for
//	for (int x = 0; x < Nx; x++) {
//		for (int y = 0; y < Ny; y++) {
//			float t; int i = x / xchunk; int j = y / ychunk;
//			//#pragma omp atomic read
//			t = temp[i][j];
//			lbm.T[lbm.index(x, y, z)] = t;
//		}
//	}
//	lbm.T.write_to_device();
}

void set_complex_zplane_temp_noomp(LBM& lbm, std::vector<std::vector<float>>& temp, int z) //пока лучше без omp
{
	lbm.update_fields();
	lbm.T.read_from_device();
	int Nx = lbm.get_Nx(); int Ny = lbm.get_Ny(); int Nz = lbm.get_Nz(); int N = lbm.get_N();
	int xsize = temp.size(); int ysize = temp[0].size();
	int xchunk = Nx / xsize; int ychunk = Ny / ysize;

	for (int x = 0; x < Nx; x++) {
		for (int y = 0; y < Ny; y++) {
			lbm.T[lbm.index(x, y, z)] = temp[x / xchunk][y / ychunk];
		}
	}
	lbm.T.write_to_device();
}

void set_complex_line_temp(LBM& lbm, std::vector<float>& temp, int cnst1, int cnst2, char var)
{
	uint3 dim = lbm.get_Dim(); ulong N = lbm.get_N();
	uint n = 0u;
	lbm.T.read_from_device();
	switch (var) {
	case'z': {
		int zsize = temp.size();
		int zchunk = dim.z / zsize; 
		for (int z = 0; z < dim.z; z++) {
			n = lbm.index(cnst1, cnst2, z);
			if (lbm.flags[n] == TYPE_T)
				lbm.T[n] = temp[z / zchunk];
		}
		break;
	}
	case'y': {
		int ysize = temp.size();
		int ychunk = dim.y / ysize;
		for (int y = 0; y < dim.y; y++) {
			n = lbm.index(cnst1, y, cnst2);
			if (lbm.flags[n] == TYPE_T)
				lbm.T[n] = temp[y / ychunk];
		}
		break;
	}
	case'x': {
		int xsize = temp.size();
		int xchunk = dim.x / xsize;
		for (int x = 0; x < dim.x; x++) {
			n = lbm.index(x, cnst1, cnst2);
			if (lbm.flags[n] == TYPE_T)
				lbm.T[n] = temp[x / xchunk];
		}
		break;
	}
	}
	lbm.T.write_to_device();
}
#endif

void print_matr_from_pointer(float* begin, int Num, int rowlen)
{
	for (int i = 0; i < Num / rowlen; i++) {
		for (int j = 0; j < rowlen; j++) {
			std::cout << *(begin + i * rowlen + j) << " ";
		}
		std::cout << std::endl;
	}
}