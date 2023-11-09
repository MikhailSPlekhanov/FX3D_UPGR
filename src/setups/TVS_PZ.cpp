
#include "TVS_PZ.hpp"

/*void main_setup() {
	std::cout << (uint)(62) << std::endl;
}/**/

void TVS_PZ() {
	std::string working_path = get_exe_path(); std::string export_path = "export/rosatom/a2/";
	int MAX_THREADS = omp_get_max_threads(); omp_set_num_threads(MAX_THREADS);
	cout.precision(10);
	
	const uint L = 256u; 
	const float Re = 3699467.18046828; const float _sqrt3 = 1 / sqrt(3); const float sqrt3 = sqrt(3);

	const float L_real = 0.179; // реальный размер шестигранника, м
	const float u_real = 3.31486676868058; // реальная сокорость, м/c 
	const float real_density = 10373.886; // реальная плотность, кг/м3
	const float real_nu = 1.60390975956358 * 0.0000001; // реальная вязкость, м2/c
	const uint vnutr_diam = ceil((float)L / 2 * sqrt(3));
	const float vnutr_rad = (float)vnutr_diam / 2;

	const float uref = 0.1f; // скорость lbm, плотность 1
	const float h_to_w = 2.02562215807249; //  отношение высоты сегмента к диаметру описанной окружности шестиграника
	const float size = (float)L * h_to_w * 0.999f;
	const float vnesh_rad = (float)L / 2;
	const uint height = ceil(L * h_to_w);
	uint zazor = ceil((float)height*0.05f);	

	Coefs coefs(L_real, (float)L*sqrt(3)/2, u_real, uref, real_density);
	float nu = coefs.LB_nu(real_nu); 

	LBM lbm(vnutr_diam, L, height, nu);

	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
	const float3 center1 = float3(lbm.center().x, lbm.center().y, lbm.center().z );
	const float3 center = lbm.center();
	float xcent = lbm.center().x; float ycent = lbm.center().y; float zcent = lbm.center().z;
	float3 p = float3((float)Nx / 2, (float)Ny / 2, (float)Nz / 2);
	lbm.voxelize_stl(working_path + "../../stl/Small/Combined.stl", center1, size);
	//lbm.voxelize_stl("C:/Users/Mikhail/source/repos/stl/Small/HexTube.stl", center,  size);
	//lbm.voxelize_stl("C:/Users/Mikhail/source/repos/stl/Small/Combined.stl", center, size);

	Clock clock, clock2; clock2.start();
	uint time = 0; uint dt_check = 1000u; const uint dt_frame = 10u;

	std::vector<uint> highs = { 2u, Nz - 3u};
	//std::vector<uint> highs = { 0u, 1u, 2u, 3u, 5u };
	
	set_geometry(lbm, zazor);
	lbm.run(0u);
	//std::cout << "nu = " << nu << std::endl;
	for ( int i = 0; i < highs.size(); i++){
		std::cout << "Fluid cells in z = " << highs[i] << ": " << fluid_cells_in_plane(lbm, 'z', highs[i]) << std::endl;
	}
	while (lbm.get_t() < 80001u) {
		//проверка параметров по срезам

		if (lbm.get_t() % 1000u == 0) {
			std::cout << "step = " << lbm.get_t() << std::endl;
			std::vector<std::vector<float>> params(highs.size(), std::vector<float>(3, 0.0f)); int timeframe = 10; 
			for (int y = 0; y < timeframe; y++){
				lbm.update_fields(); lbm.rho.read_from_device(); lbm.u.read_from_device();
				for (int i = 0; i < highs.size(); i++) {
					params[i][0] += average_rhou_plane(lbm, 'z', highs[i], false);
					params[i][1] += average_rho_plane(lbm, 'z', highs[i], false);
					params[i][2] += average_u_plane(lbm, 'z', highs[i], false);
				}
				lbm.run(1u);
				time++;	
			}
			for (int i = 0; i < highs.size(); i++) {
				for (int j = 0; j < params[i].size(); j++) {
					params[i][j] /= timeframe;
				}
				std::cout << highs[i] << ": average rhou = " << params[i][0];// *fluid_cells_in_plane(lbm, 'z', highs[i]);
				std::cout << ", rhoav = " << params[i][1] << ", uav = " << params[i][2] << std::endl;
			} 
			std::cout << "Pressure drop = " << coefs.pressure(params[0][1] - params[highs.size() - 1][1]) << " Pa" << std::endl;
			std::cout << std::endl;
			
		}

		if (lbm.get_t() == 5050u) {
			ofstream logs;
			logs.open("logs.txt");
			lbm.rho.read_from_device(); lbm.u.read_from_device(); clock2.start();
			for (int i = 0; i < Nz; i++) {
				float rhoav = average_rho_plane(lbm, 'z', i, false); 
				float uav = average_u_plane(lbm, 'z', i, false);
				float rhouav = average_rhou_plane(lbm, 'z', i, false);
				int NF = fluid_cells_in_plane(lbm, 'z', i, false);
				logs << i << ", " << std::setprecision(10) << rhoav << ", " << uav << ", " << rhouav*NF << std::endl;
			}
			logs.close();
			std::cout << "written logs, el.time = " << clock2.stop() << std::endl;
		}

		//if (time % dt_frame == 0u) export_png(lbm, export_path);

		lbm.run(1u);
		time++;
	}
	write_file(get_exe_path() + "time.txt", print_time(clock.stop()));
} /**/

void set_geometry(LBM& lbm, float zazor) {
	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
	const float _sqrt3 = 1 / sqrt(3); 
	float xcent = lbm.center().x; float ycent = lbm.center().y; float zcent = lbm.center().z;
	uint L = lbm.get_Ny(); const float h_to_w = 2.02826361766759;
	float rhoin = 1.3f; float uin = 0.1f; 
	float rhoout = 1.05f;float uout = uin*rhoin/rhoout; 
	for (ulong n = 0ull; n < N; n++) {
		uint x = 0u, y = 0u, z = 0u; lbm.coordinates(n, x, y, z);
		//xrel = ((float)x - xcent) / vnutr_rad; yrel = ((float)y - ycent) / vnutr_rad;
		//u_in = (sq(xrel) + sq(yrel) < 1) ? 0.1f * (1 - sq(xrel) - sq(yrel)) : 0.0f;
		// ########################################################################### define geometry #############################################################################################
		//if (z < zazor || z > Nz - zazor - 1u) lbm.flags[n] = TYPE_F;
		if (x == 0u || x == Nx - 1u) lbm.flags[n] = TYPE_S;
		if ((x - xcent) * _sqrt3 + (y - ycent) > (float)L / h_to_w || (x - xcent) * _sqrt3 - (y - ycent) > (float)L / h_to_w || -(x - xcent) * _sqrt3 - (y - ycent) > (float)L / h_to_w || -(x - xcent) * _sqrt3 + (y - ycent) > (float)L / h_to_w) lbm.flags[n] = TYPE_S; // границы шестигранника
		//if (sphere(x, y, z, p, 30)) lbm.flags[n] = TYPE_S;
		//if (cylinder(x, y, z, p, float3(0.0f, 0.0f, 128.0f), 30.0f)) lbm.flags[n] = TYPE_S;
		if (z == Nz - 1u && lbm.flags[n] != TYPE_S) {
			int m = lbm.index(x, y, z - 1u);
			if (lbm.flags[m] == TYPE_S) { //ячейки над трубками
				lbm.flags[n] = TYPE_S;
			}
			else { //над свободным пространством
				lbm.flags[n] = TYPE_E;
				lbm.u.z[n] = uout; //lbm.u.z[n] = 0.0085f;
				lbm.rho[n] = rhoout; //lbm.rho[n] = 0.875f;
			}
		}
		if (z == 0u && lbm.flags[n] != TYPE_S) { // нулевой слой
			lbm.flags[n] = TYPE_E;
			lbm.u.z[n] = uin;//lbm.u.z[n] = 0.0085f;
			lbm.rho[n] = rhoin;//lbm.rho[n] = 1.261f;
		}
		//if (z == Nz - 1u && lbm.flags[n] != TYPE_S) {
		//	lbm.flags[n] = TYPE_S;
		//	//lbm.rho[n] = 1.0f; 
		//	lbm.u.z[n] = 0.15f;
		//}
		//if (z == Nz - 2u && lbm.flags[n] != TYPE_S) {
		//	lbm.flags[n] = TYPE_S;
		//	//lbm.rho[n] = 1.0f;
		//	lbm.u.z[n] = 0.15f;
		//}
		if (lbm.flags[n] != TYPE_S && lbm.flags[n] != TYPE_E) { // все ячейки жидкости
			lbm.flags[n] = TYPE_F;
			lbm.u.z[n] = uin + (uout - uin)/Nz * z; //lbm.u.z[n] = 0.0085f;
			lbm.rho[n] = rhoin + (rhoout-rhoin)/Nz * z; //lbm.rho[n] = 0.85f;
		}
	}
}

void Neumann(LBM& lbm, float drho) { //	условия неймана каждый шаг
	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
	float delta_rho; uint Ne2, layNz_1, layNz_2, lay0, lay1;
	if (true) {
		delta_rho = 0.0f; Ne2 = 0u;
		for (int i = 0; i < Nx; i++) { //вычет плотности на выходе
			for (int j = 0; j < Ny; j++) {
				layNz_1 = lbm.index(i, j, Nz - 1u);
				layNz_2 = lbm.index(i, j, Nz - 1u);
				if (lbm.flags[layNz_2] != TYPE_S && lbm.flags[layNz_1] != TYPE_S && lbm.rho[layNz_2] > drho) {
					//lbm.rho[lay0] = lbm.rho[lay1]; lbm.u.y[lay0] = lbm.u.y[lay1]; lbm.u.x[lay0] = lbm.u.x[lay1]; lbm.u.z[lay0] = lbm.u.z[lay1];
					lbm.rho[layNz_1] = lbm.rho[layNz_2] - drho;
					delta_rho += drho;
					Ne2 += 1u;
				}
			}
		}
		delta_rho /= (float)Ne2;
		for (int i = 0; i < Nx; i++) { //увеличение плотности на входе
			for (int j = 0; j < Ny; j++) {
				lay0 = lbm.index(i, j, 0u);
				lay1 = lbm.index(i, j, 1u);
				if (lbm.flags[lay0] != TYPE_S && lbm.flags[lay1] != TYPE_S) {
					//lbm.rho[layNz_1] = lbm.rho[layNz_2]; lbm.u.y[layNz_1] = lbm.u.y[layNz_2]; lbm.u.x[layNz_1] = lbm.u.x[layNz_2]; lbm.u.z[layNz_1] = lbm.u.z[layNz_2];
					lbm.rho[lay0] = lbm.rho[lay1] + delta_rho;
				}
			}
		}
	}
}

void rfr_podgon(LBM& lbm) {		//подгонка РФР по краям
	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
	uint layNz_1, layNz_2, lay0, lay1;
	for (int i = 0; i < Nx; i++) { 
		for (int j = 0; j < Ny; j++) {
			lay0 = lbm.index(i, j, 0u);
			lay1 = lbm.index(i, j, 1u);
			layNz_1 = lbm.index(i, j, Nz - 1u);
			layNz_2 = lbm.index(i, j, Nz - 2u);
			if (lbm.flags[layNz_2] != TYPE_S && lbm.u.z[layNz_1] < 0.1f) {
				//lbm.rho[lay0] = lbm.rho[lay1]; lbm.u.y[lay0] = lbm.u.y[lay1]; lbm.u.x[lay0] = lbm.u.x[lay1]; lbm.u.z[lay0] = lbm.u.z[lay1];
				lbm.u.z[layNz_1] += 0.0001f;
			}
			if (lbm.flags[lay1] != TYPE_S && lbm.u.z[lay0] < 0.1f) {
				//lbm.rho[layNz_1] = lbm.rho[layNz_2]; lbm.u.y[layNz_1] = lbm.u.y[layNz_2]; lbm.u.x[layNz_1] = lbm.u.x[layNz_2]; lbm.u.z[layNz_1] = lbm.u.z[layNz_2];
				lbm.u.z[lay0] += 0.0001f;
			}
		}
	}	
}

void kill_vortices(LBM& lbm) { //костыль против вихрей на входе и выходе
	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
	uint layNz_1, layNz_2, lay0, lay1, laysr;
	lbm.update_fields(); lbm.u.read_from_device(); lbm.rho.read_from_device();
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			lay0 = lbm.index(i, j, 0u);
			lay1 = lbm.index(i, j, 5u);
			laysr = lbm.index(i, j, Nz / 2u);
			layNz_2 = lbm.index(i, j, Nz - 6u);
			layNz_1 = lbm.index(i, j, Nz - 1u);
			if (lbm.flags[lay0] != TYPE_S && lbm.flags[layNz_1] != TYPE_S) {
				//lbm.u.y[lay0] = lbm.u.y[lay1]; lbm.u.x[lay0] = lbm.u.x[lay1]; 
				lbm.u.z[lay0] = lbm.u.z[laysr] * lbm.rho[layNz_2] / lbm.rho[lay1];
				lbm.rho[lay0] = lbm.rho[lay1];
				//lbm.u.y[layNz_1] = lbm.u.y[lay1]; lbm.u.x[layNz_1] = lbm.u.x[lay1]; 
				lbm.u.z[layNz_1] = lbm.u.z[laysr];
				lbm.rho[layNz_1] = lbm.rho[layNz_2];
			}
		}
	}
	lbm.u.write_to_device(); lbm.rho.write_to_device();
}
#ifdef GRAPHICS
void export_png(LBM& lbm, std::string export_path) {
	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
	//lbm.graphics.set_camera_free(float3(-0.5f * (float)Nx, -0.5f * (float)Ny, 0.55f * (float)Nz), -133.0f, 20.0f, 100.0f);
	//key_4 = true; key_2 = false;
	//lbm.graphics.write_frame_png(get_exe_path() + "export/rosatom/4a/");
	//key_2 = true; key_4 = false;
	//lbm.graphics.write_frame_png(get_exe_path() + "export/rosatom/2a/");
	//lbm.graphics.set_camera_free(float3(-1.73f * (float)Nx, -0.47f * (float)Ny, 0.42f * (float)Nz), -166.0f, 25.0f, 100.0f);
	//key_4 = true; key_2 = false;
	//lbm.graphics.write_frame_png(get_exe_path() + "export/rosatom/4b/");
	//key_2 = true; key_4 = false;
	//lbm.graphics.write_frame_png(get_exe_path() + "export/rosatom/2b/");
	//lbm.graphics.set_camera_free(float3(-0.64f * (float)Nx, -0.6f * (float)Ny, 0.6f * (float)Nz), -130.0f, 22.0f, 100.0f);
	//lbm.graphics.write_frame_png(get_exe_path() + "export/rosatom/test/MB01/"); 
	lbm.graphics.set_camera_free(float3(-0.8f * (float)Nx, -0.8f * (float)Ny, 0.55f * (float)Nz), -133.0f, 30.0f, 100.0f);
	key_4 = true; key_2 = false; key_1 = false;
	lbm.graphics.write_frame_png(get_exe_path() + "export/rosatom/a4/");
	key_4 = false; key_2 = true; key_1 = false;
	lbm.graphics.write_frame_png(get_exe_path() + "export/rosatom/a2/");
}
#endif
void cylinder_test(LBM& lbm) { //тест с циллиндром
	uint L = lbm.get_Ny(); const float sqrt3 = sqrt(3);
	const uint vnutr_diam = ceil((float)L / 2 * sqrt3);
	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
	float xcent = lbm.center().x; float ycent = lbm.center().y; float zcent = lbm.center().z; float xrel, yrel;
	for (ulong n = 0ull; n < N; n++) { 
		uint x = 0u, y = 0u, z = 0u; lbm.coordinates(n, x, y, z);
		if ( sq(x - xcent) + sq(y - ycent) > sq(vnutr_diam / 2) ) {
			lbm.flags[n] = TYPE_S;
		}
		else {
			xrel = (x - xcent) / vnutr_diam * 2; yrel = (y - ycent) / vnutr_diam * 2;
			if (z == 0u) {
				lbm.flags[n] = TYPE_S;
				lbm.u.z[n] = 0.1f;
			}
			else if (z == Nz - 1u) {
				lbm.flags[n] = TYPE_S;
				lbm.u.z[n] = 0.1f; 
				
			}
			else {
				lbm.flags[n] = TYPE_F;
				lbm.rho[n] = 1.0f;
				lbm.u.z[n] = (1 - sq(xrel) - sq(yrel)) * 0.15f; 
			}
		}
	}
}
