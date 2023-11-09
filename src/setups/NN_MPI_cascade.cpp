//#include "mpi.h"
#include "NN_MPI.hpp"

// Rayleigh-Benard convection

void NN_MPI() {
	std::cout.precision(4);
	uint NY = 128u;
	float gr = 0.001;  float Ra = 10000.0f; float Pr = 0.71f; float delta_x = (float)NY - 2.0f; float delta_t = sqrt(gr * delta_x);
	float nu = sqrt(Pr * Ra) * delta_t / sq(delta_x);
	float k = sqrt(1.0f / Pr / Ra) * delta_t / delta_x / delta_x;

	const uint Nx = 128u, Ny = 128u, Nz = 1u; const ulong N = Nx*Ny*Nz; 

	Clock clock1, clock2, clock3;
	int MAX_THREADS = omp_get_max_threads(); omp_set_num_threads(MAX_THREADS);

	std::vector<float> temp_matr(64);
	std::vector<float> heaters(8, 1.5f);
	float down_av, Nusselt, T_down_av; float T_top = 0.5f; 
	float fy = 0.0f;

	int size, rank = 0; int e = 0; int time = 0; float umax = 0.0f;

	std::cout << "starting communication... " << std::endl;
	//MPI_Init(nullptr, nullptr); MPI_Comm_size(MPI_COMM_WORLD, &size); MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Status status;
	
	std::string working_path =  get_exe_path();
	std::string export_path =  "NN/Train_control_shift/model_1507328/";
	std::string png_path = "photo/";	
	//ofstream heaters_logs;
	//heaters_logs.open("Heaters_logs.txt");
	if (rank == 0){
		e = 0;
		while (e < 1) {
			clock2.start();
			std::cout << (e+1) << " epoch simulation started " << std::endl;

			LBM lbm(128u, 128u, 1u, nu, 0.0f, 0.0f, 0.0f, 0.0f, 0.001f, 0.1f);
			set_geometry(lbm); fy = 0.0f; lbm.run(0); Nusselt = Rayleigh_Benard_Nusselt_2D(lbm, 1.5f, 0.5f, true);;
			while(lbm.get_t() < 10000u){
				if (lbm.get_fy() > -0.001){
					fy -= 0.0000001f;
					lbm.set_fy(fy);
				}
				/*if (e == 0 && lbm->get_t() % 500 == 0) {
					key_3 = true;
					(*lbm).graphics.set_camera_free(float3(0.0f,0.0f,0.89f*128.0f), -90.0f,90.0f,90.0f);
					(*lbm).graphics.write_frame_png(working_path + export_path + png_path);
				}*/
				/*if (e == 0 && lbm->get_t() % 100 == 0) {
					heaters_logs << (lbm->get_t()/100) << ", " << Nusselt;
					for (int j = 0; j < heaters.size(); j++) { heaters_logs << ", " << std::setprecision(10) << heaters[j];}
					heaters_logs << std::endl;
					Nusselt = Rayleigh_Benard_Nusselt_2D((*lbm), 1.5f, 0.5f, true);
				}*/
				lbm.run(1u);
			}

			//T_sensors_2D_to_vector((*lbm), temp_matr, 8, 8, 0, true);
			//MPI_Send(&temp_matr[0], 64, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);
			time = 1;
			while (time < 4097) {
				//MPI_Recv(&heaters[0], 8, MPI_FLOAT, 1, 1, MPI_COMM_WORLD, &status);
				//std::cout << "received heaters: "; for(int i = 0; i < 8; i++){std::cout << heaters[i] << ", ";}
				//std::cout << ", average = " << get_av(heaters) << std::endl;
				//shift_temps(heaters, 1.5f);
				//std::cout << "shifted heaters: "; for(int i = 0; i < 8; i++){std::cout << heaters[i] << ", ";}
				//std::cout << ", average = " << get_av(heaters) << std::endl;
				//set_complex_line_temp((*lbm), heaters, 1, 0, 'x');

				/*if (e == 0 && lbm->get_t() % 500 == 0) {
					key_3 = true;
					(*lbm).graphics.set_camera_free(float3(0.0f,0.0f,0.89f*128.0f), -90.0f,90.0f,90.0f);
					(*lbm).graphics.write_frame_png(working_path + export_path + png_path);
				}
				if (e == 0) {
					heaters_logs << (lbm->get_t()/100) << ", " << Nusselt;
					for (int j = 0; j < heaters.size(); j++) { heaters_logs << ", " << std::setprecision(10) << heaters[j];}
					heaters_logs << std::endl;
				}*/

				(lbm).run(100u); (lbm).u.read_from_device(); (lbm).T.read_from_device();
				
				//T_sensors_2D_to_vector((*lbm), temp_matr, 8, 8, 0, false);
				//MPI_Send(&temp_matr[0], 64, MPI_FLOAT, 1, 2, MPI_COMM_WORLD);

				umax = find_u_max_2D_omp((lbm), false);
				//umax = 0.1f;
				if (umax < 0.22f){
					//Nusselt = Rayleigh_Benard_Nusselt_2D((lbm), 1.5f, 0.5f, false);
					//MPI_Send(&Nusselt, 1, MPI_FLOAT, 1, 3, MPI_COMM_WORLD);
				}
				else { 
					//Nusselt = 1000000.0f;
					//MPI_Send(&Nusselt, 1, MPI_FLOAT, 1, 3, MPI_COMM_WORLD);
					std::cout << "Simulation diverged " << std::endl;
					break;
				}

				time++;
			}
			
			//if (e == 0) {
			//	key_3 = true;
			//	(*lbm).graphics.set_camera_free(float3(0.0f,0.0f,0.89f*128.0f), -90.0f,90.0f,90.0f);
			//	(*lbm).graphics.write_frame_png(get_exe_path() + png_path);
			//}
			std::cout << (e+1) << " epoch simulation is over, elapsed time = " << clock2.stop() << " s"  << std::endl;
			e++;
			//heaters_logs.close();
		}
	}
	else{std::cout << "c++ got wrong rank" << std::endl;}
	//MPI_Finalize();
}

void set_geometry(LBM& lbm) {
	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz(); uint x = 0u, y = 0u, z = 0u;
	float z_rel, y_rel, x_rel;
	for (ulong n = 0ull; n < N; n++) {
		lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		//z_rel = ((float)z - (float)Nz / 2.0f) / (float)Nz * 2.0f;
		y_rel = ((float)y - (float)Ny / 2.0f) / (float)Ny * 2.0f;
		x_rel = ((float)x - (float)Nx / 2.0f) / (float)Ny * 2.0f;
		if (y == 1u) {
			lbm.T[n] = 1.5f;
			lbm.flags[n] = TYPE_T;
		}
		else if (y == Ny - 2u) {
			lbm.T[n] = 0.5f;
			lbm.flags[n] = TYPE_T;
		}
		else if (y == 0u || y == Ny - 1u) {
			lbm.flags[n] = TYPE_S; // z non periodic
		}
		else if (x == 0u || x == Nx - 1u) lbm.flags[n] = TYPE_S;
		else {
			//lbm.T[n] = 1.0f - y_rel * 0.5f;
			//lbm.rho[n] = 1.001f - y_rel * 0.0001f;
			lbm.u.y[n] = 0.05f * sin(y_rel * 4 * pif) * cos(x_rel * 4 * pif);
			lbm.u.x[n] = -0.05f * cos(y_rel * 4 * pif) * sin(x_rel * 4 * pif);
			//lbm.u.y[n] = -0.03f * sin(y_rel * pif) * cos(x_rel * pif);// +random_symmetric(0.02f);
			//lbm.u.x[n] = 0.03f * cos(y_rel * pif) * sin(x_rel * pif);// +random_symmetric(0.02f);
			//lbm.u.x[n] = random_symmetric(0.02f);
			//lbm.u.y[n] = random_symmetric(0.02f);
			//lbm.u.z[n] = random_symmetric(0.005f);
		}
	}
}/**/

void shift_temps(std::vector<float>& temps, float desire){
	float av = 0.0f; float size = temps.size();
	for (int i = 0; i < size; i++){
		av += temps[i];
	}
	av /= size;
	float delta = av - desire;
	for (int i = 0; i < size; i++){
		temps[i] -= delta;
	}
}

float get_av(std::vector<float>& vec){
	float av = 0.0f;
	for(int i = 0; i < vec.size(); i++){
		av += vec[i];	
	}
	av /= vec.size();
	return av;
}
