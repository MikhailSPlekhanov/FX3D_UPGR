//#include "mpi.h"
#include "NN_MPI.hpp"
#include <chrono>
#include <thread>

// Rayleigh-Benard convection
#ifdef TEMPERATURE
void NN_MPI() {
	using namespace std::this_thread;     // sleep_for, sleep_until
	using namespace std::chrono_literals;	// ns, us, ms, s, h, etc.
	using std::chrono::system_clock;
	std::cout.precision(7);
	uint NY = 128u;
	float gr = 0.0001;  float Ra = 10000.0f; float Pr = 0.71f; float delta_x = (float)NY - 2.0f; float delta_t = sqrt(gr * delta_x);
	float nu = sqrt(Pr * Ra) * delta_t / sq(delta_x);
	float k = sqrt(1.0f / Pr / Ra) * delta_t / delta_x / delta_x;
	float tau = 0.56f; 
	//float nu2 = (tau - 0.5f)/ 3.0f;
	float tauT = (tau - 0.5f) / Pr + 0.5; 
	//float efgr = Ra * alpha * nu2 / float(NY) / float(NY) / (float)NY;
	float beta = 0.1f;//float beta = efgr / gr;
	float delta_T = 1.0f;
	float nu2 = sqrt(Pr * beta * gr * delta_T * cb(delta_x) / Ra);
	float alpha = nu2 / Pr; // float alpha = (2.0f * tauT - 1.0f) / 4.0f;

	const uint Nx = NY, Ny = NY, Nz = 1u; const ulong N = Nx*Ny*Nz; 

	uint time = 0u;
	Clock clock1, clock2, clock3;
	int MAX_THREADS = omp_get_max_threads(); omp_set_num_threads(MAX_THREADS);

	std::vector<float> temp_matr(64);
	std::vector<float> heaters(8, 0.0f);
	float down_av, Nusselt, T_down_av; float T_top = 0.5f; 

	int size, rank; uint epoch = 0u;
	float fy = 0.0f;

	//MPI_Init(nullptr, nullptr);
	//MPI_Comm_size(MPI_COMM_WORLD, &size);
	//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//MPI_Status status;
	float umax1, umax2; int nrun = 0; uint e = 0u;
	
	while (epoch < 1u) {
		//LBM lbm(NY, NY, 1u, nu, 0.0f, 0.0f, 0.0f, 0.0f, 0.001, 0.1);
		LBM lbm(NY, NY, 1u, nu2, 0.0f, 0.0f, 0.0f, 0.0f, alpha, beta);
		std::cout << "Ra = " << 0.001 * lbm.get_beta() * lbm.get_Ny() * lbm.get_Ny() * lbm.get_Ny() / nu / lbm.get_alpha() << std::endl;
		set_geometry(lbm);
		fy = 0.0f;
		
		//T_sensors_2D_to_vector(lbm, temp_matr, 8, 8, 0, true);
		//std::cout << "time = " << lbm.get_t() << ":" << std::endl;
		//print_matr_from_pointer(&temp_matr[0], 64, 8);
		int pretime = 10000;
		float dgr = gr / (float)(pretime);
		while ((lbm).get_t() < pretime) {
			if ((lbm).get_fy() > -gr) {
				fy -= dgr;
				(lbm).set_fy(fy);
			}
			/*if (lbm.get_t() % 500 == 0) {
				key_3 = true;
				lbm.graphics.set_camera_free(float3(0.0f, 0.0f, 0.89f * (float)Nx), -90.0f, 90.0f, 90.0f);
				lbm.graphics.write_frame_png("C:/Users/Mikhail/source/repos/export/2D_RaBen_NN_T/");
			}*/
			(lbm).run(1u);
			//if (lbm.get_t() % 1000 == 0) {
			//	Nusselt = Rayleigh_Benard_Nusselt_2D(lbm, 1.5f, 0.5f, true);
			//	std::cout << "Nu = " << Nusselt << std::endl;
			//}
		}
		//key_P = false;
		for (int i = 0; i < 2048000; i++) {
			//T_sensors_2D_to_vector(lbm, temp_matr, 8, 8, 0, true);
			//std::cout << "time = " << lbm.get_t() << ":" << std::endl;
			//print_matr_from_pointer(&temp_matr[0], 64, 8);
			//key_P = false;
			/*if (lbm.get_t() % 500 == 0) {
				key_3 = true;
				lbm.graphics.set_camera_free(float3(0.0f, 0.0f, 0.89f * (float)Nx), -90.0f, 90.0f, 90.0f);
				lbm.graphics.write_frame_png("C:/Users/Mikhail/source/repos/export/2D_RaBen_NN_T/");
			}*/
			lbm.run(2000u);
			key_P = 0;
			for (int i = 0; i < 8; i++) { heaters[i] = 1.5f + random_symmetric(1.0f); }
			set_complex_line_temp((lbm), heaters, 1, 0, 'x');
			std::cout << find_u_max_2D(lbm) << std::endl;
			//sleep_for(1ms);
			//sleep_until(system_clock::now() + 1s);
			//key_P = false;
		}
		std::cout << epoch << " epoch is over" << std::endl<<std::endl;
		epoch++;
		//delete lbm;
	}
	write_file(get_exe_path() + "time.txt", print_time(clock1.stop()));
}

void set_geometry(LBM& lbm) {
	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz(); uint x = 0u, y = 0u, z = 0u;
	float z_rel, y_rel, x_rel, phi, r;
	for (ulong n = 0ull; n < N; n++) {
		lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		//z_rel = ((float)z - (float)Nz / 2.0f) / (float)Nz * 2.0f;
		y_rel = ((float)y - (float)Ny / 2.0f) / (float)Ny * 2.0f;
		x_rel = ((float)x - (float)Nx / 2.0f) / (float)Ny * 2.0f;
		r = sqrt(sq(y_rel) + sq(x_rel)); phi = atan2(y_rel, x_rel);
		if (y == 1u) {
			lbm.flags[n] = TYPE_T;
			lbm.T[n] = 1.5f;// 1.0f - cos(pif * x_rel);
		}
		else if (y == Ny - 2u) {
			lbm.flags[n] = TYPE_T;
			lbm.T[n] = 0.5f;// -0.5f * cos(pif * x_rel);
		}
		
		else if (y == 0u || y == Ny - 1u || x == 0u || x == Nx - 1u) {
			lbm.flags[n] = TYPE_S;
		}
		else {

			//lbm.T[n] = 1.5 - y_rel * 0.5f;
			//lbm.rho[n] = 1.001f - y_rel * 0.0001f;
			lbm.u.y[n] = r < 1.0f ? -0.05f * cos(phi) * r * (r - 1.0f) : 0.0f;
			lbm.u.x[n] = r < 1.0f ? 0.05f * sin(phi) * r * (r - 1.0f) : 0.0f;
			//lbm.u.y[n] = -0.01f * sin((y_rel-1) * pif /2) * cos((x_rel-1) * pif /2 );
			//lbm.u.x[n] = 0.01f * cos((y_rel-1) * pif /2 ) * sin((x_rel-1) * pif /2);
			//lbm.u.y[n] = -0.05f * sin(y_rel * pif/2) * cos(x_rel  * pif/2);
			//lbm.u.x[n] = 0.05f * cos(y_rel *pif/2) * sin(x_rel * pif/2);
			//lbm.u.y[n] = -0.05f * sin(10*y_rel * pif) * cos(10*x_rel * pif) +random_symmetric(0.002f);
			//lbm.u.x[n] = 0.05f * cos(10*y_rel * pif) * sin(10*x_rel * pif) +random_symmetric(0.002f);
			//lbm.u.x[n] = random_symmetric(0.01f);
			//lbm.u.y[n] = random_symmetric(0.01f);
			//lbm.u.z[n] = random_symmetric(0.005f);
		}
	}
}
#endif