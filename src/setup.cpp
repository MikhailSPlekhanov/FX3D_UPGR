#include "setup.hpp"
#include "stat_funcs.hpp"

#ifndef BENCHMARK

void main_setup() { // мой тест развития турбулентности
	LBM lbm(512u, 512u, 1u, 0.001f);
	// #############################################################################################################################################################################################
	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
	uint x = 0u, y = 0u, z = 0u;
	for (ulong n = 0ull; n < N; n++) {
		lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if (x >= Nx / 2) {
			lbm.u.y[n] = -0.05f + random_symmetric(0.001f);
		}
		else{
			lbm.u.y[n] = 0.0f + random_symmetric(0.001f);
		}
	}	// #########################################################################################################################################################################################
	lbm.run();
}

// Течение в бесконечном прямоугольном канале

/*void main_setup() {

	LBM lbm(128u, 355u, 226u, 0.005f, 0.00003f, 0.0f, 0.0f);
	float u = 0.15;
	const uint N = lbm.get_N(), Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
	float z_rel, y_rel, x_rel;
	float3 p1 = float3(80.0f, 177.0f, 113.0f);
	float3 p2 = float3(40.0f, 177.0f, 113.0f);
	float3 l = float3(50.0f, 100.0f, 30.0f);
	//std::cout << (float)Ny / 2 << std::endl;
	//std::cout << ((float)y_rel - (float)Ny / 2.0f) / (float)Ny * 2.0f << std::endl;
	uint x, y, z;
	for (uint n = 0u; n < N; n++) {
		lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		z_rel = ((float)z - (float)Nz / 2.0f) / (float)Nz * 2.0f; //должен быть в центре 0, на границах +-1
		y_rel = ((float)y - (float)Ny / 2.0f) / (float)Ny * 2.0f;
		x_rel = ((float)x - (float)Nx / 2.0f) / (float)Nx * 2.0f;
		if (z == 0 || z == Nz - 1 || cylinder(x, y, z, p1, float3(0.0f, 0.0f, 226.0f), 15) || cylinder(x, y, z, p2, float3(0.0f, 355.0f, 0.0f), 15) ) {
			lbm.flags[n] = TYPE_S;
		}
		else {
			//std::cout << "y = " << y << ", " << "ux = " << lbm.u.x[n] << std::endl;
			//lbm.u.x[n] = (1.0f - pow(z_rel, 4)) / 5 * (1.0f + random(0.05));
			//lbm.u.y[n] = (1.0f - pow(z_rel, 4)) / 4 * (random(0.05));
			//lbm.u.z[n] = (1.0f - pow(z_rel, 4)) / 4 * (random(0.05));
			//lbm.rho[n] = 1.0f + 0.1*exp(-sq(x_rel))+0.1*exp(-pow(y_rel-0.2,4));
		}

	};

	key_2 = true;
	Clock clock;
	lbm.run(0u);
	uint time = 0;
	const uint dt_check = 500u;
	while (time < 100000u) {
		if (time == 4000u) {
			std::cout << "trying to rewrite" << std::endl;
			lbm.update_fields();
			lbm.flags.read_from_device();
			lbm.u.read_from_device();
			lbm.rho.read_from_device();
			for (uint n = 0; n < N; n++) {
				lbm.coordinates(n, x, y, z);
				if (z != 0u && z != Nz - 1 && lbm.flags[n] == TYPE_S) {
					lbm.flags[n] = TYPE_F;
					//lbm.u.x[n] = (1.0f - pow(z_rel, 4)) / 5 * (1.0f + random(0.05));;
					//lbm.u.y[n] = (1.0f - pow(z_rel, 4)) / 4 * (random(0.05));
					//lbm.u.z[n] = (1.0f - pow(z_rel, 4)) / 4 * (random(0.05));
					//lbm.rho[n] = 1.0f;
				}
			}
			lbm.flags.write_to_device();
			lbm.rho.write_to_device();
			lbm.u.write_to_device();
		}
		if (time % dt_check == 0) {
			float umax = find_u_max(lbm);
			std::cout << umax << std::endl;
		}
		lbm.run(1u);
		time++;
	}
} /**/

// 3D Taylor-Green vortices

/*void main_setup() {
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(128u, 128u, 128u, 1u, 1u, 1u, 0.01f);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		const float A = 0.25f;
		const uint periodicity = 1u;
		const float a=(float)Nx/(float)periodicity, b=(float)Ny/(float)periodicity, c=(float)Nz/(float)periodicity;
		const float fx = (float)x+0.5f-0.5f*(float)Nx;
		const float fy = (float)y+0.5f-0.5f*(float)Ny;
		const float fz = (float)z+0.5f-0.5f*(float)Nz;
		lbm.u.x[n] =  A*cos(2.0f*pif*fx/a)*sin(2.0f*pif*fy/b)*sin(2.0f*pif*fz/c);
		lbm.u.y[n] = -A*sin(2.0f*pif*fx/a)*cos(2.0f*pif*fy/b)*sin(2.0f*pif*fz/c);
		lbm.u.z[n] =  A*sin(2.0f*pif*fx/a)*sin(2.0f*pif*fy/b)*cos(2.0f*pif*fz/c);
		lbm.rho[n] = 1.0f-sq(A)*3.0f/4.0f*(cos(4.0f*pif*fx/a)+cos(4.0f*pif*fy/b));
	}	// #########################################################################################################################################################################################
	lbm.run();
	//lbm.run(1000u); lbm.u.read_from_device(); println(lbm.u.x[lbm.index(Nx/2u, Ny/2u, Nz/2u)]); wait(); // test for binary identity
} /**/

// 2D Taylor-Green vortices

/*void main_setup() {
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################


	LBM lbm(512u, 512u, 1u, 0.01f);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); 
	uint x=0u, y=0u, z=0u;
	for(ulong n=0ull; n<N; n++) { lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		const float A = 0.05f;
		const uint periodicity = 1u;
		const float a=(float)Nx/(float)periodicity, b=(float)Ny/(float)periodicity;
		const float fx = (float)x+0.5f-0.5f*(float)Nx;
		const float fy = (float)y+0.5f-0.5f*(float)Ny;
		lbm.u.x[n] =  A*cos(2.0f*pif*fx/a)*sin(2.0f*pif*fy/b)+random_symmetric(0.01);
		lbm.u.y[n] = -A*sin(2.0f*pif*fx/a)*cos(2.0f*pif*fy/b) + random_symmetric(0.01);
		lbm.rho[n] = 1.0f-sq(A)*3.0f/4.0f*(cos(4.0f*pif*fx/a)+cos(4.0f*pif*fy/b));
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/

// Poiseuille flow validation

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint R = 63u; // channel radius (default: 63)
	const float umax = 0.1f; // maximum velocity in channel center (must be < 0.57735027f)
	const float tau = 1.0f; // relaxation time (must be > 0.5f), tau = nu*3+0.5
	const float nu = units.nu_from_tau(tau); // nu = (tau-0.5)/3
	const uint H = 2u*(R+1u);
#ifndef D2Q9
	LBM lbm(H, lcm(sq(H), WORKGROUP_SIZE)/sq(H), H, nu, 0.0f, units.f_from_u_Poiseuille_3D(umax, 1.0f, nu, R), 0.0f); // 3D
#else // D2Q9
	LBM lbm(lcm(H, WORKGROUP_SIZE)/H, H, 1u, nu, units.f_from_u_Poiseuille_2D(umax, 1.0f, nu, R), 0.0f, 0.0f); // 2D
#endif // D2Q9
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
#ifndef D2Q9
		if(!cylinder(x, y, z, lbm.center(), float3(0u, Ny, 0u), 0.5f*(float)min(Nx, Nz)-1.0f)) lbm.flags[n] = TYPE_S; // 3D
#else // D2Q9
		if(y==0u||y==Ny-1u) lbm.flags[n] = TYPE_S; // 2D
#endif // D2Q9
	}	// #########################################################################################################################################################################################
	double error_min = max_double;
	while(true) {
		lbm.run(1000);
		lbm.u.read_from_device();
		string s;
		double error_dif=0.0, error_sum=0.0;
#ifndef D2Q9
		for(uint x=0u; x<Nx; x++) {
			for(uint y=Ny/2u; y<Ny/2u+1u; y++) {
				for(uint z=0; z<Nz; z++) {
					const uint n = x+(y+z*Ny)*Nx;
					const double r = (double)sqrt(sq(x+0.5f-0.5f*(float)Nx)+sq(z+0.5f-0.5f*(float)Nz)); // radius from channel center
					if(r<R) {
						const double unum = (double)sqrt(sq(lbm.u.x[n])+sq(lbm.u.y[n])+sq(lbm.u.z[n])); // numerical velocity
						const double uref = umax*(sq(R)-sq(r))/sq(R); // theoretical velocity profile u = G*(R^2-r^2)
						error_dif += sq(unum-uref); // L2 error (Krüger p. 138)
						error_sum += sq(uref);
						s += to_string(r)+" "+to_string(unum)+" "+to_string(uref)+"\n";
					}
				}
			}
		}
#else // D2Q9
		for(uint x=Nx/2u; x<Nx/2u+1u; x++) {
			for(uint y=1u; y<Ny-1u; y++) {
				const uint n = x+(y+0u*Ny)*Nx;
				const double r = (double)(y+0.5f-0.5f*(float)Ny); // radius from channel center
				const double unum = (double)sqrt(sq(lbm.u.x[n])+sq(lbm.u.y[n])); // numerical velocity
				const double uref = umax*(sq(R)-sq(r))/sq(R); // theoretical velocity profile u = G*(R^2-r^2)
				error_dif += sq(unum-uref); // L2 error (Krüger p. 138)
				error_sum += sq(uref);
				s += to_string(r)+" "+to_string(unum)+" "+to_string(uref)+"\n";
			}
		}
#endif // D2Q9
		if(sqrt(error_dif/error_sum)>=error_min) { // stop when error has converged
			print_info("Poiseuille flow error converged after "+to_string(lbm.get_t())+" steps to "+to_string(error_min)); // typical expected L2 errors: 2-5% (Krüger p. 256)
			wait();
			exit(0);
		}
		error_min = fmin(error_min, sqrt(error_dif/error_sum));
		print_info("Poiseuille flow error after t="+to_string(lbm.get_t())+" is "+to_string(error_min)); // typical expected L2 errors: 2-5% (Krüger p. 256)
	}
} /**/

// Stokes drag validation

//void main_setup() { 
//	const uint T = 100; // check error every T steps
//	const float R = 32.0f; // sphere radius
//	const float Re = 0.01f; // Reynolds number
//	const float nu = 1.0f; // kinematic shear viscosity
//	const float rho = 1.0f; // density
//	const uint L = to_uint(8.0f*R); // simulation box size
//	const float u = units.u_from_Re(Re, 2.0f*R, nu); // velocity
//	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
//	LBM lbm(L, L, L, nu); // flow driven by equilibrium boundaries
//	// #############################################################################################################################################################################################
//	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
//		// ########################################################################### define geometry #############################################################################################
//		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E;
//		if(sphere(x, y, z, lbm.center(), R)) {
//			lbm.flags[n] = TYPE_S|TYPE_X; // flag boundary nodes for force summation additionally with TYPE_X
//		} else {
//			lbm.rho[n] = units.rho_Stokes(lbm.position(x, y, z), float3(-u, 0.0f, 0.0f), R, rho, nu);
//			const float3 un = units.u_Stokes(lbm.position(x, y, z), float3(-u, 0.0f, 0.0f), R);
//			lbm.u.x[n] = un.x;
//			lbm.u.y[n] = un.y;
//			lbm.u.z[n] = un.z;
//		}
//	}	// #########################################################################################################################################################################################
//	double E1=1000.0, E2=1000.0;
//	while(true) {
//		lbm.run(T);
//		lbm.calculate_force_on_boundaries();
//		lbm.F.read_from_device();
//		const float3 force = lbm.calculate_force_on_object(TYPE_S|TYPE_X);
//		const double F_theo = units.F_Stokes(rho, u, nu, R);
//		const double F_sim = (double)length(force);
//		const double E0 = fabs(F_sim-F_theo)/F_theo;
//		print_info(to_string(lbm.get_t())+", expected: "+to_string(F_theo, 6u)+", measured: "+to_string(F_sim, 6u)+", error = "+to_string((float)(100.0*E0), 1u)+"%");
//		if(converged(E2, E1, E0, 1E-4)) { // stop when error has sufficiently converged
//			print_info("Error converged after "+to_string(lbm.get_t())+" steps to "+to_string(100.0*E0, 1u)+"%");
//			wait();
//			break;
//		}
//		E2 = E1;
//		E1 = E0;
//	}
//} /*

// cylinder in rectangular duct

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const float Re = 25000.0f;
	const float D = 64.0f;
	const float u = rsqrt(3.0f);
	const float w = D;
	const float h = 3.0f*D;
	const float nu = units.nu_from_Re(Re, D, u);
	const float f = units.f_from_u_rectangular_duct(w, D, 1.0f, nu, u);
	LBM lbm(to_uint(w), 12u*to_uint(D), to_uint(h), nu, 0.0f, f, 0.0f);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		lbm.u.y[n] = 0.1f*u;
		if(cylinder(x, y, z, float3(lbm.center().x, 2.0f*D, lbm.center().z), float3(Nx, 0u, 0u), 0.5f*D)) lbm.flags[n] = TYPE_S;
		if(x==0u||x==Nx-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // x and z non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/

// Taylor-Couette flow

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(96u, 96u, 192u, 1u, 1u, 1u, 0.04f);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(!cylinder(x, y, z, lbm.center(), float3(0u, 0u, Nz), (float)(Nx/2u-1u))) lbm.flags[n] = TYPE_S;
		if( cylinder(x, y, z, lbm.center(), float3(0u, 0u, Nz), (float)(Nx/4u   ))) {
			const float3 relative_position = lbm.relative_position(n);
			lbm.u.x[n] =  relative_position.y;
			lbm.u.y[n] = -relative_position.x;
			lbm.u.z[n] = (1.0f-random(2.0f))*0.001f;
			lbm.flags[n] = TYPE_S;
		}
	}	// #########################################################################################################################################################################################
	lbm.run();
	//lbm.run(4000u); lbm.u.read_from_device(); println(lbm.u.x[lbm.index(Nx/4u, Ny/4u, Nz/2u)]); wait(); // test for binary identity
} /**/

//2D lid-driven cavity

/*void main_setup() {
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 1024u;
	const float Re = 1000.0f;
	const float u = 0.1f;
	LBM lbm(L, L, 1u, units.nu_from_Re(Re, (float)(L-2u), u));
	// #############################################################################################################################################################################################
	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
	uint x = 0u, y = 0u, z = 0u;
	for(ulong n = 0ull; n < N; n++) { lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if (x == 0u || y == 0u || y == Ny - 1u) {
			lbm.flags[n] = TYPE_S; // all non periodic

		}
		if (x == Nx - 1) {
			lbm.flags[n] = TYPE_S;
			lbm.u.y[n] = u;
		}
	}	// #########################################################################################################################################################################################
	
	uint time = 0u;
	float av_fluct;
	LBM_state old = LBM_state(lbm, "u", "", "");
	uint dt_check = 10000u;
	lbm.run(3u); 
	int n, m;
	float usr1, usr2, rhosr1, rhosr2;
	while (time < 1000000u) {
		//if (time % 100000 == 0) { old.copy_u_data(lbm); }
		//if (time % 100000 == 1) { 
		//	av_fluct = average_u_difference(lbm, old); 
		//	std::cout << "time = " << time << std::endl;
		//	std::cout << std::endl << "average fluctuation = " << av_fluct*1000 << "e-3" << std::endl;
		//}
		if (time % dt_check == 0u) {
			usr1 = 0.0f; usr2 = 0.0f; rhosr1 = 0.0f; rhosr2 = 0.0f;
			for (int i = 1; i < Nx - 1u; i++) {
				n = lbm.index(i, Ny / 3u, 1u);
				usr1 += lbm.u.y[n];
				rhosr1 += lbm.rho[n];

				m = lbm.index(i, 2 * Ny / 3u, 1u);
				usr2 += lbm.u.y[m];
				rhosr2 += lbm.rho[m];
			}
			usr1 /= (float)(Nx - 2u); usr2 /= (float)(Nx - 2u); rhosr1 /= (float)(Nx - 2u); rhosr2 /= (float)(Nx - 2u);
			std::cout << "time = " << time << std::endl;
			std::cout << "y = Ny/3, uy = " << usr1 << ", " << "rho = " << rhosr1 << std::endl;
			std::cout << "y = 2Ny/3, uy = " << usr2 << ", " << "rho = " << rhosr2 << std::endl << std::endl;
		}
		lbm.run(1u);
		time++;
	}
} /**/

// particle test

/*void main_setup() {
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 128u;
	const float Re = 1000.0f;
	const float u = 0.1f;
	LBM lbm(L, L, L, units.nu_from_Re(Re, (float)(L-2u), u), 0.0f, 0.0f, -0.00001f, cb(L/4u), 2.0f);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(z==Nz-1) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	for(ulong n=0ull; n<lbm.particles->length(); n++) {
		lbm.particles->x[n] = random_symmetric(0.5f*lbm.size().x/4.0f);
		lbm.particles->y[n] = random_symmetric(0.5f*lbm.size().y/4.0f);
		lbm.particles->z[n] = random_symmetric(0.5f*lbm.size().z/4.0f);
	}
	lbm.run();
} /**/

// delta wing

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 128u;
	const float Re = 10000.0f;
	const float u = 0.1f;
	LBM lbm(L, 4u*L, L, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float3 offset = float3(lbm.center().x, 0.0f, lbm.center().z);
	const float3 p0 = offset+float3(  0*(int)L/64,  5*(int)L/64,  20*(int)L/64);
	const float3 p1 = offset+float3(-20*(int)L/64, 90*(int)L/64, -10*(int)L/64);
	const float3 p2 = offset+float3(+20*(int)L/64, 90*(int)L/64, -10*(int)L/64);
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(triangle(x, y, z, p0, p1, p2)) lbm.flags[n] = TYPE_S;
		else lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/

// Concorde

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 256u;
	const float Re = 1000000.0f;
	const float u = 0.1f;
	LBM lbm(L, L*3u, L/2u, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float size = 1.75f*(float)L;
	const float3 center = float3(lbm.center().x, 0.52f*size, lbm.center().z+0.03f*size);
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(-10.0f))*float3x3(float3(0, 0, 1), radians(90.0f))*float3x3(float3(1, 0, 0), radians(90.0f));
	lbm.voxelize_stl(get_exe_path()+"../stl/Concorde.stl", center, rotation, size); // https://www.thingiverse.com/thing:1176931/files
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	lbm.run();
} /**/

// Boeing 747

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 256u;
	const float Re = 1000000.0f;
	const float u = 0.1f;
	LBM lbm(L, L*2u, L/2u, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float size = 1.0f*(float)L;
	const float3 center = float3(lbm.center().x, 0.55f*size, lbm.center().z);
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(-15.0f));
	lbm.voxelize_stl(get_exe_path()+"../stl/747.stl", center, rotation, size); // https://www.thingiverse.com/thing:2772812/files
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//lbm.graphics.set_camera_free(float3(1.0f*(float)Nx, -0.4f*(float)Ny, 2.0f*(float)Nz), -33.0f, 42.0f, 68.0f);
	//Clock clock;
	//lbm.run(0u);
	//while(lbm.get_t()<100000u) {
	//	lbm.graphics.write_frame_png();
	//	lbm.run(28u);
	//}
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
	lbm.run();
} /**/

// Boeing 757

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 912u;
	const float Re = 100000.0f;
	const float u = 0.125f;
	LBM lbm(L, 2u*L, L/2u, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float size = 1.1f*(float)L;
	const float3 center = float3(lbm.center().x, 32.0f+0.5f*size, lbm.center().z);
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(75.0f));
	lbm.voxelize_stl(get_exe_path()+"../stl/757.stl", center, rotation, size); // https://www.thingiverse.com/thing:5091064/files
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//Clock clock;
	//lbm.run(0u);
	//while(lbm.get_t()<100000u) {
	//	lbm.graphics.set_camera_free(float3(1.0f*(float)Nx, -0.4f*(float)Ny, 2.0f*(float)Nz), -33.0f, 42.0f, 68.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/t/");
	//	lbm.graphics.set_camera_free(float3(0.5f*(float)Nx, -0.35f*(float)Ny, -0.7f*(float)Nz), -35.0f, -35.0f, 100.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/b/");
	//	lbm.graphics.set_camera_free(float3(0.0f*(float)Nx, 0.51f*(float)Ny, 0.75f*(float)Nz), 90.0f, 28.0f, 80.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/f/");
	//	lbm.graphics.set_camera_free(float3(0.6f*(float)Nx, -0.15f*(float)Ny, 0.06f*(float)Nz), 0.0f, 0.0f, 100.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/s/");
	//	lbm.run(28u);
	//}
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
	lbm.run();
} /**/

// Star Wars X-wing

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 256u;
	const float Re = 100000.0f;
	const float u = 0.125f;
	LBM lbm(L, L*2u, L/2u, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float size = 1.0f*(float)L;
	const float3 center = float3(lbm.center().x, 32.0f+0.5f*size, lbm.center().z);
	const float3x3 rotation = float3x3(float3(0, 0, 1), radians(180.0f));
	lbm.voxelize_stl(get_exe_path()+"../stl/X-wing.stl", center, rotation, size); // https://www.thingiverse.com/thing:353276/files
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//Clock clock;
	//lbm.run(0u);
	//while(lbm.get_t()<50000u) {
	//	lbm.graphics.set_camera_free(float3(1.0f*(float)Nx, -0.4f*(float)Ny, 2.0f*(float)Nz), -33.0f, 42.0f, 68.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/t/");
	//	lbm.graphics.set_camera_free(float3(0.5f*(float)Nx, -0.35f*(float)Ny, -0.7f*(float)Nz), -33.0f, -40.0f, 100.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/b/");
	//	lbm.graphics.set_camera_free(float3(0.0f*(float)Nx, 0.51f*(float)Ny, 0.75f*(float)Nz), 90.0f, 28.0f, 80.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/f/");
	//	lbm.graphics.set_camera_free(float3(0.7f*(float)Nx, -0.15f*(float)Ny, 0.06f*(float)Nz), 0.0f, 0.0f, 100.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/s/");
	//	lbm.run(28u);
	//}
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
	lbm.run();
} /**/

// Star Wars TIE fighter

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 256u;
	const float Re = 100000.0f;
	const float u = 0.125f;
	LBM lbm(L, L*2u, L, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float size = 0.65f*(float)L;
	const float3 center = float3(lbm.center().x, 0.6f*size, lbm.center().z);
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(90.0f));
	Mesh* mesh = read_stl(get_exe_path()+"../stl/TIE-fighter.stl", lbm.size(), center, rotation, size); // https://www.thingiverse.com/thing:2919109/files
	lbm.voxelize_mesh_on_device(mesh);
	lbm.flags.read_from_device();
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//Clock clock;
	lbm.run(0u);
	while(lbm.get_t()<50000u) {
	//	lbm.graphics.set_camera_free(float3(1.0f*(float)Nx, -0.4f*(float)Ny, 0.63f*(float)Nz), -33.0f, 33.0f, 80.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/t/");
	//	lbm.graphics.set_camera_free(float3(0.3f*(float)Nx, -1.5f*(float)Ny, -0.45f*(float)Nz), -83.0f, -10.0f, 40.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/b/");
	//	lbm.graphics.set_camera_free(float3(0.0f*(float)Nx, 0.57f*(float)Ny, 0.7f*(float)Nz), 90.0f, 29.5f, 80.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/f/");
	//	lbm.graphics.set_camera_free(float3(2.5f*(float)Nx, 0.0f*(float)Ny, 0.0f*(float)Nz), 0.0f, 0.0f, 50.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/s/");
		lbm.run(28u);
		const float3x3 rotation = float3x3(float3(0.2f, 1.0f, 0.1f), radians(0.4032f)); // create rotation matrix to rotate mesh
		lbm.unvoxelize_mesh_on_device(mesh);
		mesh->rotate(rotation); // rotate mesh
		lbm.voxelize_mesh_on_device(mesh);
	}
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
} /**/

// radial fan (enable MOVING_BOUNDARIES and SUBGRID)

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 872u/4u;
	const float Re = 100000.0f;
	const float u = 0.12f;
	LBM lbm(L, L, L/3u, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float radius = 0.25f*(float)L;
	const float3 center = float3(lbm.center().x, lbm.center().y, 0.36f*radius);
	const uint dt = 10u;
	const float omega=u/radius, domega=omega*dt;
	Mesh* mesh = read_stl(get_exe_path()+"../stl/fan.stl", lbm.size(), center, 2.0f*radius); // https://www.thingiverse.com/thing:6113/files
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//lbm.graphics.set_camera_free(float3(0.353512f*(float)Nx, -0.150326f*(float)Ny, 1.643939f*(float)Nz), -25.0f, 61.0f, 100.0f);
	lbm.run(0u);
	uint k = 0u;
	while(lbm.get_t()<48000u) {
		lbm.voxelize_mesh_on_device(mesh, TYPE_S, center, float3(0.0f), float3(0.0f, 0.0f, omega));
		//if((k++)%4u==0u) lbm.graphics.write_frame();
		lbm.run(dt);
		mesh->rotate(float3x3(float3(0.0f, 0.0f, 1.0f), domega)); // rotate mesh
	}
} /**/

// electric ducted fan (EDF) (enable MOVING_BOUNDARIES and SUBGRID)

//void main_setup() { 
//	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
//	const uint L = 1876u/8u;
//	const float Re = 1000000.0f;
//	const float u = 0.15f;
//	LBM lbm(L, L*3u/2u, L, units.nu_from_Re(Re, (float)L, u));
//	// #############################################################################################################################################################################################
//	const float3 center = lbm.center();
//	const float3x3 rotation = float3x3(float3(0, 0, 1), radians(180.0f));
//	Mesh* rotor = read_stl(get_exe_path()+"../stl/edf_rotor.stl", lbm.size(), center, rotation, -(float)L/102.0f); // https://www.thingiverse.com/thing:3014759/files
//	Mesh* stator = read_stl(get_exe_path()+"../stl/edf_stator.stl", lbm.size(), center, rotation, -(float)L/102.0f); // https://www.thingiverse.com/thing:3014759/files
//	rotor->translate(float3(0.0f, -0.3f*lbm.size().y, 0.0f));
//	stator->translate(float3(0.0f, -0.2f*lbm.size().y, 0.0f));
//	const float radius = 0.5f*rotor->get_max_size();
//	const uint dt = 10u;
//	const float omega=u/radius, domega=omega*dt;
//	lbm.voxelize_mesh_on_device(stator, TYPE_S, center);
//	lbm.voxelize_mesh_on_device(rotor, TYPE_S, center, float3(0.0f), float3(0.0f, omega, 0.0f));
//	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
//		// ########################################################################### define geometry #############################################################################################
//		if(lbm.flags[n]==0u) lbm.u.y[n] = 0.3f*u;
//		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
//	}	// #########################################################################################################################################################################################
//	key_4 = true;
//	lbm.run(0u);
//	uint k = 0u;
//	while(lbm.get_t()<180000u) {
//		//if((k++)%10u==0u) {
//		//	lbm.graphics.set_camera_centered(-70.0f+100.0f*(float)lbm.get_t()/(float)180000u, 2.0f, 60.0f, 1.284025f);
//		//	lbm.graphics.write_frame();
//		//}
//		rotor->rotate(float3x3(float3(0.0f, 1.0f, 0.0f), domega)); // rotate mesh
//		lbm.voxelize_mesh_on_device(rotor, TYPE_S, center, float3(0.0f), float3(0.0f, omega, 0.0f));
//		lbm.run(dt);
//	}
//} /*

// F1 car

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 512u; // 2152u on 8x MI200
	const float kmh = 100.0f;
	const float si_u = kmh/3.6f;
	const float si_x = 2.0f;
	const float si_rho = 1.225f;
	const float si_nu = 1.48E-5f;
	const float Re = units.si_Re(si_x, si_u, si_nu);
	print_info("Re = "+to_string(Re));
	const float u = 0.08f;
	const float size = 1.6f*(float)L;
	units.set_m_kg_s(size*2.0f/5.5f, u, 1.0f, si_x, si_u, si_rho);
	const float nu = units.nu(si_nu);
	print_info("1s = "+to_string(units.t(1.0f)));
	LBM lbm(L, L*2u, L/2u, nu);
	// #############################################################################################################################################################################################
	const float3 center = float3(lbm.center().x, 0.525f*size, 0.116f*size);
	lbm.voxelize_stl(get_exe_path()+"../stl/Ferrari_SF71H_V5.stl", center, size); // https://www.thingiverse.com/thing:2990512/files (unfortunately, this model is not available anymore)
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==Nz-1u) lbm.flags[n] = TYPE_E;
		if(z==0u) lbm.flags[n] = TYPE_S;
		const float3 p = lbm.position(x, y, z);
		const float W = 1.05f*(0.312465f-0.179692f)*(float)Nx;
		const float R = 1.05f*0.5f*(0.361372f-0.255851f)*(float)Ny;
		const float3 FL = float3(0.247597f*(float)Nx, -0.308710f*(float)Ny, -0.260423f*(float)Nz);
		const float3 HL = float3(0.224739f*(float)Nx, 0.210758f*(float)Ny, -0.264461f*(float)Nz);
		const float3 FR = float3(-FL.x, FL.y, FL.z);
		const float3 HR = float3(-HL.x, HL.y, HL.z);
		if((lbm.flags[n]&TYPE_S) && cylinder(x, y, z, lbm.center()+FL, float3(W, 0.0f, 0.0f), R)) {
			const float3 uW = u/R*float3(0.0f, FL.z-p.z, p.y-FL.y);
			lbm.u.y[n] = uW.y;
			lbm.u.z[n] = uW.z;
		}
		if((lbm.flags[n]&TYPE_S) && cylinder(x, y, z, lbm.center()+HL, float3(W, 0.0f, 0.0f), R)) {
			const float3 uW = u/R*float3(0.0f, HL.z-p.z, p.y-HL.y);
			lbm.u.y[n] = uW.y;
			lbm.u.z[n] = uW.z;
		}
		if((lbm.flags[n]&TYPE_S) && cylinder(x, y, z, lbm.center()+FR, float3(W, 0.0f, 0.0f), R)) {
			const float3 uW = u/R*float3(0.0f, FR.z-p.z, p.y-FR.y);
			lbm.u.y[n] = uW.y;
			lbm.u.z[n] = uW.z;
		}
		if((lbm.flags[n]&TYPE_S) && cylinder(x, y, z, lbm.center()+HR, float3(W, 0.0f, 0.0f), R)) {
			const float3 uW = u/R*float3(0.0f, HR.z-p.z, p.y-HR.y);
			lbm.u.y[n] = uW.y;
			lbm.u.z[n] = uW.z;
		}
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//Clock clock;
	//lbm.run(0u);
	//while(lbm.get_t()<=units.t(1.0f)) {
	//	lbm.graphics.set_camera_free(float3(0.779346f*(float)Nx, -0.315650f*(float)Ny, 0.329444f*(float)Nz), -27.0f, 19.0f, 100.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/a/");
	//	lbm.graphics.set_camera_free(float3(0.556877f*(float)Nx, 0.228191f*(float)Ny, 1.159613f*(float)Nz), 19.0f, 53.0f, 100.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/b/");
	//	lbm.graphics.set_camera_free(float3(0.220650f*(float)Nx, -0.589529f*(float)Ny, 0.085407f*(float)Nz), -72.0f, 21.0f, 86.0f);
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/c/");
	//	lbm.run(units.t(0.5f/600.0f)); // run LBM in parallel while CPU is voxelizing the next frame
	//}
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
	lbm.run();
} /**/

// cow

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const float si_u = 1.0f;
	const float si_w = 0.8f;
	const float si_nu = 1.48E-5f;
	const float si_T = 10.0f;
	const float si_rho = 1.225f;
	const uint L = 476u;
	const float u = 0.07f;
	const float w = (float)L/3.0f;
	units.set_m_kg_s(w, u, 1.0f, si_w, si_u, si_rho);
	const float nu = units.nu(si_nu);
	print_info("Re = "+to_string(units.Re(w, u, nu), 2u));
	LBM lbm(L, 2u*L, L, nu);
	// #############################################################################################################################################################################################
	const float size = 1.3f*(float)L;
	const float3 center = float3(lbm.center().x, 0.55f*size, 0.3f*size);
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(180.0f))*float3x3(float3(0, 0, 1), radians(180.0f));
	lbm.voxelize_stl(get_exe_path()+"../stl/Cow.stl", center, rotation, size); // https://www.thingiverse.com/thing:182114/files
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(z==0u) lbm.flags[n] = TYPE_S;
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==Nz-1u) lbm.flags[n] = TYPE_E;
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//lbm.graphics.set_camera_centered(-40.0f, 20.0f, 78.0f, 1.25f);
	//Clock clock;
	//lbm.run(0u);
	//while(lbm.get_t()<=units.t(si_T)) {
	//	lbm.graphics.write_frame_png();
	//	lbm.run(units.t(si_T)/600.0f);
	//}
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
	lbm.run();
} /**/

// Space Shuttle

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 256u; // 1608u
	const float Re = 1000000.0f;
	const float u = 0.1f;
	const uint Lx=(L/8u)*8u, Ly=((L*3u)/8u)*8u, Lz=((L*4u/5u)/8u)*8u;
	print_info("("+to_string(Lx)+"x"+to_string(Ly)+"x"+to_string(Lz)+") = "+to_string((ulong)Lx*(ulong)Ly*(ulong)Lz));
	LBM lbm(Lx, Ly, Lz, 2u, 4u, 1u, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float size = 1.25f*(float)L;
	const float3 center = float3(lbm.center().x, 0.55f*size, lbm.center().z+0.05f*size);
	const float3x3 rotation = float3x3(float3(1, 0, 0), radians(-20.0f))*float3x3(float3(0, 0, 1), radians(270.0f));
	Clock clock;
	lbm.voxelize_stl(get_exe_path()+"../stl/SpaceShuttle.stl", center, rotation, size); // https://www.thingiverse.com/thing:4975964/files
	println(print_time(clock.stop()));
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.y[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//lbm.write_status();
	//Clock clock;
	//lbm.run(0u);
	//while(lbm.get_t()<=60u*60u*30u) {
	//	lbm.graphics.set_camera_free(float3(-1.435962f*(float)Nx, 0.364331f*(float)Ny, 1.344426f*(float)Nz), -205.0f, 36.0f, 74.0f); // top
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/top/");
	//	lbm.graphics.set_camera_free(float3(-1.021207f*(float)Nx, -0.518006f*(float)Ny, 0.0f*(float)Nz), -137.0f, 0.0f, 74.0f); // bottom
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/bottom/");
	//	lbm.run(30u);
	//}
	//lbm.write_status();
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
	lbm.run();
} /**/

// Starship

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 912u;
	const float Re = 1000000.0f;
	const float u = 0.05f;
	LBM lbm(L, L*2u, L*2u, 1u, 2u, 2u, units.nu_from_Re(Re, (float)L, u));
	// #############################################################################################################################################################################################
	const float size = 1.6f*(float)L;
	const float3 center = float3(lbm.center().x, lbm.center().y+0.05f*size, 0.18f*size);
	lbm.voxelize_stl(get_exe_path()+"../stl/Starship.stl", center, size); // https://www.thingiverse.com/thing:4912729/files
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(lbm.flags[n]!=TYPE_S) lbm.u.z[n] = u;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_E; // all non periodic
	}	// #########################################################################################################################################################################################
	key_4 = true;
	//lbm.write_status();
	//Clock clock;
	//lbm.run(0u);
	//while(lbm.get_t()<=20u*60u*90u) {
	//	lbm.graphics.set_camera_free(float3(2.116744f*(float)Nx, -0.775261f*(float)Ny, 1.026577f*(float)Nz), -38.0f, 37.0f, 60.0f); // top
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/top/");
	//	lbm.graphics.set_camera_free(float3(0.718942f*(float)Nx, 0.311263f*(float)Ny, -0.498366f*(float)Nz), 32.0f, -40.0f, 104.0f); // bottom
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/bottom/");
	//	lbm.graphics.set_camera_free(float3(1.748119f*(float)Nx, 0.442782f*(float)Ny, 0.087945f*(float)Nz), 24.0f, 2.0f, 92.0f); // side
	//	lbm.graphics.write_frame_png(get_exe_path()+"export/side/");
	//	lbm.run(90u);
	//}
	//lbm.write_status();
	//write_file(get_exe_path()+"time.txt", print_time(clock.stop()));
	lbm.run();
} /**/

// Test for two domains

/*void main_setup() {
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 128u; // 2152u on 8x MI200
	const float Re = 356.613440755787; const float _sqrt3 = 1 / sqrt(3); const float sqrt3 = sqrt(3);

	const uint vnutr_diam = ceil((float)L / 2 * sqrt3);
	const float vnutr_rad = (float)vnutr_diam / 2;
	const float u_over_nu = Re / vnutr_diam;
	const float uref = 0.15f;
	const float nu = uref / u_over_nu;
	const float h_to_w = 2.02826361766759;
	const float size = (float)L * h_to_w;
	const float vnesh_rad = (float)L / 2;
	const uint height = ceil(L * h_to_w * 0.99f);

	LBM lbm1(vnutr_diam, L, 80, nu, 0.0f, 0.0f, 0.01f);
	LBM lbm2(vnutr_diam, L, 80, nu, 0.0f, 0.01f, 0.001f);

	const ulong N = lbm1.get_N(); const uint Nx = lbm1.get_Nx(), Ny = lbm1.get_Ny(), Nz = lbm1.get_Nz();
	// #############################################################################################################################################################################################
	const float3 center1 = float3(lbm1.center().x, lbm1.center().y, lbm1.center().z);
	const float3 center = lbm1.center();
	float xcent = lbm1.center().x; float ycent = lbm1.center().y; float zcent = lbm1.center().z;

	const uint dt = 10u;
	float xrel, yrel, zrel;
	//lbm.voxelize_stl("C:/Users/Mikhail/source/repos/stl/Small/HexTube.stl", center,  size);
	lbm1.voxelize_stl("C:/Users/Mikhail/source/repos/stl/Small/Combined.stl", center, size);
	lbm2.voxelize_stl("C:/Users/Mikhail/source/repos/stl/Small/Combined.stl", center, size);
	uint m, n; float u_in = 0.0f;

	for (ulong n = 0ull; n < N; n++) {
		uint x = 0u, y = 0u, z = 0u; lbm1.coordinates(n, x, y, z);
		xrel = ((float)x - xcent) / vnutr_rad; yrel = ((float)y - ycent) / vnutr_rad;
		//u_in = (sq(xrel) + sq(yrel) < 1) ? 0.1f * (1 - sq(xrel) - sq(yrel)) : 0.0f;
		// ########################################################################### define geometry #############################################################################################
		if (x == 0u || x == Nx - 1u) { lbm1.flags[n] = TYPE_S; lbm2.flags[n] = TYPE_S; }
		if ((x - xcent) * _sqrt3 + (y - ycent) > (float)L / h_to_w || (x - xcent) * _sqrt3 - (y - ycent) > (float)L / h_to_w || -(x - xcent) * _sqrt3 - (y - ycent) > (float)L / h_to_w || -(x - xcent) * _sqrt3 + (y - ycent) > (float)L / h_to_w) {
			lbm1.flags[n] = TYPE_S;
			lbm2.flags[n] = TYPE_S;
		}
	}
	//################################################################################################################################################################################################

	key_2 = true;
	Clock clock;
	lbm1.run(0u);
	lbm2.run(0u);

	uint time = 0;
	uint dt_check = 200;
	const uint dt_frame = 5000u;
	const float uin = 0.15f;
	uint lay0, lay1, layNz_0, layNz_1, layNz_2;
	float coef;
	float drho = 0.00005f, delta_rho = 0.0f;
	//lbm.graphics.set_camera_free(float3(1.0f * (float)Nx, 0.2f * (float)Ny, 0.5f * (float)Nz), 0.0f, 0.0f, 100.0f);
	//lbm.graphics.set_camera_free(float3(-0.64f * (float)Nx, -0.6f * (float)Ny, 0.6f * (float)Nz), -130.0f, 22.0f, 100.0f);
	uint Ne1 = 0u, Ne2 = 0u;
	

	while (time < 20000) {
		
		if (time % dt_frame == 0u) {
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
			lbm1.graphics.set_camera_free(float3(-0.64f * (float)Nx, -0.6f * (float)Ny, 0.6f * (float)Nz), -130.0f, 22.0f, 100.0f);
			lbm1.graphics.write_frame_png(get_exe_path() + "export/multicase_test/first/");
			lbm2.graphics.set_camera_free(float3(-0.64f * (float)Nx, -0.6f * (float)Ny, 0.6f * (float)Nz), -130.0f, 22.0f, 100.0f);
			lbm2.graphics.write_frame_png(get_exe_path() + "export/multicase_test/second/");
		}
		lbm1.run(1u);
		lbm2.run(1u);
		time++;
	}
	write_file(get_exe_path() + "time.txt", print_time(clock.stop()));
} /**/

// Ламинаризатор

//void main_setup() { 
//	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
//	const uint L = 600u; // 2152u on 8x MI200
//	const uint Lp = 400u;
//	const float u = 0.1f;
//	const float size = 4*(float)Lp;
//
//	LBM lbm(Lp, Lp, L, 0.01, 0.0f, 0.0f, 0.0003f);
//	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
//	float xcent = lbm.center().x; float ycent = lbm.center().y; float zcent = lbm.center().z;
//	const float3 center1 = float3(lbm.center().x, lbm.center().y, lbm.center().z-200);
//	lbm.voxelize_stl("C:/Users/Mikhail/source/repos/ITSBRAS_FX3D/stl/test.stl", center1, size);
//	// #############################################################################################################################################################################################
//	float rsq, phi, xnorm, ynorm;
//	for (ulong n = 0ull; n < N; n++) {
//		uint x = 0u, y = 0u, z = 0u; lbm.coordinates(n, x, y, z);
//		// ########################################################################### define geometry #############################################################################################
//
//		/*if  (abs(y-50.0f) < 10 && abs(z-zcent) < 10) {
//			lbm.flags[n] = TYPE_S;
//		}
//		else{
//			lbm.flags[n] = TYPE_F;
//			lbm.u.y[n] = 0.08f;
//
//		}*/
//		if (lbm.flags[n] != TYPE_S) {
//			//lbm.u.y[n] = 0.05f;
//		}
//		if (y == 0u || y == Ny - 1u || x == 0u || x == Nx - 1u) {
//			//lbm.rho[n] = 1.0f;
//			lbm.flags[n] = TYPE_S;
//		}
//		if (sq(x - xcent) + sq(y - ycent) > sq(Lp * 0.499)){
//			lbm.flags[n] = TYPE_S;
//		}
//		if (z == Nz - 1 && (sq(x - xcent) + sq(y - ycent)) <= sq(Lp * 0.499)) {
//			lbm.flags[n] = TYPE_E;
//			//lbm.u.z[n] = 0.0f;
//			//lbm.rho[n] = 0.99f;
//		}
//		if (z == 0u && (sq(x - xcent) + sq(y - ycent)) <= sq(Lp * 0.499)) {
//			lbm.flags[n] = TYPE_E;
//			rsq = (sq(x - xcent) + sq(y - ycent)) / (sq((float)L / 2));
//			lbm.u.z[n] = 0.1f * (1-rsq);
//			//lbm.rho[n] = 1.01f;
//		}
//		
//	}	// #########################################################################################################################################################################################
//
//	uint time = 0;
//	while (time < 20000) {
//
//		lbm.run(1u);
//		time++;
//	}
//}

// Турбулентная труба

//void main_setup() { 
//	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
//	const uint L = 600u; // 2152u on 8x MI200
//	const uint Lp = 400u;
//	const float u = 0.1f;
//	const float size = 4*(float)Lp;
//	
//	LBM lbm(Lp, Lp, L, 0.01, 0.0f, 0.0f, 0.0003f);
//	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
//	float xcent = lbm.center().x; float ycent = lbm.center().y; float zcent = lbm.center().z;
//	const float3 center1 = float3(lbm.center().x, lbm.center().y, lbm.center().z-200);
//	lbm.voxelize_stl("C:/Users/Mikhail/source/repos/ITSBRAS_FX3D/stl/test.stl", center1, size);
//	// #############################################################################################################################################################################################
//	float rsq, phi, xnorm, ynorm;
//	for (ulong n = 0ull; n < N; n++) {
//		uint x = 0u, y = 0u, z = 0u; lbm.coordinates(n, x, y, z);
//		// ########################################################################### define geometry #############################################################################################
//	
//		/*if  (abs(y-50.0f) < 10 && abs(z-zcent) < 10) {
//			lbm.flags[n] = TYPE_S;
//		}
//		else{
//			lbm.flags[n] = TYPE_F;
//			lbm.u.y[n] = 0.08f;
//	
//		}*/
//		if (lbm.flags[n] != TYPE_S) {
//			//lbm.u.y[n] = 0.05f;
//		}
//		if (y == 0u || y == Ny - 1u || x == 0u || x == Nx - 1u) {
//			//lbm.rho[n] = 1.0f;
//			lbm.flags[n] = TYPE_S;
//		}
//		if (sq(x - xcent) + sq(y - ycent) > sq(Lp * 0.499)){
//			lbm.flags[n] = TYPE_S;
//		}
//		if (z == Nz - 1 && (sq(x - xcent) + sq(y - ycent)) <= sq(Lp * 0.499)) {
//			lbm.flags[n] = TYPE_E;
//			//lbm.u.z[n] = 0.0f;
//			//lbm.rho[n] = 0.99f;
//		}
//		if (z < Nz / 10u && (sq(x - xcent) + sq(y - ycent)) <= sq(Lp * 0.499) && lbm.flags[n] != TYPE_S) {
//			if (z == 0u) {
//				lbm.flags[n] = TYPE_E;
//			}
//			xnorm = (x - xcent) / ((float)Lp / 2);
//			ynorm = (y - ycent) / ((float)Lp / 2);
//			rsq = sqrt(sq(xnorm) + sq(ynorm));
//			phi = atan2(ynorm, xnorm);
//			lbm.u.z[n] = 0.1f * (1-rsq);
//			lbm.u.x[n] = -sin(phi) / 10 * (1-rsq) * rsq;
//			lbm.u.y[n] = cos(phi) / 10 * (1-rsq) * rsq;
//			//lbm.rho[n] = 1.01f;
//		}
//		
//	}	// #########################################################################################################################################################################################
//	
//	uint time = 0;
//	while (time < 20000) {
//	
//		lbm.run(1u);
//		time++;
//	}
//}

#ifdef SURFACE

// hydraulic jump

/*void main_setup() {
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(96u, 352u, 96u, 0.007f, 0.0f, 0.0f, -0.0005f);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		const uint H2 = Nz*3u/5u;
		const uint H1 = Nz*2u/5u;
		const uint P1 = Ny*1u/20u;
		const uint P3 = Ny*3u/20u;
		if(z<H2) lbm.flags[n] = TYPE_F;
		if(y<P3&&z< H1) lbm.flags[n] = TYPE_S;
		if(y<P1&&z>=H2) lbm.flags[n] = TYPE_S;
		if(y==1u&&z>=H1&&z<H2) {
			lbm.flags[n] = TYPE_E;
			lbm.rho[n] = 1.55f;
		}
		if(y==Ny-2u) {
			lbm.flags[n] = TYPE_E;
			lbm.u.y[n] = 0.2f/5.0f;
			lbm.rho[n] = 0.99f;
		}
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
	//lbm.run(1000u); lbm.u.read_from_device(); println(lbm.u.x[lbm.index(Nx/2u, Ny/4u, Nz/4u)]); wait(); // test for binary identity
} /**/

// dam break

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(128u, 256u, 256u, 0.005f, 0.0f, 0.0f, -0.0002f, 0.0001f);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(z<Nz*6u/8u && y<Ny/8u) lbm.flags[n] = TYPE_F;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/

// swirl

/*void main_setup() { 
	const uint D = 128u;
	const float f = 0.000001f;
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(D, D, D, 0.005f, 0.0f, 0.0f, -f);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		const uint H = D*5u/6u;
		const uint R = D/4u-1u;
		if(z<H) {
			lbm.flags[n] = TYPE_F;
			lbm.rho[n] = units.rho_hydrostatic(f, (float)z, (float)H);
		}
		if(cylinder(x, y, z, float3(lbm.center().x, lbm.center().y, 0.0f), float3(0u, 0u, 1u), (float)R)) {
			lbm.u.x[n] =  ((float)y+0.5f-0.5f*(float)Ny)/(float)R*0.25f;
			lbm.u.y[n] = -((float)x+0.5f-0.5f*(float)Nx)/(float)R*0.25f;
		}
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/

// liquid gallium on speaker

/*void main_setup() { 
	const uint L = 128u;
	const float rho = 1.0f;
	const float u = 0.09f; // peak velocity of speaker membrane
	const float nu = 0.01f; // make smaller
	const float sigma = 0.005f; // make larger
	const float f = 0.0005f; // make smaller
	const float frequency = 0.01f; // amplitude = u/(2.0f*pif*frequency);
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(L, L, L*3u/4u, nu, 0.0f, 0.0f, -f, sigma);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(z<Nz/3u && x>0u&&x<Nx-1u&&y>0u&&y<Ny-1u&&z>0u&&z<Nz-1u) {
			lbm.rho[n] = units.rho_hydrostatic(f, (float)z, (float)(Nz/3u));
			lbm.u.x[n] = random_symmetric(1E-9f);
			lbm.u.y[n] = random_symmetric(1E-9f);
			lbm.u.z[n] = random_symmetric(1E-9f);
			lbm.flags[n] = TYPE_F;
		}
		if(z==0u) lbm.u.z[n] = 1E-16f;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	lbm.run(0u);
	while(running) {
		lbm.u.read_from_device();
		const float uz = u*sin(2.0f*pif*frequency*(float)lbm.get_t());
		for(uint z=0u; z<1u; z++) {
			for(uint y=1u; y<Ny-1u; y++) {
				for(uint x=1u; x<Nx-1u; x++) {
					const uint n = x+(y+z*Ny)*Nx;
					lbm.u.z[n] = uz;
				}
			}
		}
		lbm.u.write_to_device();
		lbm.run(1u);
	}
} /**/

// beach

/*void main_setup() { 
	const float f = 0.001f; // make smaller
	const float u = 0.12f; // peak velocity of speaker membrane
	const float frequency = 0.0007f; // amplitude = u/(2.0f*pif*frequency);
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(128u, 640u, 96u, 0.01f, 0.0f, 0.0f, -f);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		const uint H = Nz/2u;
		if(z<H) {
			lbm.flags[n] = TYPE_F;
			lbm.rho[n] = units.rho_hydrostatic(f, (float)z, (float)H);
		}
		if(plane(x, y, z, float3(lbm.center().x, 128.0f, 0.0f), float3(0.0f, -1.0f, 8.0f))) lbm.flags[n] = TYPE_S;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
		if(y==0u && x>0u&&x<Nx-1u&&z>0u&&z<Nz-1u) lbm.flags[n] = TYPE_E;
	}	// #########################################################################################################################################################################################
	lbm.run(0u);
	while(running) {
		lbm.u.read_from_device();
		const float uy = u*sin(2.0f*pif*frequency*(float)lbm.get_t());
		const float uz = 0.5f*u*cos(2.0f*pif*frequency*(float)lbm.get_t());
		for(uint z=1u; z<Nz-1u; z++) {
			for(uint y=0u; y<1u; y++) {
				for(uint x=1u; x<Nx-1u; x++) {
					const uint n = x+(y+z*Ny)*Nx;
					lbm.u.y[n] = uy;
					lbm.u.z[n] = uz;
				}
			}
		}
		lbm.u.write_to_device();
		lbm.run(100u);
	}
} /**/

// river

//void main_setup() { 
//	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
//	LBM lbm(128u, 384u, 96u, 0.02f, 0.0f, -0.00007f, -0.0005f, 0.01f);
//	// #############################################################################################################################################################################################
//	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
//		// ########################################################################### define geometry #############################################################################################
//		const int R = 20, H = 32;
//		if(z==0) lbm.flags[n] = TYPE_S;
//		else if(z<H) {
//			lbm.flags[n] = TYPE_F;
//			lbm.u.y[n] = -0.1f;
//		}
//		if(cylinder(x, y, z, float3(Nx*2u/3u, Ny*2u/3u, Nz/2u)+0.5f, float3(0u, 0u, Nz), (float)R)) lbm.flags[n] = TYPE_S;
//		if(cuboid(x, y, z, float3(Nx/3u, Ny/3u, Nz/2u)+0.5f, float3(2u*R, 2u*R, Nz))) lbm.flags[n] = TYPE_S;
//		if(x==0u||x==Nx-1u) lbm.flags[n] = TYPE_S; // x non periodic
//	}	// #########################################################################################################################################################################################
//	lbm.run();
//} /*

// raindrop setup (#define D3Q19, SRT, VOLUME_FORCE, EQUILIBRIUM_BOUNDARIES, SURFACE)

//void main_setup() { 
//	const int box_diameter = 256; // 936 for Quadro RTX 8000, 1064 for MI200, 2122 on 8x MI200
//	float drop_diameter = box_diameter/5;
//	const int select_drop_size = 12;
//	const float alpha_sim = 20.0f;
//
//	if(drop_diameter==-1.0f) drop_diameter = 0.1f*(float)box_diameter;
//	const float scale = (float)box_diameter/(10.0f*drop_diameter); // 256.0f/400.0f; // 1.0f;
//
//	// rain drop parameters from "Effects of Altitude on Maximum Raindrop Size and Fall Velocity as Limited by Collisional Breakup, Fig. 3" in SI-units
//	float const si_nu = 1.0508E-6f; // kinematic shear viscosity [m^2/s] at 20°C and 35g/l salinity
//	const float si_rho = 1024.8103f; // fluid density [kg/m^3] at 20°C and 35g/l salinity
//	const float si_sigma = 73.81E-3f; // fluid surface tension [kg/s^2] at 20°C and 35g/l salinity
//	const float si_g = 9.81f; // gravitational acceleration [m/s^2]
//	const float alpha = alpha_sim; // impact angle [°], 0 = vertical
//
//	//                            0        1        2        3        4        5        6        7        8        9       10       11       12       13 (13 is for validation)
//	const float si_Ds[] = { 1.0E-3f, 1.5E-3f, 2.0E-3f, 2.5E-3f, 3.0E-3f, 3.5E-3f, 4.0E-3f, 4.5E-3f, 5.0E-3f, 5.5E-3f, 6.0E-3f, 6.5E-3f, 7.0E-3f, 4.1E-3f };
//	const float si_us[] = {   4.50f,   5.80f,   6.80f,   7.55f,   8.10f,   8.45f,   8.80f,   9.05f,   9.20f,   9.30f,   9.40f,   9.45f,   9.55f,   7.21f };
//	const float si_D = si_Ds[select_drop_size]; // drop diameter [m] (1-7mm)
//	const float si_u = si_us[select_drop_size]; // impact velocity [m/s] (4.50-9.55m/s)
//	const float si_H  =  4.0f*si_D*scale; // liquid pool height [m] (4*D is sufficient for deep pool)
//	const float si_Lx = 10.0f*si_D*scale; // simulation box width [m]
//	const float si_Lz =  8.5f*si_D*scale; // simulation box height [m]
//
//	// determine a length, a velocity and the mean density in simulation units
//	const float Lx = (float)box_diameter; // simulation box width
//	const float u = 0.05f; // impact velocity
//	const float rho = 1.0f; // density
//	units.set_m_kg_s(Lx, u, rho, si_Lx, si_u, si_rho); // calculate 3 independent conversion factors (m, kg, s)
//	print_info("D  = "+to_string(si_D, 6u));
//	print_info("Re = "+to_string(units.si_Re(si_D, si_u, si_nu), 6u));
//	print_info("We = "+to_string(units.si_We(si_D, si_u, si_rho, si_sigma), 6u));
//	print_info("Fr = "+to_string(units.si_Fr(si_D, si_u, si_g), 6u));
//	print_info("Ca = "+to_string(units.si_Ca(si_u, si_rho, si_nu, si_sigma), 6u));
//	print_info("Bo = "+to_string(units.si_Bo(si_D, si_rho, si_g, si_sigma), 6u));
//	print_info("10ms = "+to_string(units.t(0.01f))+" LBM steps");
//	const float nu = units.nu(si_nu); // calculate values for remaining parameters in simulation units
//	const float sigma = units.sigma(si_sigma);
//	const float f = units.f(si_rho, si_g); // force per volume
//	const float Lz = units.x(si_Lz);
//	const float H = units.x(si_H);
//	const float R = 0.5f*units.x(si_D); // drop radius
//	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
//	LBM lbm(to_uint(Lx), to_uint(Lx), to_uint(Lz), nu, 0.0f, 0.0f, -f, sigma); // largest box size on Titan Xp with FP32: 384^2, FP16: 464^3
//	// #############################################################################################################################################################################################
//	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); 
//  for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
//		// ########################################################################### define geometry #############################################################################################
//		lbm.rho[n] = rho; // set density everywhere
//		if(sphere(x, y, z, float3(0.5f*(float)Nx, 0.5f*(float)Ny-2.0f*R*tan(alpha*pif/180.0f), H+R+2.5f)+0.5f, R+2.0f)) {
//			const float b = sphere_plic(x, y, z, float3(0.5f*(float)Nx, 0.5f*(float)Ny-2.0f*R*tan(alpha*pif/180.0f)+0.5f, H+R+2.5f), R);
//			if(b!=-1.0f) {
//				lbm.u.y[n] =  sin(alpha*pif/180.0f)*u;//+random_symmetric(0.1f); // break symmetry by initializing with noise
//				lbm.u.z[n] = -cos(alpha*pif/180.0f)*u;//+random_symmetric(0.1f); // break symmetry by initializing with noise
//				if(b==1.0f) {
//					lbm.flags[n] = TYPE_F;
//					lbm.phi[n] = 1.0f;
//				} else {
//					lbm.flags[n] = TYPE_I;
//					lbm.phi[n] = b;
//				}
//			}
//		}
//		if(z==0) lbm.flags[n] = TYPE_S;
//		else if(z==to_uint(H)) {
//			lbm.flags[n] = TYPE_I;
//			lbm.phi[n] = 0.5f; // not strictly necessary, but should be clearer (phi is automatically initialized to 0.5f for TYPE_I if not initialized)
//		} else if((float)z<H) lbm.flags[n] = TYPE_F;
//		else if((x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==Nz-1u)&&(float)z>H+R) { // make drops that hit the simulation box ceiling disappear
//			lbm.rho[n] = 0.5f;
//			lbm.flags[n] = TYPE_E;
//		}
//	}	// #########################################################################################################################################################################################
//	lbm.run(0u);
//	key_1 = false; // turn off boundary
//	key_5 = true; // turn on surface rasterization
//	Clock clock;
//
//	// image
//	//lbm.run(units.t(0.0015f));
//	//print_info("compute time: "+print_time(clock.stop()));
//	//clock.start();
//	//lbm.graphics.set_camera_centered(-30.0f, 20.0f, 100.0f, 1.0f);
//	//lbm.graphics.write_frame_png();
//	//print_info("render time: "+to_string(clock.stop(), 3u));
//
//	// video
//	//while(units.si_t(lbm.get_t())<=0.003f) {
//	//	lbm.graphics.set_camera_centered(-30.0f, 20.0f, 100.0f, 1.0f);
//	//	lbm.graphics.write_frame_png(get_exe_path()+"export/new/");
//	//	lbm.graphics.set_camera_centered(10.0f, 40.0f, 100.0f, 1.0f);
//	//	lbm.graphics.write_frame_png(get_exe_path()+"export/p/");
//	//	lbm.graphics.set_camera_centered(0.0f, 0.0f, 45.0f, 1.0f);
//	//	lbm.graphics.write_frame_png(get_exe_path()+"export/o/");
//	//	lbm.graphics.set_camera_centered(0.0f, 90.0f, 45.0f, 1.0f);
//	//	lbm.graphics.write_frame_png(get_exe_path()+"export/t/");
//	//	lbm.run(units.t(0.004f/600u));
//	//}
//	//print_info("compute + render time: "+to_string(clock.stop(), 3u));
//
//	lbm.run();
//} /*

// bursting bubble setup

/*void main_setup() {
	const float d=32.0f;
	const float sigma=0.0003f;
	const float si_nu = 1E-6f; // kinematic shear viscosity (water) [m^2/s]
	const float si_rho = 1E3f; // density (water) [kg/m^3]
	const float si_sigma = 0.072f; // surface tension (water) [kg/s^2]
	const float si_d = 4E-3f; // bubble diameter [m]
	const float si_g = 9.81f; // gravitational acceleration [m/s^2]
	const float si_f = units.si_f_from_si_g(si_g, si_rho);
	const float si_rho_particles = si_rho;
	const float rho = 1.0f;
	const float m = si_d/d; // length si_x = x*[m]
	const float kg = si_rho/rho*cb(m); // density si_rho = rho*[kg/m^3]
	const float s = sqrt(sigma/si_sigma*kg); // velocity si_sigma = sigma*[kg/s^2]
	units.set_m_kg_s(m, kg, s);
	const float f = units.f(si_f);
	const float nu = units.nu(si_nu);
	const uint Lx = to_uint(4.0f*d);
	const uint Ly = to_uint(4.0f*d);
	const uint Lz = d==112.0f ? to_uint(2.5f*d) : to_uint(3.0f*d);
	const uint H = to_uint(2.0f*d);
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(Lx, Ly, Lz, nu, 0.0f, 0.0f, -f, sigma);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz();
	uint x = 0u, y = 0u, z = 0u;
	for(ulong n=0ull; n<N; n++) { lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(z < H) lbm.flags[n] = TYPE_F;
		const float r = 0.5f*d;
		if(sphere(x, y, z, float3(lbm.center().x, lbm.center().y, (float)H-0.5f*d), r+1.0f)) { // bubble
			const float b = clamp(sphere_plic(x, y, z, float3(lbm.center().x, lbm.center().y, (float)H-0.5f*d), r), 0.0f, 1.0f);
			if(b == 1.0f) {
				lbm.flags[n] = TYPE_G;
				lbm.phi[n] = 0.0f;
			} else {
				lbm.flags[n] = TYPE_I;
				lbm.phi[n] = (1.0f-b);
			}
		}
		if(z == 0) lbm.flags[n] = TYPE_S;
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/

// cube with changing gravity

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(96u, 96u, 96u, 0.02f, 0.0f, 0.0f, -0.001f, 0.001f);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(x<Nx*2u/3u&&y<Ny*2u/3u) lbm.flags[n] = TYPE_F;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S;
	}	// #########################################################################################################################################################################################
	while(running) {
		lbm.set_f(0.0f, 0.0f, -0.001f);
		lbm.run(2500u);
		lbm.set_f(0.0f, +0.001f, 0.0f);
		lbm.run(2500u);
		lbm.set_f(0.0f, 0.0f, +0.001f);
		lbm.run(2500u);
		lbm.set_f(0.0f, -0.001f, 0.0f);
		lbm.run(2000u);
		lbm.set_f(0.0f, 0.0f, 0.0f);
		lbm.run(3000u);
	}
} /**/

// periodic faucet mass conservation test

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(96u, 192u, 128u, 0.02f, 0.0f, 0.0f, -0.001f);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); 
	for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(y>Ny*5/6) lbm.flags[n] = TYPE_F;
		const uint D = max(Nx, Nz);
		const uint r = D/6;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u) lbm.flags[n] = TYPE_S; // x and y non periodic
		if((z==0u||z==Nz-1u) && sq(x-Nx/2)+sq(y-Nx/2)>sq(r)) lbm.flags[n] = TYPE_S; // z non periodic
		if(y<=Nx/2+2*r && torus_x(x, y, z, float3(Nx/2, Nx/2+r, Nz)+0.5f, (float)r, (float)r)) lbm.flags[n] = TYPE_S;
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/

// two colliding droplets in force field

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(256u, 256u, 128u, 0.014f, 0.0f, 0.0f, 0.0f, 0.0001f);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); for(ulong n=0ull; n<N; n++) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if(sphere(x, y, z, lbm.center()-float3(0u, 10u, 0u), 32.0f)) {
			lbm.flags[n] = TYPE_F;
			lbm.u.y[n] = 0.025f;
		}
		if(sphere(x, y, z, lbm.center()+float3(30u, 40u, 0u), 12.0f)) {
			lbm.flags[n] = TYPE_F;
			lbm.u.y[n] = -0.2f;
		}
		lbm.F.x[n] = -0.001f*lbm.relative_position(n).x;
		lbm.F.y[n] = -0.001f*lbm.relative_position(n).y;
		lbm.F.z[n] = -0.0005f*lbm.relative_position(n).z;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
	}	// #########################################################################################################################################################################################
	lbm.run();
} /**/



#endif // SURFACE
#ifdef TEMPERATURE

// Rayleigh-Benard convection

void main_setup() {
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(256u, 256u, 4u, 0.005f, 0.0f, 0.0f, 0.0f, 0.0f, 0.01f, 0.1f);
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); uint x=0u, y=0u, z=0u;
	float z_rel, y_rel, x_rel;
	for(ulong n=0ull; n<N; n++) { lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		//lbm.u.x[n] = random_symmetric(0.02f);
		//lbm.u.y[n] = random_symmetric(0.02f);
		//lbm.u.z[n] = random_symmetric(0.005f);
		z_rel = ((float)z - (float)Nz / 2.0f) / (float)Nz * 2.0f;
		y_rel = ((float)y - (float)Ny / 2.0f) / (float)Ny * 2.0f;
		x_rel = ((float)x - (float)Nx / 2.0f) / (float)Ny * 2.0f;
		if(y == 1u) {
			lbm.T[n] = 1.5f;
			lbm.flags[n] = TYPE_T;
		} else if(y == Ny - 2u) {
			lbm.T[n] = 0.5f;
			lbm.flags[n] = TYPE_T;
		}
		else if (y == 0u || y == Ny - 1u) {
			lbm.flags[n] = TYPE_S; // z non periodic
			//lbm.u.x[n] = random(0.02f) - 0.01f;	lbm.u.y[n] = random(0.02f) - 0.01f;
		}
		else if (x == 0u || x == Nx - 1u) lbm.flags[n] = TYPE_S;
		else {
			lbm.T[n] = 1.0f - y_rel * 0.5f;
			//lbm.rho[n] = 1.001f - y_rel * 0.0001f;
			//lbm.u.y[n] = 0.05f * sin((y_rel-1.0f)/2 * 3*pif ) * cos((x_rel-1.0f)/2 * 3*pif ) ;
			//lbm.u.x[n] = -0.05f * cos((y_rel-1.0f)/2 * 3*pif ) * sin((x_rel-1.0f)/2 * 3*pif  ) ;
			lbm.u.y[n] = 0.05f * sin(y_rel * pif) * cos(x_rel * pif);
			lbm.u.x[n] = -0.05f * cos(y_rel * pif) * sin(x_rel * pif);
		}
		//else lbm.rho[n] = 1 - (z_rel - 1.0f);
		//if(z==0u||z==Nz-1u||x==0u||x==Nx-1u||y==0u||y==Ny-1u) lbm.flags[n] = TYPE_S;
	}	// #########################################################################################################################################################################################
	uint time = 0u;
	Clock clock1, clock2;
	double for_time = 0.0;
	int times = 0;
	float av_dif, nusselt1, nusselt2;
	//LBM_state old = LBM_state(lbm, "u", "", "");
	//int MAX_THREADS = omp_get_max_threads(); omp_set_num_threads(MAX_THREADS);
	std::vector<std::vector<float>> utemp = { { 1.0f } };
	std::vector<std::vector<float>> dtemp = { { 1.0f } };
	//set_complex_zplane_temp(lbm, utemp, (int)(Nz - 2u));
	//set_complex_zplane_temp(lbm, dtemp, (int)(1u));
	float fy = 0.0f, delta_T = 0.0f;
	//lbm.run(0u);
	while (time < 1000001) {
		if (lbm.get_fy() > -0.001f) {
			fy -= 0.00000001f;
			lbm.set_fy(fy);
		}
		/*if (delta_T < 0.3f) {
			delta_T += 0.0001f;
			dtemp[0][0] = 1 + delta_T;
			utemp[0][0] = 1 - delta_T;
			set_complex_plane_temp(lbm, dtemp, 1, 'y');
			set_complex_plane_temp(lbm, utemp, Ny-2, 'y');
		}*/
		/*if (time % 100 == 0) {
			float umax = find_u_max(lbm);
			std::cout << umax << std::endl;
			if (umax > 0.22) {
				key_P = false;
			}
		}*/
		/*if (time % 1000 == 0) {
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Nz; j++) {
					int n = lbm.index(i, 1, j);
					std::cout << lbm.T[n] << " ";
				}
				std::cout << std::endl;
			}
		}*/
		//if (time % 2 == 0) { old.copy_u_data(lbm); }
		//if (time % 2 == 1) {
		//	av_dif = average_u_difference(lbm, old);
		//	//std::cout << "average difference = " << av_dif * 1000 << " e-3" << std::endl;
		//}
		//old.copy_u_data(lbm);
		//
		//av_dif = average_u_difference(lbm, old);
		// times++;
		//if (time % 100u == 0) {
		//	clock2.start();
		//	nusselt1 = Rayleigh_Benard_Nusselt(lbm);
		//	std::cout << "nusselt number with omp = " << nusselt1 << ", elapsed time: " << clock2.stop()*1000 << " ms" << std::endl;

		//	clock2.start();
		//	nusselt2 = Rayleigh_Benard_Nusselt_noomp(lbm);
		//	std::cout << "nusselt number with no omp = " << nusselt2 << ", elapsed time: " << clock2.stop() * 1000 << " ms" << std::endl;
		//	//break;
		//}
		lbm.run(1u);
		//key_P = false;
		//clock2.start();
		//set_complex_zplane_temp_noomp(lbm, temp, (int)(Nz - 2u));
		//std::cout << "No omp time: " << clock2.stop() * 1000 << " ms " << std::endl;
		//lbm.run(10u);
		//clock2.start();
		//set_complex_zplane_temp(lbm, temp, (int)(Nz - 2u));
		//std::cout << "omp time: " << clock2.stop() * 1000 << " ms " << std::endl;
		time++;
	}
	//std::cout << "copy time is equal to " << for_time * 1000 / (float)times << " e-3 seconds" << std::endl;
	write_file(get_exe_path() + "time.txt", print_time(clock1.stop()));
} /**/

// Пуазейль с традиентом температур

//void main_setup() { 
//	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
//	LBM lbm(256u, 64u, 256u, 0.02f, 0.0f, 0.0f, -0.00001f, 0.0f, 0.1f, 0.1f);
//	// #############################################################################################################################################################################################
//	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz();
//	uint x=0u, y=0u, z=0u;
//	float x_rel, y_rel;
//	for(ulong n=0ull; n<N; n++) { lbm.coordinates(n, x, y, z);
//		// ########################################################################### define geometry #############################################################################################
//		//y_rel = ((float)y - (float)Ny / 2.0f) / (float)Ny * 2.0f;
//		//x_rel = ((float)x - (float)Nx / 2.0f) / (float)Nx * 2.0f;
//		lbm.u.x[n] = random_symmetric(0.005f);
//		lbm.u.y[n] = random_symmetric(0.005f);
//		lbm.u.z[n] = random_symmetric(0.005f);
//		if (y == 0u || y == Ny - 1u) {
//			lbm.flags[n] = TYPE_S;
//		}
//		else if(y == 1u) {
//			lbm.T[n] = 1.2f;
//			lbm.flags[n] = TYPE_T;
//		} else if(y == Nz - 2u) {
//			lbm.T[n] = 0.8f;
//			lbm.flags[n] = TYPE_T;
//		}
//		else {
//			lbm.flags[n] = TYPE_F;
//			//lbm.rho[n] = 1.0f + 0.01 * exp(-sq(x_rel)) + 0.01f * exp(-sq(y_rel));
//		}
//		//if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
//
//	}	// #########################################################################################################################################################################################
//	lbm.run();
//} /*

// TEMPERATURE test

/*void main_setup() { 
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	LBM lbm(32u, 196u, 60u, 0.05f, 0.0f, 0.0f, -0.001f, 0.0f, 0.1f, 0.1f);
	float y_rel;
	// #############################################################################################################################################################################################
	const ulong N=lbm.get_N(); const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz();
	uint x = 0u, y = 0u, z = 0u;
	for(ulong n=0ull; n<N; n++) { lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		
		if(y == 1) {
			lbm.T[n] = 1.5f;
			lbm.flags[n] = TYPE_T;
		} else if(y == Ny - 2) {
			lbm.T[n] = 0.5f;
			lbm.flags[n] = TYPE_T;
		}
		else if( x == 0u || x == Nx - 1u || y == 0u || y == Ny - 1u || z == 0u || z == Nz - 1u ) lbm.flags[n] = TYPE_S; // all non periodic
		else if (z < Nz * 0.8f) {
			lbm.flags[n] = TYPE_F;

		}
	}	// #########################################################################################################################################################################################
	lbm.run();
	//lbm.run(1000u); lbm.u.read_from_device(); println(lbm.u.x[lbm.index(Nx/2u, Ny/2u, Nz/2u)]); wait(); // test for binary identity
} /**/

#endif // TEMPERATURE

// BENCHMARK
#else 
#include "info.hpp"
//void main_setup() { // benchmark
//	uint mlups = 0u;
//	{ // ######################################################## define simulation box size, viscosity and volume force ###########################################################################
//		//LBM lbm( 32u,  32u,  32u, 1.0f);
//		//LBM lbm( 48u,  48u,  48u, 1.0f);
//		//LBM lbm( 64u,  64u,  64u, 1.0f);
//		//LBM lbm( 96u,  96u,  96u, 1.0f);
//		//LBM lbm(128u, 128u, 128u, 1.0f);
//		//LBM lbm(192u, 192u, 192u, 1.0f);
//		LBM lbm(256u, 256u, 256u, 1.0f);
//		//LBM lbm(384u, 384u, 384u, 1.0f);
//		//LBM lbm(464u, 464u, 464u, 1.0f);
//		//LBM lbm(480u, 480u, 480u, 1.0f);
//		//LBM lbm(512u, 512u, 512u, 1.0f);
//
//		//const uint memory = 1488u; // in MB
//		//const uint L = ((uint)cbrt(fmin((float)memory*1048576.0f/(19.0f*(float)sizeof(fpxx)+17.0f), (float)max_uint))/2u)*2u;
//		//LBM lbm(1u*L, 1u*L, 1u*L, 1u, 1u, 1u, 1.0f); // 1 GPU
//		//LBM lbm(2u*L, 1u*L, 1u*L, 2u, 1u, 1u, 1.0f); // 2 GPUs
//		//LBM lbm(2u*L, 2u*L, 1u*L, 2u, 2u, 1u, 1.0f); // 4 GPUs
//		//LBM lbm(2u*L, 2u*L, 2u*L, 2u, 2u, 2u, 1.0f); // 8 GPUs
//		// #########################################################################################################################################################################################
//		for(uint i=0u; i<1000u; i++) {
//			lbm.run(10u);
//			mlups = max(mlups, to_uint((double)lbm.get_N()*1E-6/info.dt_smooth));
//		}
//	} // make lbm object go out of scope to free its memory
//	print_info("Peak MLUPs/s = "+to_string(mlups));
//#if defined(_WIN32)
//	wait();
//#endif // Windows
//} /**/
#endif // BENCHMARK