#include "mercedes.hpp"

void mercedes() { // Μεπρεδερ
	// ######################################################### define simulation box size, viscosity and volume force ############################################################################
	const uint L = 256u; // 2152u on 8x MI200
	const float Re = 100000.0f;
	const float u = 0.1f;
	const float size = (float)L;
	LBM lbm(L, L * 3u, L / 2u, units.nu_from_Re(Re, (float)(L/2u), u));

	// #############################################################################################################################################################################################
	const float3 center = float3(lbm.center().x, 0.6f * size, 0.116f * size);
	uint xcent = lbm.center().x; uint ycent = lbm.center().y; uint zcent = lbm.center().z;
	//lbm.voxelize_stl(get_exe_path() + "../../stl/Mercedes-Benz_clk-gtr.stl", center + (float3)(0, 0, 4), size);
	lbm.voxelize_stl("C:/Users/Mikhail/source/repos/stl/Mercedes-Benz_clk-gtr.stl", center + (float3)(0, 0, (float)L/130), size);

	const ulong N = lbm.get_N(); const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
	for (ulong n = 0ull; n < N; n++) { uint x = 0u, y = 0u, z = 0u; lbm.coordinates(n, x, y, z);
		// ########################################################################### define geometry #############################################################################################
		if (lbm.flags[n] != TYPE_S) lbm.u.y[n] = u;
		if (x == 0u || x == Nx - 1u || y == 0u || y == Ny - 1u || z == Nz - 1u || z == 0u) {
			lbm.flags[n] = TYPE_E; // all non periodic
		}
	}	// #########################################################################################################################################################################################

	//key_4 = true;
	Clock clock;
	lbm.run(0u);
	uint time = 0;
	uint dt = 20;
	string export_path = "export/mercedes/MultiGPU/";
	//lbm.graphics.set_camera_free(float3(1.5f * (float)Nx, -0.05f * (float)Ny, 0.0f * (float)Nz), 0.0f, 0.0f, 86.0f);
	//lbm.graphics.set_camera_free(float3(0.02f * (float)Nx, -0.23f * (float)Ny, 1.57f * (float)Nz), 90.0f, 82.0f, 90.0f);
	while (time < 20001) {//lbm.get_t()<=units.t(1.0f)) {
#ifdef GRAPHICS
		//lbm.graphics.set_camera_free(float3(0.6f*(float)Nx, -0.499f*(float)Ny, 0.3f*(float)Nz), -50.0f, 15.0f, 100.0f);
		//if (time % dt == 0) lbm.graphics.write_frame_png(get_exe_path() + export_path + "/a/");
		lbm.graphics.set_camera_free(float3(0.02f * (float)Nx, -0.23f * (float)Ny, 1.57f * (float)Nz), 90.0f, 82.0f, 90.0f);
		if (time % dt == 0) lbm.graphics.write_frame_png(get_exe_path() + export_path + "/b/");
		lbm.graphics.set_camera_free(float3(1.5f * (float)Nx, -0.05f * (float)Ny, 0.0f * (float)Nz), 0.0f, 0.0f, 86.0f);
		if (time % dt == 0) lbm.graphics.write_frame_png(get_exe_path() + export_path + "/c/");
#endif
		/////////////////////////////////////lbm.run(units.t(0.5f/600.0f)); // run LBM in parallel while CPU is voxelizing the next frame
		lbm.run(1u);
		time++;
	}
	write_file(get_exe_path() + export_path + "time.txt", print_time(clock.stop()));
	//lbm.run();
}