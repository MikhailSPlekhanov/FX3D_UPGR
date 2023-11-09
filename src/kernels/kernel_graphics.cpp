#include "../kernel.hpp"

string kernel_graphics() {
	return R( 
		// ################################################## graphics code ##################################################
	) + "#ifdef GRAPHICS" + R(

	) + "#ifndef FORCE_FIELD" + R( // render flags as grid
	) + R(kernel void graphics_flags(const global uchar * flags, const global float* camera, global int* bitmap, global int* zbuffer) {
	) + "#else" + R( // FORCE_FIELD
	) + R(kernel void graphics_flags(const global uchar * flags, const global float* camera, global int* bitmap, global int* zbuffer, const global float* F) {
		) + "#endif" + R( // FORCE_FIELD
			const uint n = get_global_id(0);
		if (n >= (uint)def_N || is_halo(n)) return; // don't execute graphics_flags() on halo
		const uchar flagsn = flags[n]; // cache flags
		const uchar flagsn_bo = flagsn & TYPE_BO; // extract boundary flags
		if (flagsn == 0u || flagsn == TYPE_G) return; // don't draw regular fluid nodes
		//if(flagsn&TYPE_SU) return; // don't draw surface
		float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
		for (uint i = 0u; i < 15u; i++) camera_cache[i] = camera[i];
		uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
		calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
		const uint3 xyz = coordinates(n);
		const float3 p = position(xyz);
		const int c =  // coloring scheme
			flagsn_bo == TYPE_S ? COLOR_S : // solid boundary
			((flagsn & TYPE_T) && flagsn_bo == TYPE_E) ? color_mix(COLOR_T, COLOR_E, 0.5f) : // both temperature boundary and equilibrium boundary
			((flagsn & TYPE_T) && flagsn_bo == TYPE_MS) ? color_mix(COLOR_T, COLOR_M, 0.5f) : // both temperature boundary and moving boundary
			flagsn & TYPE_T ? COLOR_T : // temperature boundary
			flagsn_bo == TYPE_E ? COLOR_E : // equilibrium boundary
			flagsn_bo == TYPE_MS ? COLOR_M : // moving boundary
			flagsn & TYPE_F ? COLOR_F : // fluid
			flagsn & TYPE_I ? COLOR_I : // interface
			flagsn & TYPE_X ? COLOR_X : // reserved type X
			flagsn & TYPE_Y ? COLOR_Y : // reserved type Y
			COLOR_0; // regular or gas node
		//draw_point(p, c, camera_cache, bitmap, zbuffer); // draw one pixel for every boundary node
		uint t;
		t = xp + y0 + z0; const bool not_xp = xyz.x < def_Nx - 1u && flagsn == flags[t] && !is_halo(t); // +00
		t = xm + y0 + z0; const bool not_xm = xyz.x > 0u && flagsn == flags[t] && !is_halo(t); // -00
		t = x0 + yp + z0; const bool not_yp = xyz.y < def_Ny - 1u && flagsn == flags[t] && !is_halo(t); // 0+0
		t = x0 + ym + z0; const bool not_ym = xyz.y > 0u && flagsn == flags[t] && !is_halo(t); // 0-0
		t = x0 + y0 + zp; const bool not_zp = xyz.z < def_Nz - 1u && flagsn == flags[t] && !is_halo(t); // 00+
		t = x0 + y0 + zm; const bool not_zm = xyz.z > 0u && flagsn == flags[t] && !is_halo(t); // 00-
		const float3 p0 = (float3)(p.x - 0.5f, p.y - 0.5f, p.z - 0.5f); // ---
		const float3 p1 = (float3)(p.x + 0.5f, p.y + 0.5f, p.z + 0.5f); // +++
		const float3 p2 = (float3)(p.x - 0.5f, p.y - 0.5f, p.z + 0.5f); // --+
		const float3 p3 = (float3)(p.x + 0.5f, p.y + 0.5f, p.z - 0.5f); // ++-
		const float3 p4 = (float3)(p.x - 0.5f, p.y + 0.5f, p.z - 0.5f); // -+-
		const float3 p5 = (float3)(p.x + 0.5f, p.y - 0.5f, p.z + 0.5f); // +-+
		const float3 p6 = (float3)(p.x + 0.5f, p.y - 0.5f, p.z - 0.5f); // +--
		const float3 p7 = (float3)(p.x - 0.5f, p.y + 0.5f, p.z + 0.5f); // -++
		if (!(not_xm || not_ym)) draw_line(p0, p2, c, camera_cache, bitmap, zbuffer); // to draw the entire surface, replace || by &&
		if (!(not_xm || not_zm)) draw_line(p0, p4, c, camera_cache, bitmap, zbuffer);
		if (!(not_ym || not_zm)) draw_line(p0, p6, c, camera_cache, bitmap, zbuffer);
		if (!(not_xp || not_yp)) draw_line(p1, p3, c, camera_cache, bitmap, zbuffer);
		if (!(not_xp || not_zp)) draw_line(p1, p5, c, camera_cache, bitmap, zbuffer);
		if (!(not_yp || not_zp)) draw_line(p1, p7, c, camera_cache, bitmap, zbuffer);
		if (!(not_ym || not_zp)) draw_line(p2, p5, c, camera_cache, bitmap, zbuffer);
		if (!(not_xm || not_zp)) draw_line(p2, p7, c, camera_cache, bitmap, zbuffer);
		if (!(not_yp || not_zm)) draw_line(p3, p4, c, camera_cache, bitmap, zbuffer);
		if (!(not_xp || not_zm)) draw_line(p3, p6, c, camera_cache, bitmap, zbuffer);
		if (!(not_xm || not_yp)) draw_line(p4, p7, c, camera_cache, bitmap, zbuffer);
		if (!(not_xp || not_ym)) draw_line(p5, p6, c, camera_cache, bitmap, zbuffer);
		) + "#ifdef FORCE_FIELD" + R(
			if (flagsn_bo == TYPE_S) {
				const float3 Fn = def_scale_F * (float3)(F[n], F[def_N + (ulong)n], F[2ul * def_N + (ulong)n]);
				const float Fnl = length(Fn);
				if (Fnl > 0.0f) {
					const int c = iron_color(255.0f * Fnl); // color boundaries depending on the force on them
					draw_line(p, p + 5.0f * Fn, c, camera_cache, bitmap, zbuffer); // draw colored force vectors
				}
			}
		) + "#endif" + R( // FORCE_FIELD
	}/**/

	/*) + "#ifndef FORCE_FIELD" + R( // render solid boundaries with marching-cubes
	)+R(kernel void graphics_flags(const global uchar* flags, const global float* camera, global int* bitmap, global int* zbuffer) {
	)+"#else"+R( // FORCE_FIELD
	) + R(kernel void graphics_flags(const global uchar * flags, const global float* camera, global int* bitmap, global int* zbuffer, const global float* F) {
		) + "#endif" + R( // FORCE_FIELD
			const uint n = get_global_id(0);
		if (n >= (uint)def_N || is_halo(n)) return; // don't execute graphics_flags() on halo
		const uint3 xyz = coordinates(n); const float xcent = (float)def_Nx * 0.5f - 0.5f; const float ycent = (float)def_Ny * 0.5f - 0.5f; const float h_to_w = 2.02826361766759; const uint vnutr_diam = ceil((float)def_Ny / 2 * sqrt(3.0f));
		if (xyz.x >= def_Nx - 1u || xyz.y >= def_Ny - 1u || xyz.z >= def_Nz - 1u) return;
		//if (xyz.z < def_Nz / 5u) return;
		//if(xyz.x==0u||xyz.y==0u||xyz.z==0u||xyz.x>=def_Nx-2u||xyz.y>=def_Ny-2u||xyz.z>=def_Nz-2u) return;
		//if (((float)xyz.x - xcent) / sqrt(3.0f) - ((float)xyz.y - ycent) > (float)def_Ny / h_to_w) return;
		if ((xcent - (float)xyz.x) / sqrt(3.0f) - ((float)xyz.y - ycent) > (float)def_Ny / h_to_w || xyz.x < 3u || xyz.z > def_Nz - 3u) return; //для шестиугольника
		//if (xyz.y < def_Ny / 2u) return;
		uint j[8];
		const uint x0 = xyz.x; // cube stencil
		const uint xp = xyz.x + 1u;
		const uint y0 = xyz.y * def_Nx;
		const uint yp = (xyz.y + 1u) * def_Nx;
		const uint z0 = xyz.z * def_Ny * def_Nx;
		const uint zp = (xyz.z + 1u) * def_Ny * def_Nx;
		j[0] = n; // 000
		j[1] = xp + y0 + z0; // +00
		j[2] = xp + y0 + zp; // +0+
		j[3] = x0 + y0 + zp; // 00+
		j[4] = x0 + yp + z0; // 0+0
		j[5] = xp + yp + z0; // ++0
		j[6] = xp + yp + zp; // +++
		j[7] = x0 + yp + zp; // 0++
		float v[8];
		for (uint i = 0u; i < 8u; i++) v[i] = (float)((flags[j[i]] & TYPE_BO) == TYPE_S);
		float3 triangles[15]; // maximum of 5 triangles with 3 vertices each
		const uint tn = marching_cubes(v, 0.5f, triangles); // run marching cubes algorithm
		if (tn == 0u) return;
		float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
		for (uint i = 0u; i < 15u; i++) camera_cache[i] = camera[i];
		const float3 offset = (float3)((float)xyz.x + 0.5f - 0.5f * (float)def_Nx, (float)xyz.y + 0.5f - 0.5f * (float)def_Ny, (float)xyz.z + 0.5f - 0.5f * (float)def_Nz);
		for (uint i = 0u; i < tn; i++) {
			const float3 p0 = triangles[3u * i] + offset;
			const float3 p1 = triangles[3u * i + 1u] + offset;
			const float3 p2 = triangles[3u * i + 2u] + offset;
			const float3 p = (p0 + p1 + p2) / 3.0f, normal = cross(p1 - p0, p2 - p0);
			const int c = lighting(191 << 16 | 191 << 8 | 191, p, normal, camera_cache);
			draw_triangle(p0, p1, p2, c, camera_cache, bitmap, zbuffer);
		}
		) + "#ifdef FORCE_FIELD" + R(
			const uchar flagsn_bo = flags[n] & TYPE_BO;
		const float3 p = position(xyz);
		if (flagsn_bo == TYPE_S) {
			const float3 Fn = def_scale_F * (float3)(F[n], F[def_N + (ulong)n], F[2ul * def_N + (ulong)n]);
			const float Fnl = length(Fn);
			if (Fnl > 0.0f) {
				const int c = iron_color(255.0f * Fnl); // color boundaries depending on the force on them
				draw_line(p, p + 5.0f * Fn, c, camera_cache, bitmap, zbuffer); // draw colored force vectors
			}
		}
		) + "#endif" + R( // FORCE_FIELD
	}/**/

	) + R(kernel void graphics_field(const global uchar * flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer) {
		const uint n = get_global_id(0);
		if (n >= (uint)def_N || is_halo(n)) return; // don't execute graphics_field() on halo
		) + "#ifndef MOVING_BOUNDARIES" + R(
			if (flags[n] & (TYPE_S | TYPE_E | TYPE_I | TYPE_G)) return;
		) + "#else" + R( // EQUILIBRIUM_BOUNDARIES
			if (flags[n] & (TYPE_I | TYPE_G)) return;
		) + "#endif" + R( // EQUILIBRIUM_BOUNDARIES
			float3 un = load_u(n, u); // cache velocity
		const float ul = length(un);
		if (def_scale_u * ul < 0.1f) return; // don't draw lattice points where the velocity is lower than this threshold
		float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
		for (uint i = 0u; i < 15u; i++) camera_cache[i] = camera[i];
		const float3 p = position(coordinates(n));
		const int c = rainbow_color(255.0f * def_scale_u * ul); // coloring by velocity
		draw_line(p - (0.5f / ul) * un, p + (0.5f / ul) * un, c, camera_cache, bitmap, zbuffer);
	}

	) + "#ifndef GRAPHICS_TEMPERATURE" + R(
	) + R(kernel void graphics_streamline(const global uchar * flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer) {
		) + "#else" + R( // GRAPHICS_TEMPERATURE
		) + R(kernel void graphics_streamline(const global uchar * flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer, const global float* T) {
			) + "#endif" + R( // GRAPHICS_TEMPERATURE
				const uint n = get_global_id(0);
			) + "#ifndef D2Q9" + R(
				if (n >= (def_Nx / def_streamline_sparse) * (def_Ny / def_streamline_sparse) * (def_Nz / def_streamline_sparse)) return;
			const uint z = n / ((def_Nx / def_streamline_sparse) * (def_Ny / def_streamline_sparse)); // disassemble 1D index to 3D coordinates
			const uint t = n % ((def_Nx / def_streamline_sparse) * (def_Ny / def_streamline_sparse));
			const uint y = t / (def_Nx / def_streamline_sparse);
			const uint x = t % (def_Nx / def_streamline_sparse);
			float3 p = (float)def_streamline_sparse * ((float3)((float)x + 0.5f, (float)y + 0.5f, (float)z + 0.5f)) - 0.5f * ((float3)((float)def_Nx, (float)def_Ny, (float)def_Nz));
			) + "#else" + R( // D2Q9
				if (n >= (def_Nx / def_streamline_sparse) * (def_Ny / def_streamline_sparse)) return;
			const uint y = n / (def_Nx / def_streamline_sparse); // disassemble 1D index to 3D coordinates
			const uint x = n % (def_Nx / def_streamline_sparse);
			float3 p = ((float3)((float)def_streamline_sparse * ((float)x + 0.5f), (float)def_streamline_sparse * ((float)y + 0.5f), 0.5f)) - 0.5f * ((float3)((float)def_Nx, (float)def_Ny, (float)def_Nz));
			) + "#endif" + R( // D2Q9
				float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
			for (uint i = 0u; i < 15u; i++) camera_cache[i] = camera[i];
			const float hLx = 0.5f * (float)(def_Nx - 2u * (def_Dx > 1u)), hLy = 0.5f * (float)(def_Ny - 2u * (def_Dy > 1u)), hLz = 0.5f * (float)(def_Nz - 2u * (def_Dz > 1u));
			//draw_circle(p, 0.5f*def_streamline_sparse, 0xFFFFFF, camera_cache, bitmap, zbuffer);
			for (float dt = -1.0f; dt <= 1.0f; dt += 2.0f) { // integrate forward and backward in time
				float3 p0, p1 = p;
				for (uint l = 0u; l < def_streamline_length / 2u; l++) {
					const uint x = (uint)(p1.x + 1.5f * (float)def_Nx) % def_Nx;
					const uint y = (uint)(p1.y + 1.5f * (float)def_Ny) % def_Ny;
					const uint z = (uint)(p1.z + 1.5f * (float)def_Nz) % def_Nz;
					const uint n = x + (y + z * def_Ny) * def_Nx;
					if (flags[n] & (TYPE_S | TYPE_E | TYPE_I | TYPE_G)) return;
					const float3 un = load_u(n, u); // interpolate_u(p1, u)
					const float ul = length(un);
					p0 = p1;
					p1 += (dt / ul) * un; // integrate forward in time
					if (def_scale_u * ul < 0.1f || p1.x<-hLx || p1.x>hLx || p1.y<-hLy || p1.y>hLy || p1.z<-hLz || p1.z>hLz) break;
					) + "#ifndef GRAPHICS_TEMPERATURE" + R(
						const int c = rainbow_color(255.0f * def_scale_u * ul);
					) + "#else" + R( // GRAPHICS_TEMPERATURE
						const int c = iron_color(167.0f + 255.0f * (T[n] - def_T_avg));
					) + "#endif" + R( // GRAPHICS_TEMPERATURE
						draw_line(p0, p1, c, camera_cache, bitmap, zbuffer);
				}
			}
		}

		) + R(kernel void graphics_q_field(const global uchar * flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer) {
			const uint n = get_global_id(0);
			if (n >= (uint)def_N || is_halo(n)) return; // don't execute graphics_q_field() on halo
			if (flags[n] & (TYPE_S | TYPE_E | TYPE_I | TYPE_G)) return;
			float3 un = load_u(n, u); // cache velocity
			const float ul = length(un);
			const float Q = calculate_Q(n, u);
			if (Q < def_scale_Q_min || ul == 0.0f) return; // don't draw lattice points where the velocity is very low
			float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
			for (uint i = 0u; i < 15u; i++) camera_cache[i] = camera[i];
			const float3 p = position(coordinates(n));
			const int c = rainbow_color(255.0f * def_scale_u * ul); // coloring by velocity
			draw_line(p - (0.5f / ul) * un, p + (0.5f / ul) * un, c, camera_cache, bitmap, zbuffer);
		}

		) + R(kernel void graphics_q(const global uchar * flags, const global float* u, const global float* camera, global int* bitmap, global int* zbuffer) {
			const uint n = get_global_id(0);
			if (is_halo_q(n)) return; // don't execute graphics_q_field() on marching-cubes halo
			const uint3 xyz = coordinates(n);
			if (xyz.x >= def_Nx - 1u || xyz.y >= def_Ny - 1u || xyz.z >= def_Nz - 1u) return;
			const uint x0 = xyz.x; // cube stencil
			const uint xp = xyz.x + 1u;
			const uint y0 = xyz.y * def_Nx;
			const uint yp = (xyz.y + 1u) * def_Nx;
			const uint z0 = xyz.z * def_Ny * def_Nx;
			const uint zp = (xyz.z + 1u) * def_Ny * def_Nx;
			const uint xq = (xyz.x + 2u) % def_Nx; // central difference stencil on each cube corner point
			const uint xm = (xyz.x + def_Nx - 1u) % def_Nx;
			const uint yq = ((xyz.y + 2u) % def_Ny) * def_Nx;
			const uint ym = ((xyz.y + def_Ny - 1u) % def_Ny) * def_Nx;
			const uint zq = ((xyz.z + 2u) % def_Nz) * def_Ny * def_Nx;
			const uint zm = ((xyz.z + def_Nz - 1u) % def_Nz) * def_Ny * def_Nx;
			uint j[32];
			j[0] = n; // 000 // cube stencil
			j[1] = xp + y0 + z0; // +00
			j[2] = xp + y0 + zp; // +0+
			j[3] = x0 + y0 + zp; // 00+
			j[4] = x0 + yp + z0; // 0+0
			j[5] = xp + yp + z0; // ++0
			j[6] = xp + yp + zp; // +++
			j[7] = x0 + yp + zp; // 0++
			j[8] = xm + y0 + z0; // -00 // central difference stencil on each cube corner point
			j[9] = x0 + ym + z0; // 0-0
			j[10] = x0 + y0 + zm; // 00-
			j[11] = xq + y0 + z0; // #00
			j[12] = xp + ym + z0; // +-0
			j[13] = xp + y0 + zm; // +0-
			j[14] = xq + y0 + zp; // #0+
			j[15] = xp + ym + zp; // +-+
			j[16] = xp + y0 + zq; // +0#
			j[17] = xm + y0 + zp; // -0+
			j[18] = x0 + ym + zp; // 0-+
			j[19] = x0 + y0 + zq; // 00#
			j[20] = xm + yp + z0; // -+0
			j[21] = x0 + yq + z0; // 0#0
			j[22] = x0 + yp + zm; // 0+-
			j[23] = xq + yp + z0; // #+0
			j[24] = xp + yq + z0; // +#0
			j[25] = xp + yp + zm; // ++-
			j[26] = xq + yp + zp; // #++
			j[27] = xp + yq + zp; // +#+
			j[28] = xp + yp + zq; // ++#
			j[29] = xm + yp + zp; // -++
			j[30] = x0 + yq + zp; // 0#+
			j[31] = x0 + yp + zq; // 0+#
			uchar flags_cell = 0u;
			for (uint i = 0u; i < 32u; i++) flags_cell |= flags[j[i]];
			if (flags_cell & (TYPE_S | TYPE_E | TYPE_I | TYPE_G)) return;
			float3 uj[8], u0[6], u1[6], u2[6], u3[6], u4[6], u5[6], u6[6], u7[6]; // don't load any velocity twice from global memory
			for (uint i = 0u; i < 8u; i++) uj[i] = load_u(j[i], u);
			u0[0] = uj[1]; u0[1] = load_u(j[8], u); u0[2] = uj[4]; u0[3] = load_u(j[9], u); u0[4] = uj[3]; u0[5] = load_u(j[10], u);
			u1[0] = load_u(j[11], u); u1[1] = uj[0]; u1[2] = uj[5]; u1[3] = load_u(j[12], u); u1[4] = uj[2]; u1[5] = load_u(j[13], u);
			u2[0] = load_u(j[14], u); u2[1] = uj[3]; u2[2] = uj[6]; u2[3] = load_u(j[15], u); u2[4] = load_u(j[16], u); u2[5] = uj[1];
			u3[0] = uj[2]; u3[1] = load_u(j[17], u); u3[2] = uj[7]; u3[3] = load_u(j[18], u); u3[4] = load_u(j[19], u); u3[5] = uj[0];
			u4[0] = uj[5]; u4[1] = load_u(j[20], u); u4[2] = load_u(j[21], u); u4[3] = uj[0]; u4[4] = uj[7]; u4[5] = load_u(j[22], u);
			u5[0] = load_u(j[23], u); u5[1] = uj[4]; u5[2] = load_u(j[24], u); u5[3] = uj[1]; u5[4] = uj[6]; u5[5] = load_u(j[25], u);
			u6[0] = load_u(j[26], u); u6[1] = uj[7]; u6[2] = load_u(j[27], u); u6[3] = uj[2]; u6[4] = load_u(j[28], u); u6[5] = uj[5];
			u7[0] = uj[6]; u7[1] = load_u(j[29], u); u7[2] = load_u(j[30], u); u7[3] = uj[3]; u7[4] = load_u(j[31], u); u7[5] = uj[4];
			float v[8];
			v[0] = calculate_Q_cached(u0);
			v[1] = calculate_Q_cached(u1);
			v[2] = calculate_Q_cached(u2);
			v[3] = calculate_Q_cached(u3);
			v[4] = calculate_Q_cached(u4);
			v[5] = calculate_Q_cached(u5);
			v[6] = calculate_Q_cached(u6);
			v[7] = calculate_Q_cached(u7);
			float3 triangles[15]; // maximum of 5 triangles with 3 vertices each
			const uint tn = marching_cubes(v, def_scale_Q_min, triangles); // run marching cubes algorithm
			if (tn == 0u) return;
			float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
			for (uint i = 0u; i < 15u; i++) camera_cache[i] = camera[i];
			const float3 offset = (float3)((float)xyz.x + 0.5f - 0.5f * (float)def_Nx, (float)xyz.y + 0.5f - 0.5f * (float)def_Ny, (float)xyz.z + 0.5f - 0.5f * (float)def_Nz);
			for (uint i = 0u; i < tn; i++) {
				const float3 p0 = triangles[3u * i]; // triangle coordinates in [0,1] (local cell)
				const float3 p1 = triangles[3u * i + 1u];
				const float3 p2 = triangles[3u * i + 2u];
				const float3 normal = cross(p1 - p0, p2 - p0);
				int c0, c1, c2;
				{
					const float x1 = p0.x, y1 = p0.y, z1 = p0.z, x0 = 1.0f - x1, y0 = 1.0f - y1, z0 = 1.0f - z1; // calculate interpolation factors
					const float3 ui = (x0 * y0 * z0) * uj[0] + (x1 * y0 * z0) * uj[1] + (x1 * y0 * z1) * uj[2] + (x0 * y0 * z1) * uj[3] + (x0 * y1 * z0) * uj[4] + (x1 * y1 * z0) * uj[5] + (x1 * y1 * z1) * uj[6] + (x0 * y1 * z1) * uj[7]; // perform trilinear interpolation
					c0 = lighting(rainbow_color(255.0f * def_scale_u * length(ui)), p0 + offset, normal, camera_cache); // rainbow_color(255.0f*def_scale_u*length(ui));
				} {
					const float x1 = p1.x, y1 = p1.y, z1 = p1.z, x0 = 1.0f - x1, y0 = 1.0f - y1, z0 = 1.0f - z1; // calculate interpolation factors
					const float3 ui = (x0 * y0 * z0) * uj[0] + (x1 * y0 * z0) * uj[1] + (x1 * y0 * z1) * uj[2] + (x0 * y0 * z1) * uj[3] + (x0 * y1 * z0) * uj[4] + (x1 * y1 * z0) * uj[5] + (x1 * y1 * z1) * uj[6] + (x0 * y1 * z1) * uj[7]; // perform trilinear interpolation
					c1 = lighting(rainbow_color(255.0f * def_scale_u * length(ui)), p1 + offset, normal, camera_cache); // rainbow_color(255.0f*def_scale_u*length(ui));
				} {
					const float x1 = p2.x, y1 = p2.y, z1 = p2.z, x0 = 1.0f - x1, y0 = 1.0f - y1, z0 = 1.0f - z1; // calculate interpolation factors
					const float3 ui = (x0 * y0 * z0) * uj[0] + (x1 * y0 * z0) * uj[1] + (x1 * y0 * z1) * uj[2] + (x0 * y0 * z1) * uj[3] + (x0 * y1 * z0) * uj[4] + (x1 * y1 * z0) * uj[5] + (x1 * y1 * z1) * uj[6] + (x0 * y1 * z1) * uj[7]; // perform trilinear interpolation
					c2 = lighting(rainbow_color(255.0f * def_scale_u * length(ui)), p2 + offset, normal, camera_cache); // rainbow_color(255.0f*def_scale_u*length(ui));
				}
				draw_triangle_interpolated(p0 + offset, p1 + offset, p2 + offset, c0, c1, c2, camera_cache, bitmap, zbuffer); // draw triangle with interpolated colors
			}
		}

		) + "#ifdef SURFACE" + R(
		) + R(kernel void graphics_rasterize_phi(const global float* phi, const global float* camera, global int* bitmap, global int* zbuffer) { // marching cubes
			const uint n = get_global_id(0);
			const uint3 xyz = coordinates(n);
			if (xyz.x >= def_Nx - 1u || xyz.y >= def_Ny - 1u || xyz.z >= def_Nz - 1u) return;
			uint j[8];
			const uint x0 = xyz.x; // cube stencil
			const uint xp = xyz.x + 1u;
			const uint y0 = xyz.y * def_Nx;
			const uint yp = (xyz.y + 1u) * def_Nx;
			const uint z0 = xyz.z * def_Ny * def_Nx;
			const uint zp = (xyz.z + 1u) * def_Ny * def_Nx;
			j[0] = n; // 000
			j[1] = xp + y0 + z0; // +00
			j[2] = xp + y0 + zp; // +0+
			j[3] = x0 + y0 + zp; // 00+
			j[4] = x0 + yp + z0; // 0+0
			j[5] = xp + yp + z0; // ++0
			j[6] = xp + yp + zp; // +++
			j[7] = x0 + yp + zp; // 0++
			float v[8];
			for (uint i = 0u; i < 8u; i++) v[i] = phi[j[i]];
			float3 triangles[15]; // maximum of 5 triangles with 3 vertices each
			const uint tn = marching_cubes(v, 0.5f, triangles); // run marching cubes algorithm
			if (tn == 0u) return;
			float camera_cache[15]; // cache camera parameters in case the kernel draws more than one shape
			for (uint i = 0u; i < 15u; i++) camera_cache[i] = camera[i];
			const float3 offset = (float3)((float)xyz.x + 0.5f - 0.5f * (float)def_Nx, (float)xyz.y + 0.5f - 0.5f * (float)def_Ny, (float)xyz.z + 0.5f - 0.5f * (float)def_Nz);
			for (uint i = 0u; i < tn; i++) {
				const float3 p0 = triangles[3u * i] + offset;
				const float3 p1 = triangles[3u * i + 1u] + offset;
				const float3 p2 = triangles[3u * i + 2u] + offset;
				const float3 p = (p0 + p1 + p2) / 3.0f, normal = cross(p1 - p0, p2 - p0);
				const int c = lighting(55 << 16 | 155 << 8 | 255, p, normal, camera_cache);
				draw_triangle(p0, p1, p2, c, camera_cache, bitmap, zbuffer);
				//draw_line(p0, p1, c, camera_cache, bitmap, zbuffer); // wireframe rendering
				//draw_line(p0, p2, c, camera_cache, bitmap, zbuffer);
				//draw_line(p1, p2, c, camera_cache, bitmap, zbuffer);
			}
		}

		) + R(int raytrace_phi_next_ray(const ray reflection, const ray transmission, const int pixelcolor, const float reflectivity, const global float* phi, const global uchar * flags, const global int* skybox) {
			int color_reflect = pixelcolor, color_transmit = pixelcolor;
			ray reflection_next, transmission_next;
			float reflection_reflectivity, transmission_reflectivity;
			if (raytrace_phi(reflection, &reflection_next, &transmission_next, &reflection_reflectivity, phi, flags, skybox, def_Nx, def_Ny, def_Nz)) {
				color_reflect = last_ray_reflectivity(reflection_next, transmission_next, color_reflect, reflection_reflectivity, skybox);
			}
			else {
				color_reflect = skybox_color(reflection, skybox);
			}
			if (raytrace_phi(transmission, &reflection_next, &transmission_next, &transmission_reflectivity, phi, flags, skybox, def_Nx, def_Ny, def_Nz)) {
				color_transmit = last_ray_reflectivity(reflection_next, transmission_next, color_transmit, transmission_reflectivity, skybox);
			}
			else {
				color_transmit = skybox_color(transmission, skybox);
			}
			return color_mix(color_reflect, color_transmit, reflectivity);
		}

		) + R(int raytrace_phi_next_ray_mirror(const ray reflection, const int pixelcolor, const global float* phi, const global uchar * flags, const global int* skybox) {
			int color_reflect = pixelcolor;
			ray reflection_next;
			if (raytrace_phi_mirror(reflection, &reflection_next, phi, flags, skybox, def_Nx, def_Ny, def_Nz)) {
				color_reflect = skybox_color(reflection_next, skybox);
			}
			else {
				color_reflect = skybox_color(reflection, skybox);
			}
			return color_reflect;
		}

		) + R(kernel void graphics_raytrace_phi(const global float* phi, const global uchar * flags, const global int* skybox, const global float* camera, global int* bitmap) { // marching cubes
			const uint gid = get_global_id(0); // workgroup size alignment is critical
			const uint lid = get_local_id(0); // make workgropus not horizontal stripes of pixels, but 8x8 rectangular (close to square) tiles
			const uint lsi = get_local_size(0); // (50% performance boost due to more coalesced memory access)
			const uint tile_width = 8u, tile_height = lsi / tile_width, tiles_x = def_screen_width / tile_width;
			const int lx = lid % tile_width, ly = lid / tile_width;
			const int tx = (gid / lsi) % tiles_x, ty = (gid / lsi) / tiles_x;
			const int x = tx * tile_width + lx, y = ty * tile_height + ly;
			const uint n = x + y * def_screen_width;
			float camera_cache[15]; // cache parameters in case the kernel draws more than one shape
			for (uint i = 0u; i < 15u; i++) camera_cache[i] = camera[i];
			ray camray = get_camray(x, y, camera_cache);
			int pixelcolor = 0;
			const float distance = intersect_cuboid(camray, (float3)(0.0f, 0.0f, 0.0f), (float)def_Nx, (float)def_Ny, (float)def_Nz);
			camray.origin = camray.origin + fmax(distance, 0.0f) * camray.direction;
			ray reflection, transmission; // reflection and transmission
			float reflectivity;
			if (raytrace_phi(camray, &reflection, &transmission, &reflectivity, phi, flags, skybox, def_Nx, def_Ny, def_Nz)) {
				pixelcolor = last_ray_reflectivity(reflection, transmission, pixelcolor, reflectivity, skybox); // 1 ray pass
				//pixelcolor = raytrace_phi_next_ray(reflection, transmission, pixelcolor, reflectivity, phi, flags, skybox); // 2 ray passes
			}
			else {
				pixelcolor = skybox_color(camray, skybox);
			}
			//if(raytrace_phi_mirror(camray, &reflection, phi, flags, skybox, def_Nx, def_Ny, def_Nz)) { // reflection only
			//	//pixelcolor = skybox_color(reflection, skybox); // 1 ray pass
			//	pixelcolor = raytrace_phi_next_ray_mirror(reflection, pixelcolor, phi, flags, skybox); // 2 ray passes
			//} else {
			//	pixelcolor = skybox_color(camray, skybox);
			//}
			bitmap[n] = pixelcolor; // no zbuffer required
		}
		) + "#endif" + R( // SURFACE

		) + "#ifdef PARTICLES" + R(
		) + R(kernel void graphics_particles(const global float* particles, const global float* camera, global int* bitmap, global int* zbuffer) {
			const uint n = get_global_id(0);
			if (n >= (uint)def_particles_N) return;
			float camera_cache[15]; // cache parameters in case the kernel draws more than one shape
			for (uint i = 0u; i < 15u; i++) camera_cache[i] = camera[i];
			const int c = COLOR_P; // coloring scheme
			const float3 p = (float3)(particles[n], particles[def_particles_N + (ulong)n], particles[2ul * def_particles_N + (ulong)n]);
			//draw_circle(p, 0.5f, c, camera_cache, bitmap, zbuffer);
			draw_point(p, c, camera_cache, bitmap, zbuffer);
		}
		) + "#endif" + R( // PARTICLES
		) + "#endif" + R( // GRAPHICS
	);
}