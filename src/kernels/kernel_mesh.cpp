#include "../kernel.hpp"



string kernel_mesh() {
	return R(

	) + R(kernel void voxelize_mesh(const uint direction, global fpxx * fi, global float* u, global uchar * flags, const ulong t, const uchar flag, const global float* p0, const global float* p1, const global float* p2, const global float* bbu) { // voxelize triangle mesh
		const uint a = get_global_id(0), A = get_area(direction); // a = domain area index for each side, A = area of the domain boundary
		if (a >= A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
		const uint triangle_number = as_uint(bbu[0]);
		const float x0 = bbu[1], y0 = bbu[2], z0 = bbu[3], x1 = bbu[4], y1 = bbu[5], z1 = bbu[6];
		const float cx = bbu[7], cy = bbu[8], cz = bbu[9], ux = bbu[10], uy = bbu[11], uz = bbu[12], rx = bbu[13], ry = bbu[14], rz = bbu[15];

		const uint3 xyz = direction == 0u ? (uint3)((uint)max(0, (int)x0 - def_Ox), a % def_Ny, a / def_Ny) : direction == 1u ? (uint3)(a / def_Nz, (uint)max(0, (int)y0 - def_Oy), a % def_Nz) : (uint3)(a % def_Nx, a / def_Nx, (uint)max(0, (int)z0 - def_Oz));
		const float3 p = position(xyz);
		const float3 offset = (float3)(0.5f * (float)((def_Nx - 2u * (def_Dx > 1u)) * def_Dx) - 0.5f, 0.5f * (float)((def_Ny - 2u * (def_Dy > 1u)) * def_Dy) - 0.5f, 0.5f * (float)((def_Nz - 2u * (def_Dz > 1u)) * def_Dz) - 0.5f) + (float3)(def_domain_offset_x, def_domain_offset_y, def_domain_offset_z);
		const float3 r_origin = p + offset;
		const float3 r_direction = (float3)((float)(direction == 0u), (float)(direction == 1u), (float)(direction == 2u));
		uint intersections = 0u, intersections_check = 0u;
		ushort distances[64]; // allow up to 64 mesh intersections
		const bool condition = direction == 0u ? r_origin.y < y0 || r_origin.z < z0 || r_origin.y >= y1 || r_origin.z >= z1 : direction == 1u ? r_origin.x < x0 || r_origin.z < z0 || r_origin.x >= x1 || r_origin.z >= z1 : r_origin.x < x0 || r_origin.y < y0 || r_origin.x >= x1 || r_origin.y >= y1;

		volatile local uint workgroup_condition; // use local memory optimization (~25% faster)
		workgroup_condition = 1u;
		barrier(CLK_LOCAL_MEM_FENCE);
		atomic_and(&workgroup_condition, (uint)condition);
		barrier(CLK_LOCAL_MEM_FENCE);
		const bool workgroup_all = (bool)workgroup_condition;
		if (workgroup_all) return; // return only if the entire workgroup is outside of the bounding-box of the mesh
		const uint lid = get_local_id(0);
		local float3 cache_p0[def_workgroup_size];
		local float3 cache_p1[def_workgroup_size];
		local float3 cache_p2[def_workgroup_size];
		for (uint i = 0u; i < triangle_number; i += def_workgroup_size) {
			const uint tx = 3u * (i + lid), ty = tx + 1u, tz = ty + 1u;
			cache_p0[lid] = (float3)(p0[tx], p0[ty], p0[tz]);
			cache_p1[lid] = (float3)(p1[tx], p1[ty], p1[tz]);
			cache_p2[lid] = (float3)(p2[tx], p2[ty], p2[tz]);
			barrier(CLK_LOCAL_MEM_FENCE);
			for (int j = 0; j < def_workgroup_size && i + j < triangle_number; j++) {
				const float3 p0i = cache_p0[j], p1i = cache_p1[j], p2i = cache_p2[j];
				const float3 u = p1i - p0i, v = p2i - p0i, w = r_origin - p0i, h = cross(r_direction, v), q = cross(w, u); // bidirectional ray-triangle intersection (Moeller-Trumbore algorithm)
				const float f = 1.0f / dot(u, h), s = f * dot(w, h), t = f * dot(r_direction, q), d = f * dot(v, q);
				if (s >= 0.0f && s < 1.0f && t >= 0.0f && s + t < 1.0f) { // ray-triangle intersection ahead or behind
					if (d > 0.0f) { // ray-triangle intersection ahead
						if (intersections < 64u && d < 65536.0f) distances[intersections] = (ushort)d; // store distance to intersection in array as ushort
						intersections++;
					}
					else { // ray-triangle intersection behind
						intersections_check++; // cast a second ray to check if starting point is really inside (error correction)
					}
				}
			}
			barrier(CLK_LOCAL_MEM_FENCE);
		}
		if (condition) return; // extra workgroup threads outside of the bounding-box are not needed anymore, so return /**/

		/*if(condition) return; // don't use local memory (this also runs on old OpenCL 1.0 GPUs)
		for(uint i=0u; i<triangle_number; i++) {
			const float3 p0i = (float3)(p0[3u*i], p0[3u*i+1u], p0[3u*i+2u]);
			const float3 p1i = (float3)(p1[3u*i], p1[3u*i+1u], p1[3u*i+2u]);
			const float3 p2i = (float3)(p2[3u*i], p2[3u*i+1u], p2[3u*i+2u]);
			const float3 u=p1i-p0i, v=p2i-p0i, w=r_origin-p0i, h=cross(r_direction, v), q=cross(w, u); // bidirectional ray-triangle intersection (Moeller-Trumbore algorithm)
			const float f=1.0f/dot(u, h), s=f*dot(w, h), t=f*dot(r_direction, q), d=f*dot(v, q);
			if(s>=0.0f&&s<1.0f&&t>=0.0f&&s+t<1.0f) { // ray-triangle intersection ahead or behind
				if(d>0.0f) { // ray-triangle intersection ahead
					if(intersections<64u&&d<65536.0f) distances[intersections] = (ushort)d; // store distance to intersection in array as ushort
					intersections++;
				} else { // ray-triangle intersection behind
					intersections_check++; // cast a second ray to check if starting point is really inside (error correction)
				}
			}
		}/**/

		for (int i = 1; i < (int)intersections; i++) { // insertion-sort distances
			ushort t = distances[i];
			int j = i - 1;
			while (distances[j] > t && j >= 0) {
				distances[j + 1] = distances[j];
				j--;
			}
			distances[j + 1] = t;
		}

		bool inside = (intersections % 2u) && (intersections_check % 2u);
		const bool set_u = sq(ux) + sq(uy) + sq(uz) + sq(rx) + sq(ry) + sq(rz) > 0.0f;
		uint intersection = intersections % 2u != intersections_check % 2u; // iterate through column, start with 0 regularly, start with 1 if forward and backward intersection count evenness differs (error correction)
		const uint h0 = direction == 0u ? xyz.x : direction == 1u ? xyz.y : xyz.z;
		const uint hmax = direction == 0u ? (uint)clamp((int)x1 - def_Ox, 0, (int)def_Nx - 1) : direction == 1u ? (uint)clamp((int)y1 - def_Oy, 0, (int)def_Ny - 1) : (uint)clamp((int)z1 - def_Oz, 0, (int)def_Nz - 1);
		const uint hmesh = h0 + (uint)distances[intersections - 1u];
		for (uint h = h0; h < hmax; h++) {
			while (intersection<intersections && h>h0 + (uint)distances[intersection]) {
				inside = !inside; // passed mesh intersection, so switch inside/outside state
				intersection++;
			}
			inside &= (intersection < intersections&& h < hmesh); // point must be outside if there are no more ray-mesh intersections ahead (error correction)
			const ulong n = index((uint3)(direction == 0u ? h : xyz.x, direction == 1u ? h : xyz.y, direction == 2u ? h : xyz.z));
			uchar flagsn = flags[n];
			if (inside) {
				flagsn = (flagsn & ~TYPE_BO) | flag;
				if (set_u) {
					const float3 p = position(coordinates(n)) + offset;
					const float3 un = (float3)(ux, uy, uz) + cross((float3)(cx, cy, cz) - p, (float3)(rx, ry, rz));
					u[n] = un.x;
					u[def_N + (ulong)n] = un.y;
					u[2ul * def_N + (ulong)n] = un.z;
				}
			}
			else {
				if (set_u) {
					const float3 un = (float3)(u[n], u[def_N + (ulong)n], u[2ul * def_N + (ulong)n]); // for velocity voxelization, only clear moving boundaries
					if ((flagsn & TYPE_BO) == TYPE_S) { // reconstruct DDFs when boundary point is converted to fluid
						uint j[def_velocity_set]; // neighbor indices
						neighbors(n, j); // calculate neighbor indices
						float feq[def_velocity_set]; // f_equilibrium
						calculate_f_eq(1.0f, un.x, un.y, un.z, feq);
						store_f(n, feq, fi, j, t); // write to fi
					}
					if (sq(un.x) + sq(un.y) + sq(un.z) > 0.0f) {
						flagsn = (flagsn & TYPE_BO) == TYPE_MS ? flagsn & ~TYPE_MS : flagsn & ~flag;
					}
				}
				else {
					flagsn = (flagsn & TYPE_BO) == TYPE_MS ? flagsn & ~TYPE_MS : flagsn & ~flag;
				}
			}
			flags[n] = flagsn;
		}
	} // voxelize_mesh()

	) + R(kernel void unvoxelize_mesh(global uchar * flags, const uchar flag, float x0, float y0, float z0, float x1, float y1, float z1) { // remove voxelized triangle mesh
		const uint n = get_global_id(0);
		const float3 p = position(coordinates(n)) + (float3)(0.5f * (float)((def_Nx - 2u * (def_Dx > 1u)) * def_Dx) - 0.5f, 0.5f * (float)((def_Ny - 2u * (def_Dy > 1u)) * def_Dy) - 0.5f, 0.5f * (float)((def_Nz - 2u * (def_Dz > 1u)) * def_Dz) - 0.5f) + (float3)(def_domain_offset_x, def_domain_offset_y, def_domain_offset_z);
		if (p.x >= x0 - 1.0f && p.y >= y0 - 1.0f && p.z >= z0 - 1.0f && p.x <= x1 + 1.0f && p.y <= y1 + 1.0f && p.z <= z1 + 1.0f) flags[n] &= ~flag;
	} // unvoxelize_mesh()

	);
}