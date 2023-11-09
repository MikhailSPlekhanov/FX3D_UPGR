#include "../kernel.hpp"

string kernel_transfer() {
	return R(

	) + R(uint get_area(const uint direction) {
		const uint A[3] = { def_Ax, def_Ay, def_Az };
		return A[direction];
	}
	) + R(uint index_extract_p(const uint a, const uint direction) {
		const uint3 coordinates[3] = { (uint3)(def_Nx - 2u, a % def_Ny, a / def_Ny), (uint3)(a / def_Nz, def_Ny - 2u, a % def_Nz), (uint3)(a % def_Nx, a / def_Nx, def_Nz - 2u) };
		return index(coordinates[direction]);
	}
	) + R(uint index_extract_m(const uint a, const uint direction) {
		const uint3 coordinates[3] = { (uint3)(1u, a % def_Ny, a / def_Ny), (uint3)(a / def_Nz,        1u, a % def_Nz), (uint3)(a % def_Nx, a / def_Nx,        1u) };
		return index(coordinates[direction]);
	}
	) + R(uint index_insert_p(const uint a, const uint direction) {
		const uint3 coordinates[3] = { (uint3)(def_Nx - 1u, a % def_Ny, a / def_Ny), (uint3)(a / def_Nz, def_Ny - 1u, a % def_Nz), (uint3)(a % def_Nx, a / def_Nx, def_Nz - 1u) };
		return index(coordinates[direction]);
	}
	) + R(uint index_insert_m(const uint a, const uint direction) {
		const uint3 coordinates[3] = { (uint3)(0u, a % def_Ny, a / def_Ny), (uint3)(a / def_Nz,        0u, a % def_Nz), (uint3)(a % def_Nx, a / def_Nx,        0u) };
		return index(coordinates[direction]);
	}

	) + R(uint index_transfer(const uint side_i) {
		const uchar index_transfer_data[2u * def_dimensions * def_transfers] = {
	) + "#if defined(D2Q9)" + R(
			1,  5,  7, // xp
			2,  6,  8, // xm
			3,  5,  8, // yp
			4,  6,  7  // ym
	) + "#elif defined(D3Q15)" + R(
			1,  7,  9, 11, 14, // xp
			2,  8, 10, 12, 13, // xm
			3,  7,  9, 12, 13, // yp
			4,  8, 10, 11, 14, // ym
			5,  7, 10, 11, 13, // zp
			6,  8,  9, 12, 14  // zm
	) + "#elif defined(D3Q19)" + R(
			1,  7,  9, 13, 15, // xp
			2,  8, 10, 14, 16, // xm
			3,  7, 11, 14, 17, // yp
			4,  8, 12, 13, 18, // ym
			5,  9, 11, 16, 18, // zp
			6, 10, 12, 15, 17  // zm
	) + "#elif defined(D3Q27)" + R(
			1,  7,  9, 13, 15, 19, 21, 23, 26, // xp
			2,  8, 10, 14, 16, 20, 22, 24, 25, // xm
			3,  7, 11, 14, 17, 19, 21, 24, 25, // yp
			4,  8, 12, 13, 18, 20, 22, 23, 26, // ym
			5,  9, 11, 16, 18, 19, 22, 23, 25, // zp
			6, 10, 12, 15, 17, 20, 21, 24, 26  // zm
	) + "#endif" + R( // D3Q27
		};
		return (uint)index_transfer_data[side_i];
	}
	) + R(void extract_fi(const uint a, const uint A, const uint n, const uint side, const ulong t, global fpxx_copy * transfer_buffer, const global fpxx_copy * fi) {
		uint j[def_velocity_set]; // neighbor indices
		neighbors(n, j); // calculate neighbor indices
		for (uint b = 0u; b < def_transfers; b++) {
			const uint i = index_transfer(side * def_transfers + b);
			const ulong index = index_f(i % 2u ? j[i] : n, t % 2ul ? (i % 2u ? i + 1u : i - 1u) : i); // Esoteric-Pull: standard store, or streaming part 1/2
			transfer_buffer[b * A + a] = fi[index]; // fpxx_copy allows direct copying without decompression+compression
		}
	}
	) + R(void insert_fi(const uint a, const uint A, const uint n, const uint side, const ulong t, const global fpxx_copy * transfer_buffer, global fpxx_copy * fi) {
		uint j[def_velocity_set]; // neighbor indices
		neighbors(n, j); // calculate neighbor indices
		for (uint b = 0u; b < def_transfers; b++) {
			const uint i = index_transfer(side * def_transfers + b);
			const ulong index = index_f(i % 2u ? n : j[i - 1u], t % 2ul ? i : (i % 2u ? i + 1u : i - 1u)); // Esoteric-Pull: standard load, or streaming part 2/2
			fi[index] = transfer_buffer[b * A + a]; // fpxx_copy allows direct copying without decompression+compression
		}
	}
	) + R(kernel void transfer_extract_fi(const uint direction, const ulong t, global fpxx_copy * transfer_buffer_p, global fpxx_copy * transfer_buffer_m, const global fpxx_copy * fi) {
		const uint a = get_global_id(0), A = get_area(direction); // a = domain area index for each side, A = area of the domain boundary
		if (a >= A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
		extract_fi(a, A, index_extract_p(a, direction), 2u * direction + 0u, t, transfer_buffer_p, fi);
		extract_fi(a, A, index_extract_m(a, direction), 2u * direction + 1u, t, transfer_buffer_m, fi);
	}
	) + R(kernel void transfer__insert_fi(const uint direction, const ulong t, const global fpxx_copy * transfer_buffer_p, const global fpxx_copy * transfer_buffer_m, global fpxx_copy * fi) {
		const uint a = get_global_id(0), A = get_area(direction); // a = domain area index for each side, A = area of the domain boundary
		if (a >= A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
		insert_fi(a, A, index_insert_p(a, direction), 2u * direction + 0u, t, transfer_buffer_p, fi);
		insert_fi(a, A, index_insert_m(a, direction), 2u * direction + 1u, t, transfer_buffer_m, fi);
	}

	) + R(void extract_rho_u_flags(const uint a, const uint A, const uint n, global char* transfer_buffer, const global float* rho, const global float* u, const global uchar * flags) {
		((global float*)transfer_buffer)[a] = rho[n];
		((global float*)transfer_buffer)[A + a] = u[n];
		((global float*)transfer_buffer)[2u * A + a] = u[def_N + (ulong)n];
		((global float*)transfer_buffer)[3u * A + a] = u[2ul * def_N + (ulong)n];
		((global uchar*)transfer_buffer)[16u * A + a] = flags[n];
	}
	) + R(void insert_rho_u_flags(const uint a, const uint A, const uint n, const global char* transfer_buffer, global float* rho, global float* u, global uchar * flags) {
		rho[n] = ((const global float*)transfer_buffer)[a];
		u[n] = ((const global float*)transfer_buffer)[A + a];
		u[def_N + (ulong)n] = ((const global float*)transfer_buffer)[2u * A + a];
		u[2ul * def_N + (ulong)n] = ((const global float*)transfer_buffer)[3u * A + a];
		flags[n] = ((const global uchar*)transfer_buffer)[16u * A + a];
	}
	) + R(kernel void transfer_extract_rho_u_flags(const uint direction, const ulong t, global char* transfer_buffer_p, global char* transfer_buffer_m, const global float* rho, const global float* u, const global uchar * flags) {
		const uint a = get_global_id(0), A = get_area(direction); // a = domain area index for each side, A = area of the domain boundary
		if (a >= A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
		extract_rho_u_flags(a, A, index_extract_p(a, direction), transfer_buffer_p, rho, u, flags);
		extract_rho_u_flags(a, A, index_extract_m(a, direction), transfer_buffer_m, rho, u, flags);
	}
	) + R(kernel void transfer__insert_rho_u_flags(const uint direction, const ulong t, const global char* transfer_buffer_p, const global char* transfer_buffer_m, global float* rho, global float* u, global uchar * flags) {
		const uint a = get_global_id(0), A = get_area(direction); // a = domain area index for each side, A = area of the domain boundary
		if (a >= A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
		insert_rho_u_flags(a, A, index_insert_p(a, direction), transfer_buffer_p, rho, u, flags);
		insert_rho_u_flags(a, A, index_insert_m(a, direction), transfer_buffer_m, rho, u, flags);
	}

	) + "#ifdef SURFACE" + R(
	) + R(kernel void transfer_extract_flags(const uint direction, const ulong t, global uchar * transfer_buffer_p, global uchar * transfer_buffer_m, const global uchar * flags) {
		const uint a = get_global_id(0), A = get_area(direction); // a = domain area index for each side, A = area of the domain boundary
		if (a >= A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
		transfer_buffer_p[a] = flags[index_extract_p(a, direction)];
		transfer_buffer_m[a] = flags[index_extract_m(a, direction)];
	}
	) + R(kernel void transfer__insert_flags(const uint direction, const ulong t, const global uchar * transfer_buffer_p, const global uchar * transfer_buffer_m, global uchar * flags) {
		const uint a = get_global_id(0), A = get_area(direction); // a = domain area index for each side, A = area of the domain boundary
		if (a >= A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
		flags[index_insert_p(a, direction)] = transfer_buffer_p[a];
		flags[index_insert_m(a, direction)] = transfer_buffer_m[a];
	}

	) + R(void extract_phi_massex_flags(const uint a, const uint A, const uint n, global char* transfer_buffer, const global float* phi, const global float* massex, const global uchar * flags) {
		((global float*)transfer_buffer)[a] = phi[n];
		((global float*)transfer_buffer)[A + a] = massex[n];
		((global uchar*)transfer_buffer)[8u * A + a] = flags[n];
	}
	) + R(void insert_phi_massex_flags(const uint a, const uint A, const uint n, const global char* transfer_buffer, global float* phi, global float* massex, global uchar * flags) {
		phi[n] = ((global float*)transfer_buffer)[a];
		massex[n] = ((global float*)transfer_buffer)[A + a];
		flags[n] = ((global uchar*)transfer_buffer)[8u * A + a];
	}
	) + R(kernel void transfer_extract_phi_massex_flags(const uint direction, const ulong t, global char* transfer_buffer_p, global char* transfer_buffer_m, const global float* phi, const global float* massex, const global uchar * flags) {
		const uint a = get_global_id(0), A = get_area(direction); // a = domain area index for each side, A = area of the domain boundary
		if (a >= A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
		extract_phi_massex_flags(a, A, index_extract_p(a, direction), transfer_buffer_p, phi, massex, flags);
		extract_phi_massex_flags(a, A, index_extract_m(a, direction), transfer_buffer_m, phi, massex, flags);
	}
	) + R(kernel void transfer__insert_phi_massex_flags(const uint direction, const ulong t, const global char* transfer_buffer_p, const global char* transfer_buffer_m, global float* phi, global float* massex, global uchar * flags) {
		const uint a = get_global_id(0), A = get_area(direction); // a = domain area index for each side, A = area of the domain boundary
		if (a >= A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
		insert_phi_massex_flags(a, A, index_insert_p(a, direction), transfer_buffer_p, phi, massex, flags);
		insert_phi_massex_flags(a, A, index_insert_m(a, direction), transfer_buffer_m, phi, massex, flags);
	}
	) + "#endif" + R( // SURFACE

	) + "#ifdef TEMPERATURE" + R(
	) + R(void extract_gi(const uint a, const uint n, const uint side, const ulong t, global fpxx_copy * transfer_buffer, const global fpxx_copy * gi) {
		uint j7[7u]; // neighbor indices
		neighbors_temperature(n, j7); // calculate neighbor indices
		const uint i = side + 1u;
		const ulong index = index_f(i % 2u ? j7[i] : n, t % 2ul ? (i % 2u ? i + 1u : i - 1u) : i); // Esoteric-Pull: standard store, or streaming part 1/2
		transfer_buffer[a] = gi[index]; // fpxx_copy allows direct copying without decompression+compression
	}
	) + R(void insert_gi(const uint a, const uint n, const uint side, const ulong t, const global fpxx_copy * transfer_buffer, global fpxx_copy * gi) {
		uint j7[7u]; // neighbor indices
		neighbors_temperature(n, j7); // calculate neighbor indices
		const uint i = side + 1u;
		const ulong index = index_f(i % 2u ? n : j7[i - 1u], t % 2ul ? i : (i % 2u ? i + 1u : i - 1u)); // Esoteric-Pull: standard load, or streaming part 2/2
		gi[index] = transfer_buffer[a]; // fpxx_copy allows direct copying without decompression+compression
	}
	) + R(kernel void transfer_extract_gi(const uint direction, const ulong t, global fpxx_copy * transfer_buffer_p, global fpxx_copy * transfer_buffer_m, const global fpxx_copy * gi) {
		const uint a = get_global_id(0), A = get_area(direction); // a = domain area index for each side, A = area of the domain boundary
		if (a >= A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
		extract_gi(a, index_extract_p(a, direction), 2u * direction + 0u, t, transfer_buffer_p, gi);
		extract_gi(a, index_extract_m(a, direction), 2u * direction + 1u, t, transfer_buffer_m, gi);
	}
	) + R(kernel void transfer__insert_gi(const uint direction, const ulong t, const global fpxx_copy * transfer_buffer_p, const global fpxx_copy * transfer_buffer_m, global fpxx_copy * gi) {
		const uint a = get_global_id(0), A = get_area(direction); // a = domain area index for each side, A = area of the domain boundary
		if (a >= A) return; // area might not be a multiple of def_workgroup_size, so return here to avoid writing in unallocated memory space
		insert_gi(a, index_insert_p(a, direction), 2u * direction + 0u, t, transfer_buffer_p, gi);
		insert_gi(a, index_insert_m(a, direction), 2u * direction + 1u, t, transfer_buffer_m, gi);
	}
	) + "#endif" + R( // TEMPERATURE

	);
}