#include "../kernel.hpp"



string kernel_utils() {
	return R( 

	) + R(float sq(const float x) {
		return x * x;
	}
	) + R(float cb(const float x) {
		return x * x * x;
	}
	) + R(float angle(const float3 v1, const float3 v2) {
		return acos(dot(v1, v2) / (length(v1) * length(v2)));
	}
	) + R(float fast_rsqrt(const float x) { // slightly fastwer approximation
		return as_float(0x5F37642F - (as_int(x) >> 1));
	}
	) + R(float fast_asin(const float x) { // slightly fastwer approximation
		return x * fma(0.5702f, sq(sq(sq(x))), 1.0f); // 0.5707964f = (pi-2)/2
	}
	) + R(float fast_acos(const float x) { // slightly fastwer approximation
		return fma(fma(-0.5702f, sq(sq(sq(x))), -1.0f), x, 1.5712963f); // 0.5707964f = (pi-2)/2
	}
	) + R(void lu_solve(float* M, float* x, float* b, const int N, const int Nsol) { // solves system of N linear equations M*x=b within dimensionality Nsol<=N
		for (int i = 0; i < Nsol; i++) { // decompose M in M=L*U
			for (int j = i + 1; j < Nsol; j++) {
				M[N * j + i] /= M[N * i + i];
				for (int k = i + 1; k < Nsol; k++) M[N * j + k] -= M[N * j + i] * M[N * i + k];
			}
		}
		for (int i = 0; i < Nsol; i++) { // find solution of L*y=b
			x[i] = b[i];
			for (int k = 0; k < i; k++) x[i] -= M[N * i + k] * x[k];
		}
		for (int i = Nsol - 1; i >= 0; i--) { // find solution of U*x=y
			for (int k = i + 1; k < Nsol; k++) x[i] -= M[N * i + k] * x[k];
			x[i] /= M[N * i + i];
		}
	}
	//bool workgroup_any(const bool condition) { // returns true if any thread within the workgroup enters true
//	volatile local uint workgroup_condition; // does not work on AMD GPUs (error: non-kernel function variable cannot be declared in local address space)
//	workgroup_condition = 0u;
//	barrier(CLK_LOCAL_MEM_FENCE);
//	atomic_or(&workgroup_condition, (uint)condition);
//	barrier(CLK_LOCAL_MEM_FENCE);
//	return (bool)workgroup_condition;
//}
//bool workgroup_all(const bool condition) { // returns true if all threads within the workgroup enter true
//	volatile local uint workgroup_condition; // does not work on AMD GPUs (error: non-kernel function variable cannot be declared in local address space)
//	workgroup_condition = 1u;
//	barrier(CLK_LOCAL_MEM_FENCE);
//	atomic_and(&workgroup_condition, (uint)condition);
//	barrier(CLK_LOCAL_MEM_FENCE);
//	return (bool)workgroup_condition;
//}
	) + R(void atomic_add_f(volatile global float* addr, const float val) { // not deterministic because the order of addition can vary: (a+b)+c is rounded differently than a+(b+c)
		union { // https://streamhpc.com/blog/2016-02-09/atomic-operations-for-floats-in-opencl-improved/
			uint  u32;
			float f32;
		} next, expected, current;
		current.f32 = *addr;
		do {
			next.f32 = (expected.f32 = current.f32) + val; // ...*val for atomic_mul_f()
			current.u32 = atomic_cmpxchg((volatile global uint*)addr, expected.u32, next.u32);
		} while (current.u32 != expected.u32);
	}
	//)+"#ifdef cl_khr_int64_base_atomics"+R( // OpenCL C defines don't work in R() stringification macro
	//)+R(void atomic_add_d(volatile global double* addr, const double val) { // not deterministic because the order of addition can vary: (a+b)+c is rounded differently than a+(b+c)
	//	union { // https://streamhpc.com/blog/2016-02-09/atomic-operations-for-floats-in-opencl-improved/
	//		ulong  u64;
	//		double f64;
	//	} next, expected, current;
	//	current.f64 = *addr;
	//	do {
	//		next.f64 = (expected.f64=current.f64)+val; // ...*val for atomic_mul_d()
	//		current.u64 = atom_cmpxchg((volatile global ulong*)addr, expected.u64, next.u64); // does not work on some older GPUs
	//	} while(current.u64!=expected.u64);
	//}
	//)+"#endif"+R( // cl_khr_int64_base_atomics



	// ################################################## Line3D code ##################################################

	// Line3D OpenCL C version 
	// draw_point(...)    : draw 3D pixel
	// draw_circle(...)   : draw 3D circle
	// draw_line(...)     : draw 3D line
	// draw_triangle(...) : draw 3D triangle
	// iron_color(...)    : convert float in [0,255] to iron spectrum int color
	// graphics_clear()   : kernel to reset bitmap and zbuffer
	);
}