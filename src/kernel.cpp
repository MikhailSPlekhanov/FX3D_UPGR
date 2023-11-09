#include "kernel.hpp" // note: unbalanced round brackets () are not allowed and string literals can't be arbitrarily long, so periodically interrupt with )+R(
string get_opencl_c_code() {
	string r =
		kernel_utils() +
		kernel_color_utils() +
		kernel_graphics_utils() +
		opencl_c_container() +
		kernel_transfer() +
		kernel_mesh() +
		kernel_graphics();
	r = replace(r, " ", "\n"); // replace all spaces by new lines
	r = replace(r, "#ifdef\n", "#ifdef "); // except for the arguments after some preprocessor options that need to be in the same line
	r = replace(r, "#ifndef\n", "#ifndef ");
	r = replace(r, "#define\n", "#define "); // #define with two arguments will not work
	r = replace(r, "#if\n", "#if "); // don't leave any spaces in arguments
	r = replace(r, "#elif\n", "#elif "); // don't leave any spaces in arguments
	r = replace(r, "#pragma\n", "#pragma ");
	return "\n" + r;
}
string opencl_c_container() { return R( // ########################## begin of OpenCL C code ####################################################################

// ################################################## LBM code ##################################################

)+R(uint3 coordinates(const uint n) { // disassemble 1D index to 3D coordinates (n -> x,y,z)
	const uint t = n%(def_Nx*def_Ny);
	return (uint3)(t%def_Nx, t/def_Nx, n/(def_Nx*def_Ny)); // n = x+(y+z*Ny)*Nx
}
)+R(uint index(const uint3 xyz) { // assemble 1D index from 3D coordinates (x,y,z -> n)
	return xyz.x+(xyz.y+xyz.z*def_Ny)*def_Nx; // n = x+(y+z*Ny)*Nx
}
)+R(float3 position(const uint3 xyz) { // 3D coordinates to 3D position
	return (float3)((float)xyz.x+0.5f-0.5f*(float)def_Nx, (float)xyz.y+0.5f-0.5f*(float)def_Ny, (float)xyz.z+0.5f-0.5f*(float)def_Nz);
}
)+R(float3 mirror_position(const float3 p) { // mirror position into periodic boundaries
	float3 r;
	r.x = sign(p.x)*(fmod(fabs(p.x)+0.5f*(float)def_Nx, (float)def_Nx)-0.5f*(float)def_Nx);
	r.y = sign(p.y)*(fmod(fabs(p.y)+0.5f*(float)def_Ny, (float)def_Ny)-0.5f*(float)def_Ny);
	r.z = sign(p.z)*(fmod(fabs(p.z)+0.5f*(float)def_Nz, (float)def_Nz)-0.5f*(float)def_Nz);
	return r;
}
)+R(float3 mirror_distance(const float3 d) { // mirror distance vector into periodic boundaries
	return mirror_position(d);
}
)+R(bool is_halo(const uint n) {
	const uint3 xyz = coordinates(n);
	return ((def_Dx>1u)&(xyz.x==0u||xyz.x>=def_Nx-1u))||((def_Dy>1u)&(xyz.y==0u||xyz.y>=def_Ny-1u))||((def_Dz>1u)&(xyz.z==0u||xyz.z>=def_Nz-1u));
}
)+R(bool is_halo_q(const uint n) {
	const uint3 xyz = coordinates(n);
)+"#ifndef SURFACE"+R(
	return ((def_Dx>1u)&(xyz.x<=1u||xyz.x>=def_Nx-3u))||((def_Dy>1u)&(xyz.y<=1u||xyz.y>=def_Ny-3u))||((def_Dz>1u)&(xyz.z<=1u||xyz.z>=def_Nz-3u));
)+"#else"+R( // SURFACE
	return ((def_Dx>1u)&(xyz.x==0u||xyz.x>=def_Nx-2u))||((def_Dy>0u)&(xyz.y==1u||xyz.y>=def_Ny-2u))||((def_Dz>1u)&(xyz.z==0u||xyz.z>=def_Nz-2u)); // halo data is kept up-to-date with SURFACE extension, so allow using halo data for rendering
)+"#endif"+R( // SURFACE
} // is_halo_q()

)+R(float half_to_float_custom(const ushort x) { // custom 16-bit floating-point format, 1-4-11, exp-15, +-1.99951168, +-6.10351562E-5, +-2.98023224E-8, 3.612 digits
	const uint e = (x&0x7800)>>11; // exponent
	const uint m = (x&0x07FF)<<12; // mantissa
	const uint v = as_uint((float)m)>>23; // evil log2 bit hack to count leading zeros in denormalized format
	return as_float((x&0x8000)<<16 | (e!=0)*((e+112)<<23|m) | ((e==0)&(m!=0))*((v-37)<<23|((m<<(150-v))&0x007FF000))); // sign : normalized : denormalized
}
)+R(ushort float_to_half_custom(const float x) { // custom 16-bit floating-point format, 1-4-11, exp-15, +-1.99951168, +-6.10351562E-5, +-2.98023224E-8, 3.612 digits
	const uint b = as_uint(x)+0x00000800; // round-to-nearest-even: add last bit after truncated mantissa
	const uint e = (b&0x7F800000)>>23; // exponent
	const uint m = b&0x007FFFFF; // mantissa; in line below: 0x007FF800 = 0x00800000-0x00000800 = decimal indicator flag - initial rounding
	return (b&0x80000000)>>16 | (e>112)*((((e-112)<<11)&0x7800)|m>>12) | ((e<113)&(e>100))*((((0x007FF800+m)>>(124-e))+1)>>1); // sign : normalized : denormalized (assume [-2,2])
}

)+R(ulong index_f(const uint n, const uint i) { // 64-bit indexing (maximum 2^32 lattice points (1624^3 lattice resolution, 225GB)
	return (ulong)i*def_N+(ulong)n; // SoA (229% faster on GPU)
}
)+R(float c(const uint i) { // avoid constant keyword by encapsulating data in function which gets inlined by compiler
	const float c[3u*def_velocity_set] = {
)+"#if defined(D2Q9)"+R(
		0, 1,-1, 0, 0, 1,-1, 1,-1, // x
		0, 0, 0, 1,-1, 1,-1,-1, 1, // y
		0, 0, 0, 0, 0, 0, 0, 0, 0  // z
)+"#elif defined(D3Q15)"+R(
		0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, // x
		0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1, // y
		0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1  // z
)+"#elif defined(D3Q19)"+R(
		0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, // x
		0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, // y
		0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1  // z
)+"#elif defined(D3Q27)"+R(
		0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, // x
		0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1, // y
		0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1  // z
)+"#endif"+R( // D3Q27
	};
	return c[i];
}
)+R(float w(const uint i) { // avoid constant keyword by encapsulating data in function which gets inlined by compiler
	const float w[def_velocity_set] = { def_w0, // velocity set weights
)+"#if defined(D2Q9)"+R(
		def_ws, def_ws, def_ws, def_ws, def_we, def_we, def_we, def_we
)+"#elif defined(D3Q15)"+R(
		def_ws, def_ws, def_ws, def_ws, def_ws, def_ws,
		def_wc, def_wc, def_wc, def_wc, def_wc, def_wc, def_wc, def_wc
)+"#elif defined(D3Q19)"+R(
		def_ws, def_ws, def_ws, def_ws, def_ws, def_ws,
		def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we
)+"#elif defined(D3Q27)"+R(
		def_ws, def_ws, def_ws, def_ws, def_ws, def_ws,
		def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we, def_we,
		def_wc, def_wc, def_wc, def_wc, def_wc, def_wc, def_wc, def_wc
)+"#endif"+R( // D3Q27
	};
	return w[i];
}
)+R(void calculate_indices(const uint n, uint* x0, uint* xp, uint* xm, uint* y0, uint* yp, uint* ym, uint* z0, uint* zp, uint* zm) {
	const uint3 xyz = coordinates(n);
	*x0 =   xyz.x; // pre-calculate indices (periodic boundary conditions)
	*xp =  (xyz.x       +1u)%def_Nx;
	*xm =  (xyz.x+def_Nx-1u)%def_Nx;
	*y0 =   xyz.y                   *def_Nx;
	*yp = ((xyz.y       +1u)%def_Ny)*def_Nx;
	*ym = ((xyz.y+def_Ny-1u)%def_Ny)*def_Nx;
	*z0 =   xyz.z                   *def_Ny*def_Nx;
	*zp = ((xyz.z       +1u)%def_Nz)*def_Ny*def_Nx;
	*zm = ((xyz.z+def_Nz-1u)%def_Nz)*def_Ny*def_Nx;
} // calculate_indices()
)+R(void neighbors(const uint n, uint* j) { // calculate neighbor indices
	uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
	calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
	j[0] = n;
)+"#if defined(D2Q9)"+R(
	j[ 1] = xp+y0; j[ 2] = xm+y0; // +00 -00
	j[ 3] = x0+yp; j[ 4] = x0+ym; // 0+0 0-0
	j[ 5] = xp+yp; j[ 6] = xm+ym; // ++0 --0
	j[ 7] = xp+ym; j[ 8] = xm+yp; // +-0 -+0
)+"#elif defined(D3Q15)"+R(
	j[ 1] = xp+y0+z0; j[ 2] = xm+y0+z0; // +00 -00
	j[ 3] = x0+yp+z0; j[ 4] = x0+ym+z0; // 0+0 0-0
	j[ 5] = x0+y0+zp; j[ 6] = x0+y0+zm; // 00+ 00-
	j[ 7] = xp+yp+zp; j[ 8] = xm+ym+zm; // +++ ---
	j[ 9] = xp+yp+zm; j[10] = xm+ym+zp; // ++- --+
	j[11] = xp+ym+zp; j[12] = xm+yp+zm; // +-+ -+-
	j[13] = xm+yp+zp; j[14] = xp+ym+zm; // -++ +--
)+"#elif defined(D3Q19)"+R(
	j[ 1] = xp+y0+z0; j[ 2] = xm+y0+z0; // +00 -00
	j[ 3] = x0+yp+z0; j[ 4] = x0+ym+z0; // 0+0 0-0
	j[ 5] = x0+y0+zp; j[ 6] = x0+y0+zm; // 00+ 00-
	j[ 7] = xp+yp+z0; j[ 8] = xm+ym+z0; // ++0 --0
	j[ 9] = xp+y0+zp; j[10] = xm+y0+zm; // +0+ -0-
	j[11] = x0+yp+zp; j[12] = x0+ym+zm; // 0++ 0--
	j[13] = xp+ym+z0; j[14] = xm+yp+z0; // +-0 -+0
	j[15] = xp+y0+zm; j[16] = xm+y0+zp; // +0- -0+
	j[17] = x0+yp+zm; j[18] = x0+ym+zp; // 0+- 0-+
)+"#elif defined(D3Q27)"+R(
	j[ 1] = xp+y0+z0; j[ 2] = xm+y0+z0; // +00 -00
	j[ 3] = x0+yp+z0; j[ 4] = x0+ym+z0; // 0+0 0-0
	j[ 5] = x0+y0+zp; j[ 6] = x0+y0+zm; // 00+ 00-
	j[ 7] = xp+yp+z0; j[ 8] = xm+ym+z0; // ++0 --0
	j[ 9] = xp+y0+zp; j[10] = xm+y0+zm; // +0+ -0-
	j[11] = x0+yp+zp; j[12] = x0+ym+zm; // 0++ 0--
	j[13] = xp+ym+z0; j[14] = xm+yp+z0; // +-0 -+0
	j[15] = xp+y0+zm; j[16] = xm+y0+zp; // +0- -0+
	j[17] = x0+yp+zm; j[18] = x0+ym+zp; // 0+- 0-+
	j[19] = xp+yp+zp; j[20] = xm+ym+zm; // +++ ---
	j[21] = xp+yp+zm; j[22] = xm+ym+zp; // ++- --+
	j[23] = xp+ym+zp; j[24] = xm+yp+zm; // +-+ -+-
	j[25] = xm+yp+zp; j[26] = xp+ym+zm; // -++ +--
)+"#endif"+R( // D3Q27
} // neighbors()

)+R(float3 load_u(const uint n, const global float* u) {
	return (float3)(u[n], u[def_N+(ulong)n], u[2ul*def_N+(ulong)n]);
}
)+R(float3 closest_u(const float3 p, const global float* u) { // return velocity of closest lattice point to point p
	const uint x = (uint)(p.x+1.5f*(float)def_Nx)%def_Nx;
	const uint y = (uint)(p.y+1.5f*(float)def_Ny)%def_Ny;
	const uint z = (uint)(p.z+1.5f*(float)def_Nz)%def_Nz;
	const uint n = x+(y+z*def_Ny)*def_Nx;
	return load_u(n, u);
} // closest_u()
)+R(float3 interpolate_u(const float3 p, const global float* u) { // trilinear interpolation of velocity at point p
	const float xa=p.x-0.5f+1.5f*def_Nx, ya=p.y-0.5f+1.5f*def_Ny, za=p.z-0.5f+1.5f*def_Nz; // subtract lattice offsets
	const uint xb=(uint)xa, yb=(uint)ya, zb=(uint)za; // integer casting to find bottom left corner
	const float x1=xa-(float)xb, y1=ya-(float)yb, z1=za-(float)zb, x0=1.0f-x1, y0=1.0f-y1, z0=1.0f-z1; // calculate interpolation factors
	float3 un[8]; // velocities of unit cube corner points
	for(uint c=0u; c<8u; c++) { // count over eight corner points
		const uint i=(c&0x04u)>>2, j=(c&0x02u)>>1, k=c&0x01u; // disassemble c into corner indices ijk
		const uint x=(xb+i)%def_Nx, y=(yb+j)%def_Ny, z=(zb+k)%def_Nz; // calculate corner lattice positions
		const uint n = x+(y+z*def_Ny)*def_Nx; // calculate lattice linear index
		un[c] = load_u(n, u); // load velocity from lattice point
	}
	return (x0*y0*z0)*un[0]+(x0*y0*z1)*un[1]+(x0*y1*z0)*un[2]+(x0*y1*z1)*un[3]+(x1*y0*z0)*un[4]+(x1*y0*z1)*un[5]+(x1*y1*z0)*un[6]+(x1*y1*z1)*un[7]; // perform trilinear interpolation
} // interpolate_u()
)+R(float calculate_Q_cached(const float3* uj) { // Q-criterion
	const float duxdx=uj[0].x-uj[1].x, duydx=uj[0].y-uj[1].y, duzdx=uj[0].z-uj[1].z; // du/dx = (u2-u0)/2
	const float duxdy=uj[2].x-uj[3].x, duydy=uj[2].y-uj[3].y, duzdy=uj[2].z-uj[3].z;
	const float duxdz=uj[4].x-uj[5].x, duydz=uj[4].y-uj[5].y, duzdz=uj[4].z-uj[5].z;
	const float omega_xy=duxdy-duydx, omega_xz=duxdz-duzdx, omega_yz=duydz-duzdy; // antisymmetric tensor, omega_xx = omega_yy = omega_zz = 0
	const float s_xx2=duxdx, s_yy2=duydy, s_zz2=duzdz; // s_xx2 = s_xx/2, s_yy2 = s_yy/2, s_zz2 = s_zz/2
	const float s_xy=duxdy+duydx, s_xz=duxdz+duzdx, s_yz=duydz+duzdy; // symmetric tensor
	const float omega2 = sq(omega_xy)+sq(omega_xz)+sq(omega_yz); // ||omega||_2^2
	const float s2 = 2.0f*(sq(s_xx2)+sq(s_yy2)+sq(s_zz2))+sq(s_xy)+sq(s_xz)+sq(s_yz); // ||s||_2^2
	return 0.25f*(omega2-s2); // Q = 1/2*(||omega||_2^2-||s||_2^2), addidional factor 1/2 from cental finite differences of velocity
} // calculate_Q_cached()
)+R(float calculate_Q(const uint n, const global float* u) { // Q-criterion
	uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
	calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
	uint j[6];
	j[0] = xp+y0+z0; j[1] = xm+y0+z0; // +00 -00
	j[2] = x0+yp+z0; j[3] = x0+ym+z0; // 0+0 0-0
	j[4] = x0+y0+zp; j[5] = x0+y0+zm; // 00+ 00-
	float3 uj[6];
	for(uint i=0u; i<6u; i++) uj[i] = load_u(j[i], u);
	return calculate_Q_cached(uj);
} // calculate_Q()

)+R(void calculate_f_eq(const float rho, float ux, float uy, float uz, float* feq) { // calculate f_equilibrium from density and velocity field (perturbation method / DDF-shifting)
	const float c3=-3.0f*(sq(ux)+sq(uy)+sq(uz)), rhom1=rho-1.0f; // c3 = -2*sq(u)/(2*sq(c)), rhom1 is arithmetic optimization to minimize digit extinction
	ux *= 3.0f;
	uy *= 3.0f;
	uz *= 3.0f;
	feq[ 0] = def_w0*fma(rho, 0.5f*c3, rhom1); // 000 (identical for all velocity sets)
)+"#if defined(D2Q9)"+R(
	const float u0=ux+uy, u1=ux-uy; // these pre-calculations make manual unrolling require less FLOPs
	const float rhos=def_ws*rho, rhoe=def_we*rho, rhom1s=def_ws*rhom1, rhom1e=def_we*rhom1;
	feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
	feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
	feq[ 5] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), u0), rhom1e); feq[ 6] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), -u0), rhom1e); // ++0 --0
	feq[ 7] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), u1), rhom1e); feq[ 8] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), -u1), rhom1e); // +-0 -+0
)+"#elif defined(D3Q15)"+R(
	const float u0=ux+uy+uz, u1=ux+uy-uz, u2=ux-uy+uz, u3=-ux+uy+uz;
	const float rhos=def_ws*rho, rhoc=def_wc*rho, rhom1s=def_ws*rhom1, rhom1c=def_wc*rhom1;
	feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
	feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
	feq[ 5] = fma(rhos, fma(0.5f, fma(uz, uz, c3), uz), rhom1s); feq[ 6] = fma(rhos, fma(0.5f, fma(uz, uz, c3), -uz), rhom1s); // 00+ 00-
	feq[ 7] = fma(rhoc, fma(0.5f, fma(u0, u0, c3), u0), rhom1c); feq[ 8] = fma(rhoc, fma(0.5f, fma(u0, u0, c3), -u0), rhom1c); // +++ ---
	feq[ 9] = fma(rhoc, fma(0.5f, fma(u1, u1, c3), u1), rhom1c); feq[10] = fma(rhoc, fma(0.5f, fma(u1, u1, c3), -u1), rhom1c); // ++- --+
	feq[11] = fma(rhoc, fma(0.5f, fma(u2, u2, c3), u2), rhom1c); feq[12] = fma(rhoc, fma(0.5f, fma(u2, u2, c3), -u2), rhom1c); // +-+ -+-
	feq[13] = fma(rhoc, fma(0.5f, fma(u3, u3, c3), u3), rhom1c); feq[14] = fma(rhoc, fma(0.5f, fma(u3, u3, c3), -u3), rhom1c); // -++ +--
)+"#elif defined(D3Q19)"+R(
	const float u0=ux+uy, u1=ux+uz, u2=uy+uz, u3=ux-uy, u4=ux-uz, u5=uy-uz;
	const float rhos=def_ws*rho, rhoe=def_we*rho, rhom1s=def_ws*rhom1, rhom1e=def_we*rhom1;
	feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
	feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
	feq[ 5] = fma(rhos, fma(0.5f, fma(uz, uz, c3), uz), rhom1s); feq[ 6] = fma(rhos, fma(0.5f, fma(uz, uz, c3), -uz), rhom1s); // 00+ 00-
	feq[ 7] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), u0), rhom1e); feq[ 8] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), -u0), rhom1e); // ++0 --0
	feq[ 9] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), u1), rhom1e); feq[10] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), -u1), rhom1e); // +0+ -0-
	feq[11] = fma(rhoe, fma(0.5f, fma(u2, u2, c3), u2), rhom1e); feq[12] = fma(rhoe, fma(0.5f, fma(u2, u2, c3), -u2), rhom1e); // 0++ 0--
	feq[13] = fma(rhoe, fma(0.5f, fma(u3, u3, c3), u3), rhom1e); feq[14] = fma(rhoe, fma(0.5f, fma(u3, u3, c3), -u3), rhom1e); // +-0 -+0
	feq[15] = fma(rhoe, fma(0.5f, fma(u4, u4, c3), u4), rhom1e); feq[16] = fma(rhoe, fma(0.5f, fma(u4, u4, c3), -u4), rhom1e); // +0- -0+
	feq[17] = fma(rhoe, fma(0.5f, fma(u5, u5, c3), u5), rhom1e); feq[18] = fma(rhoe, fma(0.5f, fma(u5, u5, c3), -u5), rhom1e); // 0+- 0-+
)+"#elif defined(D3Q27)"+R(
	const float u0=ux+uy, u1=ux+uz, u2=uy+uz, u3=ux-uy, u4=ux-uz, u5=uy-uz, u6=ux+uy+uz, u7=ux+uy-uz, u8=ux-uy+uz, u9=-ux+uy+uz;
	const float rhos=def_ws*rho, rhoe=def_we*rho, rhoc=def_wc*rho, rhom1s=def_ws*rhom1, rhom1e=def_we*rhom1, rhom1c=def_wc*rhom1;
	feq[ 1] = fma(rhos, fma(0.5f, fma(ux, ux, c3), ux), rhom1s); feq[ 2] = fma(rhos, fma(0.5f, fma(ux, ux, c3), -ux), rhom1s); // +00 -00
	feq[ 3] = fma(rhos, fma(0.5f, fma(uy, uy, c3), uy), rhom1s); feq[ 4] = fma(rhos, fma(0.5f, fma(uy, uy, c3), -uy), rhom1s); // 0+0 0-0
	feq[ 5] = fma(rhos, fma(0.5f, fma(uz, uz, c3), uz), rhom1s); feq[ 6] = fma(rhos, fma(0.5f, fma(uz, uz, c3), -uz), rhom1s); // 00+ 00-
	feq[ 7] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), u0), rhom1e); feq[ 8] = fma(rhoe, fma(0.5f, fma(u0, u0, c3), -u0), rhom1e); // ++0 --0
	feq[ 9] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), u1), rhom1e); feq[10] = fma(rhoe, fma(0.5f, fma(u1, u1, c3), -u1), rhom1e); // +0+ -0-
	feq[11] = fma(rhoe, fma(0.5f, fma(u2, u2, c3), u2), rhom1e); feq[12] = fma(rhoe, fma(0.5f, fma(u2, u2, c3), -u2), rhom1e); // 0++ 0--
	feq[13] = fma(rhoe, fma(0.5f, fma(u3, u3, c3), u3), rhom1e); feq[14] = fma(rhoe, fma(0.5f, fma(u3, u3, c3), -u3), rhom1e); // +-0 -+0
	feq[15] = fma(rhoe, fma(0.5f, fma(u4, u4, c3), u4), rhom1e); feq[16] = fma(rhoe, fma(0.5f, fma(u4, u4, c3), -u4), rhom1e); // +0- -0+
	feq[17] = fma(rhoe, fma(0.5f, fma(u5, u5, c3), u5), rhom1e); feq[18] = fma(rhoe, fma(0.5f, fma(u5, u5, c3), -u5), rhom1e); // 0+- 0-+
	feq[19] = fma(rhoc, fma(0.5f, fma(u6, u6, c3), u6), rhom1c); feq[20] = fma(rhoc, fma(0.5f, fma(u6, u6, c3), -u6), rhom1c); // +++ ---
	feq[21] = fma(rhoc, fma(0.5f, fma(u7, u7, c3), u7), rhom1c); feq[22] = fma(rhoc, fma(0.5f, fma(u7, u7, c3), -u7), rhom1c); // ++- --+
	feq[23] = fma(rhoc, fma(0.5f, fma(u8, u8, c3), u8), rhom1c); feq[24] = fma(rhoc, fma(0.5f, fma(u8, u8, c3), -u8), rhom1c); // +-+ -+-
	feq[25] = fma(rhoc, fma(0.5f, fma(u9, u9, c3), u9), rhom1c); feq[26] = fma(rhoc, fma(0.5f, fma(u9, u9, c3), -u9), rhom1c); // -++ +--
)+"#endif"+R( // D3Q27
} // calculate_f_eq()

)+R(void calculate_rho_u(const float* f, float* rhon, float* uxn, float* uyn, float* uzn) { // calculate density and velocity fields from fi
	float rho=f[0], ux, uy, uz;
	for(uint i=1u; i<def_velocity_set; i++) rho += f[i]; // calculate density from fi
	rho += 1.0f; // add 1.0f last to avoid digit extinction effects when summing up fi (perturbation method / DDF-shifting)
)+"#if defined(D2Q9)"+R(
	ux = f[1]-f[2]+f[5]-f[6]+f[7]-f[8]; // calculate velocity from fi (alternating + and - for best accuracy)
	uy = f[3]-f[4]+f[5]-f[6]+f[8]-f[7];
	uz = 0.0f;
)+"#elif defined(D3Q15)"+R(
	ux = f[ 1]-f[ 2]+f[ 7]-f[ 8]+f[ 9]-f[10]+f[11]-f[12]+f[14]-f[13]; // calculate velocity from fi (alternating + and - for best accuracy)
	uy = f[ 3]-f[ 4]+f[ 7]-f[ 8]+f[ 9]-f[10]+f[12]-f[11]+f[13]-f[14];
	uz = f[ 5]-f[ 6]+f[ 7]-f[ 8]+f[10]-f[ 9]+f[11]-f[12]+f[13]-f[14];
)+"#elif defined(D3Q19)"+R(
	ux = f[ 1]-f[ 2]+f[ 7]-f[ 8]+f[ 9]-f[10]+f[13]-f[14]+f[15]-f[16]; // calculate velocity from fi (alternating + and - for best accuracy)
	uy = f[ 3]-f[ 4]+f[ 7]-f[ 8]+f[11]-f[12]+f[14]-f[13]+f[17]-f[18];
	uz = f[ 5]-f[ 6]+f[ 9]-f[10]+f[11]-f[12]+f[16]-f[15]+f[18]-f[17];
)+"#elif defined(D3Q27)"+R(
	ux = f[ 1]-f[ 2]+f[ 7]-f[ 8]+f[ 9]-f[10]+f[13]-f[14]+f[15]-f[16]+f[19]-f[20]+f[21]-f[22]+f[23]-f[24]+f[26]-f[25]; // calculate velocity from fi (alternating + and - for best accuracy)
	uy = f[ 3]-f[ 4]+f[ 7]-f[ 8]+f[11]-f[12]+f[14]-f[13]+f[17]-f[18]+f[19]-f[20]+f[21]-f[22]+f[24]-f[23]+f[25]-f[26];
	uz = f[ 5]-f[ 6]+f[ 9]-f[10]+f[11]-f[12]+f[16]-f[15]+f[18]-f[17]+f[19]-f[20]+f[22]-f[21]+f[23]-f[24]+f[25]-f[26];
)+"#endif"+R( // D3Q27
	*rhon = rho;
	*uxn = ux/rho;
	*uyn = uy/rho;
	*uzn = uz/rho;
} // calculate_rho_u()

)+"#ifdef VOLUME_FORCE"+R(
)+R(void calculate_forcing_terms(const float ux, const float uy, const float uz, const float fx, const float fy, const float fz, float* Fin) { // calculate volume force terms Fin from velocity field (Guo forcing, Krueger p.233f)
)+"#ifdef D2Q9"+R(
	const float uF = -0.33333334f*fma(ux, fx, uy*fy); // 2D
)+"#else"+R( // D2Q9
	const float uF = -0.33333334f*fma(ux, fx, fma(uy, fy, uz*fz)); // 3D
)+"#endif"+R( // D2Q9
	Fin[0] = 9.0f*def_w0*uF ; // 000 (identical for all velocity sets)
	for(uint i=1u; i<def_velocity_set; i++) { // loop is entirely unrolled by compiler, no unnecessary FLOPs are happening
		Fin[i] = 9.0f*w(i)*fma(c(i)*fx+c(def_velocity_set+i)*fy+c(2u*def_velocity_set+i)*fz, c(i)*ux+c(def_velocity_set+i)*uy+c(2u*def_velocity_set+i)*uz+0.33333334f, uF);
	}
} // calculate_forcing_terms()
)+"#endif"+R( // VOLUME_FORCE

)+"#ifdef MOVING_BOUNDARIES"+R(
)+R(void apply_moving_boundaries(float* fhn, const uint* j, const global float* u, const global uchar* flags) { // apply Dirichlet velocity boundaries if necessary (Krueger p.180, rho_solid=1)
	uint ji; // reads velocities of only neighboring boundary nodes, which do not change during simulation
	for(uint i=1u; i<def_velocity_set; i+=2u) { // loop is entirely unrolled by compiler, no unnecessary memory access is happening
		const float w6 = -6.0f*w(i); // w(i) = w(i+1) if i is odd
		ji = j[i+1u]; fhn[i   ] = (flags[ji]&TYPE_BO)==TYPE_S ? fma(w6, c(i+1u)*u[ji]+c(def_velocity_set+i+1u)*u[def_N+(ulong)ji]+c(2u*def_velocity_set+i+1u)*u[2ul*def_N+(ulong)ji], fhn[i   ]) : fhn[i   ]; // boundary : regular
		ji = j[i   ]; fhn[i+1u] = (flags[ji]&TYPE_BO)==TYPE_S ? fma(w6, c(i   )*u[ji]+c(def_velocity_set+i   )*u[def_N+(ulong)ji]+c(2u*def_velocity_set+i   )*u[2ul*def_N+(ulong)ji], fhn[i+1u]) : fhn[i+1u];
	}
} // apply_moving_boundaries()
)+"#endif"+R( // MOVING_BOUNDARIES

)+"#ifdef SURFACE"+R(
)+R(void average_neighbors_non_gas(const uint n, const global float* rho, const global float* u, const global uchar* flags, float* rhon, float* uxn, float* uyn, float* uzn) { // calculate average density and velocity of neighbors of node n
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	float rhot=0.0f, uxt=0.0f, uyt=0.0f, uzt=0.0f, counter=0.0f; // average over all fluid/interface neighbors
	for(uint i=1u; i<def_velocity_set; i++) {
		const uchar flagsji_sus = flags[j[i]]&(TYPE_SU|TYPE_S); // extract SURFACE flags
		if(flagsji_sus==TYPE_F||flagsji_sus==TYPE_I||flagsji_sus==TYPE_IF) { // fluid or interface or (interface->fluid) neighbor
			counter += 1.0f;
			rhot += rho[               j[i]];
			uxt  += u[                 j[i]];
			uyt  += u[    def_N+(ulong)j[i]];
			uzt  += u[2ul*def_N+(ulong)j[i]];
		}
	}
	*rhon = counter>0.0f ? rhot/counter : 1.0f;
	*uxn  = counter>0.0f ? uxt /counter : 0.0f;
	*uyn  = counter>0.0f ? uyt /counter : 0.0f;
	*uzn  = counter>0.0f ? uzt /counter : 0.0f;
}
)+R(void average_neighbors_fluid(const uint n, const global float* rho, const global float* u, const global uchar* flags, float* rhon, float* uxn, float* uyn, float* uzn) { // calculate average density and velocity of neighbors of node n
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	float rhot=0.0f, uxt=0.0f, uyt=0.0f, uzt=0.0f, counter=0.0f; // average over all fluid/interface neighbors
	for(uint i=1u; i<def_velocity_set; i++) {
		const uchar flagsji_su = flags[j[i]]&TYPE_SU;
		if(flagsji_su==TYPE_F) { // fluid neighbor
			counter += 1.0f;
			rhot += rho[               j[i]];
			uxt  += u[                 j[i]];
			uyt  += u[    def_N+(ulong)j[i]];
			uzt  += u[2ul*def_N+(ulong)j[i]];
		}
	}
	*rhon = counter>0.0f ? rhot/counter : 1.0f;
	*uxn  = counter>0.0f ? uxt /counter : 0.0f;
	*uyn  = counter>0.0f ? uyt /counter : 0.0f;
	*uzn  = counter>0.0f ? uzt /counter : 0.0f;
}
)+R(float calculate_phi(const float rhon, const float massn, const uchar flagsn) { // calculate fill level
	return flagsn&TYPE_F ? 1.0f : flagsn&TYPE_I ? rhon>0.0f ? clamp(massn/rhon, 0.0f, 1.0f) : 0.5f : 0.0f;
}
)+R(float3 calculate_normal_py(const float* phij) { // calculate surface normal vector (Parker-youngs approximation, more accurate, works only for D3Q27 neighborhood)
	float3 n; // normal vector
)+"#ifdef D2Q9"+R(
	n.x = 2.0f*(phij[2]-phij[1])+phij[6]-phij[5]+phij[8]-phij[7];
	n.y = 2.0f*(phij[4]-phij[3])+phij[6]-phij[5]+phij[7]-phij[8];
	n.z = 0.0f;
)+"#else"+R( // D2Q9
	n.x = 4.0f*(phij[ 2]-phij[ 1])+2.0f*(phij[ 8]-phij[ 7]+phij[10]-phij[ 9]+phij[14]-phij[13]+phij[16]-phij[15])+phij[20]-phij[19]+phij[22]-phij[21]+phij[24]-phij[23]+phij[25]-phij[26];
	n.y = 4.0f*(phij[ 4]-phij[ 3])+2.0f*(phij[ 8]-phij[ 7]+phij[12]-phij[11]+phij[13]-phij[14]+phij[18]-phij[17])+phij[20]-phij[19]+phij[22]-phij[21]+phij[23]-phij[24]+phij[26]-phij[25];
	n.z = 4.0f*(phij[ 6]-phij[ 5])+2.0f*(phij[10]-phij[ 9]+phij[12]-phij[11]+phij[15]-phij[16]+phij[17]-phij[18])+phij[20]-phij[19]+phij[21]-phij[22]+phij[24]-phij[23]+phij[26]-phij[25];
)+"#endif"+R( // D2Q9
	return normalize(n);
}
)+R(float plic_cube_reduced(const float V, const float n1, const float n2, const float n3) { // optimized solution from SZ and Kawano, source: https://doi.org/10.3390/computation10020021
	const float n12=n1+n2, n3V=n3*V;
	if(n12<=2.0f*n3V) return n3V+0.5f*n12; // case (5)
	const float sqn1=sq(n1), n26=6.0f*n2, v1=sqn1/n26; // after case (5) check n2>0 is true
	if(v1<=n3V && n3V<v1+0.5f*(n2-n1)) return 0.5f*(n1+sqrt(sqn1+8.0f*n2*(n3V-v1))); // case (2)
	const float V6 = n1*n26*n3V;
	if(n3V<v1) return cbrt(V6); // case (1)
	const float v3 = n3<n12 ? (sq(n3)*(3.0f*n12-n3)+sqn1*(n1-3.0f*n3)+sq(n2)*(n2-3.0f*n3))/(n1*n26) : 0.5f*n12; // after case (2) check n1>0 is true
	const float sqn12=sqn1+sq(n2), V6cbn12=V6-cb(n1)-cb(n2);
	const bool case34 = n3V<v3; // true: case (3), false: case (4)
	const float a = case34 ? V6cbn12 : 0.5f*(V6cbn12-cb(n3));
	const float b = case34 ?   sqn12 : 0.5f*(sqn12+sq(n3));
	const float c = case34 ?     n12 : 0.5f;
	const float t = sqrt(sq(c)-b);
	return c-2.0f*t*sin(0.33333334f*asin((cb(c)-0.5f*a-1.5f*b*c)/cb(t)));
}
)+R(float plic_cube(const float V0, const float3 n) { // unit cube - plane intersection: volume V0 in [0,1], normal vector n -> plane offset d0
	const float ax=fabs(n.x), ay=fabs(n.y), az=fabs(n.z), V=0.5f-fabs(V0-0.5f), l=ax+ay+az; // eliminate symmetry cases, normalize n using L1 norm
	const float n1 = fmin(fmin(ax, ay), az)/l;
	const float n3 = fmax(fmax(ax, ay), az)/l;
	const float n2 = fdim(1.0f, n1+n3); // ensure n2>=0
	const float d = plic_cube_reduced(V, n1, n2, n3); // calculate PLIC with reduced symmetry
	return l*copysign(0.5f-d, V0-0.5f); // rescale result and apply symmetry for V0>0.5
}
)+R(void get_remaining_neighbor_phij(const uint n, const float* phit, const global float* phi, float* phij) { // get remaining phij for D3Q27 neighborhood
)+"#ifndef D3Q27"+R(
	uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
	calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
)+"#endif"+R( // D3Q27
)+"#if defined(D3Q15)"+R(
	uint j[12]; // calculate neighbor indices
	j[ 0] = xp+yp+z0; j[ 1] = xm+ym+z0; // ++0 --0
	j[ 2] = xp+y0+zp; j[ 3] = xm+y0+zm; // +0+ -0-
	j[ 4] = x0+yp+zp; j[ 5] = x0+ym+zm; // 0++ 0--
	j[ 6] = xp+ym+z0; j[ 7] = xm+yp+z0; // +-0 -+0
	j[ 8] = xp+y0+zm; j[ 9] = xm+y0+zp; // +0- -0+
	j[10] = x0+yp+zm; j[11] = x0+ym+zp; // 0+- 0-+
	for(uint i=0u; i<7u; i++) phij[i] = phit[i];
	for(uint i=7u; i<19u; i++) phij[i] = phi[j[i-7u]];
	for(uint i=19u; i<27u; i++) phij[i] = phit[i-12u];
)+"#elif defined(D3Q19)"+R(
	uint j[8]; // calculate remaining neighbor indices
	j[0] = xp+yp+zp; j[1] = xm+ym+zm; // +++ ---
	j[2] = xp+yp+zm; j[3] = xm+ym+zp; // ++- --+
	j[4] = xp+ym+zp; j[5] = xm+yp+zm; // +-+ -+-
	j[6] = xm+yp+zp; j[7] = xp+ym+zm; // -++ +--
	for(uint i=0u; i<19u; i++) phij[i] = phit[i];
	for(uint i=19u; i<27u; i++) phij[i] = phi[j[i-19u]];
)+"#elif defined(D3Q27)"+R(
	for(uint i=0u; i<def_velocity_set; i++) phij[i] = phit[i];
)+"#endif"+R( // D3Q27
}
)+R(float c_D3Q27(const uint i) { // avoid constant keyword by encapsulating data in function which gets inlined by compiler
	const float c[3*27] = {
		0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, // x
		0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1, // y
		0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1  // z
	};
	return c[i];
}
)+R(float curvature_calculation(const uint n, const float* phit, const global float* phi) { // calculate surface curvature, always use D3Q27 stencil here, source: https://doi.org/10.3390/computation10020021
)+"#ifndef D2Q9"+R(
	float phij[27];
	get_remaining_neighbor_phij(n, phit, phi, phij); // complete neighborhood from whatever velocity set is selected to D3Q27
	const float3 bz = calculate_normal_py(phij); // new coordinate system: bz is normal to surface, bx and by are tangent to surface
	const float3 rn = (float3)(0.56270900f, 0.32704452f, 0.75921047f); // random normalized vector that is just by random chance not collinear with bz
	const float3 by = normalize(cross(bz, rn)); // normalize() is necessary here because bz and rn are not perpendicular
	const float3 bx = cross(by, bz);
	uint number = 0; // number of neighboring interface points
	float3 p[24]; // number of neighboring interface points is less or equal than than 26 minus 1 gas and minus 1 fluid point = 24
	const float center_offset = plic_cube(phij[0], bz); // calculate z-offset PLIC of center point only once
	for(uint i=1u; i<27u; i++) { // iterate over neighbors, no loop unrolling here (50% better perfoemance without loop unrolling)
		if(phij[i]>0.0f&&phij[i]<1.0f) { // limit neighbors to interface nodes
			const float3 ei = (float3)(c_D3Q27(i), c_D3Q27(27u+i), c_D3Q27(2u*27u+i)); // assume neighbor normal vector is the same as center normal vector
			const float offset = plic_cube(phij[i], bz)-center_offset;
			p[number++] = (float3)(dot(ei, bx), dot(ei, by), dot(ei, bz)+offset); // do coordinate system transformation into (x, y, f(x,y)) and apply PLIC pffsets
		}
	}
	float M[25], x[5]={0.0f,0.0f,0.0f,0.0f,0.0f}, b[5]={0.0f,0.0f,0.0f,0.0f,0.0f};
	for(uint i=0u; i<25u; i++) M[i] = 0.0f;
	for(uint i=0u; i<number; i++) { // f(x,y)=A*x2+B*y2+C*x*y+H*x+I*y, x=(A,B,C,H,I), Q=(x2,y2,x*y,x,y), M*x=b, M=Q*Q^T, b=Q*z
		const float x=p[i].x, y=p[i].y, z=p[i].z, x2=x*x, y2=y*y, x3=x2*x, y3=y2*y;
		/**/M[ 0]+=x2*x2; M[ 1]+=x2*y2; M[ 2]+=x3*y ; M[ 3]+=x3   ; M[ 4]+=x2*y ; b[0]+=x2   *z;
		/*M[ 5]+=x2*y2;*/ M[ 6]+=y2*y2; M[ 7]+=x *y3; M[ 8]+=x *y2; M[ 9]+=   y3; b[1]+=   y2*z;
		/*M[10]+=x3*y ; M[11]+=x *y3;*/ M[12]+=x2*y2; M[13]+=x2*y ; M[14]+=x *y2; b[2]+=x *y *z;
		/*M[15]+=x3   ; M[16]+=x *y2; M[17]+=x2*y ;*/ M[18]+=x2   ; M[19]+=x *y ; b[3]+=x    *z;
		/*M[20]+=x2*y ; M[21]+=   y3; M[22]+=x *y2; M[23]+=x *y ;*/ M[24]+=   y2; b[4]+=   y *z;
	}
	for(uint i=1u; i<5u; i++) { // use symmetry of matrix to save arithmetic operations
		for(uint j=0u; j<i; j++) M[i*5+j] = M[j*5+i];
	}
	if(number>=5u) lu_solve(M, x, b, 5, 5);
	else lu_solve(M, x, b, 5, min(5u, number)); // cannot do loop unrolling here -> slower -> extra if-else to avoid slowdown
	const float A=x[0], B=x[1], C=x[2], H=x[3], I=x[4];
	const float K = (A*(I*I+1.0f)+B*(H*H+1.0f)-C*H*I)*cb(rsqrt(H*H+I*I+1.0f)); // mean curvature of Monge patch (x, y, f(x, y))
)+"#else"+R( // D2Q9
	const float3 by = calculate_normal_py(phit); // new coordinate system: bz is normal to surface, bx and by are tangent to surface
	const float3 bx = cross(by, (float3)(0.0f, 0.0f, 1.0f)); // normalize() is necessary here because bz and rn are not perpendicular
	uint number = 0u; // number of neighboring interface points
	float2 p[6]; // number of neighboring interface points is less or equal than than 8 minus 1 gas and minus 1 fluid point = 6
	const float center_offset = plic_cube(phit[0], by); // calculate z-offset PLIC of center point only once
	for(uint i=1u; i<9u; i++) { // iterate over neighbors, no loop unrolling here (50% better perfoemance without loop unrolling)
		if(phit[i]>0.0f&&phit[i]<1.0f) { // limit neighbors to interface nodes
			const float3 ei = (float3)(c_D2Q9(i), c_D2Q9(9+i), 0.0f); // assume neighbor normal vector is the same as center normal vector
			const float offset = plic_cube(phit[i], by)-center_offset;
			p[number++] = (float2)(dot(ei, bx), dot(ei, by)+offset); // do coordinate system transformation into (x, f(x)) and apply PLIC pffsets
		}
	}
	float M[4]={0.0f,0.0f,0.0f,0.0f}, x[2]={0.0f,0.0f}, b[2]={0.0f,0.0f};
	for(uint i=0u; i<number; i++) { // f(x,y)=A*x2+H*x, x=(A,H), Q=(x2,x), M*x=b, M=Q*Q^T, b=Q*z
		const float x=p[i].x, y=p[i].y, x2=x*x, x3=x2*x;
		/**/M[0]+=x2*x2; M[1]+=x3; b[0]+=x2*y;
		/*M[2]+=x3   ;*/ M[3]+=x2; b[1]+=x *y;
	}
	M[2] = M[1]; // use symmetry of matrix to save arithmetic operations
	if(number>=2u) lu_solve(M, x, b, 2, 2);
	else lu_solve(M, x, b, 2, min(2u, number)); // cannot do loop unrolling here -> slower -> extra if-else to avoid slowdown
	const float A=x[0], H=x[1];
	const float K = 2.0f*A*cb(rsqrt(H*H+1.0f)); // mean curvature of Monge patch (x, f(x)), note that curvature definition in 2D is different than 3D (additional factor 2)
)+"#endif"+R( // D2Q9
	return clamp(K, -1.0f, 1.0f); // prevent extreme pressures in the case of almost degenerate matrices
}
)+"#endif"+R( // SURFACE

)+"#ifdef TEMPERATURE"+R(
)+R(void neighbors_temperature(const uint n, uint* j7) { // calculate neighbor indices
	uint x0, xp, xm, y0, yp, ym, z0, zp, zm;
	calculate_indices(n, &x0, &xp, &xm, &y0, &yp, &ym, &z0, &zp, &zm);
	j7[0] = n;
	j7[1] = xp+y0+z0; j7[2] = xm+y0+z0; // +00 -00
	j7[3] = x0+yp+z0; j7[4] = x0+ym+z0; // 0+0 0-0
	j7[5] = x0+y0+zp; j7[6] = x0+y0+zm; // 00+ 00-
}
)+R(void calculate_g_eq(const float T, const float ux, const float uy, const float uz, float* geq) { // calculate g_equilibrium from density and velocity field (perturbation method / DDF-shifting)
	const float wsT4=0.5f*T, wsTm1=0.125f*(T-1.0f); // 0.125f*T*4.0f (straight directions in D3Q7), wsTm1 is arithmetic optimization to minimize digit extinction, lattice speed of sound is 1/2 for D3Q7 and not 1/sqrt(3)
	geq[0] = fma(0.25f, T, -0.25f); // 000
	geq[1] = fma(wsT4, ux, wsTm1); geq[2] = fma(wsT4, -ux, wsTm1); // +00 -00, source: http://dx.doi.org/10.1016/j.ijheatmasstransfer.2009.11.014
	geq[3] = fma(wsT4, uy, wsTm1); geq[4] = fma(wsT4, -uy, wsTm1); // 0+0 0-0
	geq[5] = fma(wsT4, uz, wsTm1); geq[6] = fma(wsT4, -uz, wsTm1); // 00+ 00-
}
)+R(void load_g(const uint n, float* ghn, const global fpxx* gi, const uint* j7, const ulong t) {
	ghn[0] = load(gi, index_f(n, 0u)); // Esoteric-Pull
	for(uint i=1u; i<7u; i+=2u) {
		ghn[i   ] = load(gi, index_f(n    , t%2ul ? i    : i+1u));
		ghn[i+1u] = load(gi, index_f(j7[i], t%2ul ? i+1u : i   ));
	}
}
)+R(void store_g(const uint n, const float* ghn, global fpxx* gi, const uint* j7, const ulong t) {
	store(gi, index_f(n, 0u), ghn[0]); // Esoteric-Pull
	for(uint i=1u; i<7u; i+=2u) {
		store(gi, index_f(j7[i], t%2ul ? i+1u : i   ), ghn[i   ]);
		store(gi, index_f(n    , t%2ul ? i    : i+1u), ghn[i+1u]);
	}
}
)+"#endif"+R( // TEMPERATURE

)+R(void load_f(const uint n, float* fhn, const global fpxx* fi, const uint* j, const ulong t) {
	fhn[0] = load(fi, index_f(n, 0u)); // Esoteric-Pull
	for(uint i=1u; i<def_velocity_set; i+=2u) {
		fhn[i   ] = load(fi, index_f(n   , t%2ul ? i    : i+1u));
		fhn[i+1u] = load(fi, index_f(j[i], t%2ul ? i+1u : i   ));
	}
}
)+R(void store_f(const uint n, const float* fhn, global fpxx* fi, const uint* j, const ulong t) {
	store(fi, index_f(n, 0u), fhn[0]); // Esoteric-Pull
	for(uint i=1u; i<def_velocity_set; i+=2u) {
		store(fi, index_f(j[i], t%2ul ? i+1u : i   ), fhn[i   ]);
		store(fi, index_f(n   , t%2ul ? i    : i+1u), fhn[i+1u]);
	}
}

)+"#ifdef SURFACE"+R(
)+R(void load_f_outgoing(const uint n, float* fon, const global fpxx* fi, const uint* j, const ulong t) { // load outgoing DDFs, even: 1:1 like stream-out odd, odd: 1:1 like stream-out even
	for(uint i=1u; i<def_velocity_set; i+=2u) { // Esoteric-Pull
		fon[i   ] = load(fi, index_f(j[i], t%2ul ? i    : i+1u));
		fon[i+1u] = load(fi, index_f(n   , t%2ul ? i+1u : i   ));
	}
}
)+R(void store_f_reconstructed(const uint n, const float* fhn, global fpxx* fi, const uint* j, const ulong t, const uchar* flagsj_su) { // store reconstructed gas DDFs, even: 1:1 like stream-in even, odd: 1:1 like stream-in odd
	for(uint i=1u; i<def_velocity_set; i+=2u) { // Esoteric-Pull
		if(flagsj_su[i+1u]==TYPE_G) store(fi, index_f(n   , t%2ul ? i    : i+1u), fhn[i   ]); // only store reconstructed gas DDFs to locations from which
		if(flagsj_su[i   ]==TYPE_G) store(fi, index_f(j[i], t%2ul ? i+1u : i   ), fhn[i+1u]); // they are going to be streamed in during next stream_collide()
	}
}
)+"#endif"+R( // SURFACE



)+R(kernel void initialize)+"("+R(global fpxx* fi, const global float* rho, global float* u, global uchar* flags // ) { // initialize LBM
)+"#ifdef SURFACE"+R(
	, global float* mass, global float* massex, global float* phi // argument order is important
)+"#endif"+R( // SURFACE
)+"#ifdef TEMPERATURE"+R(
	, global fpxx* gi, const global float* T // argument order is important
)+"#endif"+R( // TEMPERATURE
)+") {"+R( // initialize()
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute initialize() on halo
	uchar flagsn = flags[n];
	const uchar flagsn_bo = flagsn&TYPE_BO; // extract boundary flags
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	uchar flagsj[def_velocity_set]; // cache neighbor flags for multiple readings
	for(uint i=1u; i<def_velocity_set; i++) flagsj[i] = flags[j[i]];
	if(flagsn_bo==TYPE_S) { // node is solid
		bool TYPE_ONLY_S = true; // has only solid neighbors
		for(uint i=1u; i<def_velocity_set; i++) TYPE_ONLY_S = TYPE_ONLY_S&&(flagsj[i]&TYPE_BO)==TYPE_S;
		if(TYPE_ONLY_S) {
			u[                 n] = 0.0f; // reset velocity for solid lattice points with only boundary neighbors
			u[    def_N+(ulong)n] = 0.0f;
			u[2ul*def_N+(ulong)n] = 0.0f;
		}
)+"#ifndef MOVING_BOUNDARIES"+R(
		if(flagsn_bo==TYPE_S) {
			u[                 n] = 0.0f; // reset velocity for all solid lattice points
			u[    def_N+(ulong)n] = 0.0f;
			u[2ul*def_N+(ulong)n] = 0.0f;
		}
)+"#else"+R( // MOVING_BOUNDARIES
	} else if(flagsn_bo!=TYPE_E) { // local lattice point is not solid and not equilibrium boundary
		bool next_to_moving_boundary = false;
		for(uint i=1u; i<def_velocity_set; i++) {
			next_to_moving_boundary = next_to_moving_boundary||((flagsj[i]&TYPE_BO)==TYPE_S&&(u[j[i]]!=0.0f||u[def_N+(ulong)j[i]]!=0.0f||u[2ul*def_N+(ulong)j[i]]!=0.0f));
		}
		flags[n] = flagsn = next_to_moving_boundary ? flagsn|TYPE_MS : flagsn&~TYPE_MS; // mark/unmark nodes next to TYPE_S nodes with velocity!=0 with TYPE_MS
)+"#endif"+R( // MOVING_BOUNDARIES
	}
	float feq[def_velocity_set]; // f_equilibrium
	calculate_f_eq(rho[n], u[n], u[def_N+(ulong)n], u[2ul*def_N+(ulong)n], feq);
)+"#ifdef SURFACE"+R( // automatically generate the interface layer between fluid and gas
	{ // separate block to avoid variable name conflicts
		float phin = phi[n];
		if(!(flagsn&(TYPE_S|TYPE_E|TYPE_T|TYPE_F|TYPE_I))) flagsn = (flagsn&~TYPE_SU)|TYPE_G; // change all non-fluid and non-interface flags to gas
		if((flagsn&TYPE_SU)==TYPE_G) { // node with updated flags is gas
			bool change = false; // check if node has to be changed to interface
			for(uint i=1u; i<def_velocity_set; i++) change = change||(flagsj[i]&TYPE_SU)==TYPE_F; // if neighbor flag fluid is set, the node must be interface
			if(change) { // create interface automatically if phi has not explicitely defined for the interface layer
				flagsn = (flagsn&~TYPE_SU)|TYPE_I; // node must be interface
				phin = 0.5f;
				float rhon, uxn, uyn, uzn; // initialize interface nodes with average density/velocity of fluid neighbors
				average_neighbors_fluid(n, rho, u, flags, &rhon, &uxn, &uyn, &uzn); // get average rho/u from all fluid neighbors
				calculate_f_eq(rhon, uxn, uyn, uzn, feq); // calculate equilibrium DDFs
			}
		}
		if((flagsn&TYPE_SU)==TYPE_G) { // node with updated flags is still gas
			u[                 n] = 0.0f; // reset velocity for gas nodes
			u[    def_N+(ulong)n] = 0.0f;
			u[2ul*def_N+(ulong)n] = 0.0f;
			phin = 0.0f;
		} else if((flagsn&TYPE_SU)==TYPE_I && (phin<0.0f||phin>1.0f)) {
			phin = 0.5f; // node should be interface, but phi was invalid
		} else if((flagsn&TYPE_SU)==TYPE_F) {
			phin = 1.0f;
		}
		phi[n] = phin;
		mass[n] = phin*rho[n];
		massex[n] = 0.0f; // reset excess mass
		flags[n] = flagsn;
	}
)+"#endif"+R( // SURFACE
)+"#ifdef TEMPERATURE"+R(
	{ // separate block to avoid variable name conflicts
		float geq[7];
		calculate_g_eq(T[n], u[n], u[def_N+(ulong)n], u[2ul*def_N+(ulong)n], geq);
		uint j7[7]; // neighbors of D3Q7 subset
		neighbors_temperature(n, j7);
		store_g(n, geq, gi, j7, 1ul);
	}
)+"#endif"+R( // TEMPERATURE
	store_f(n, feq, fi, j, 1ul); // write to fi
} // initialize()

)+"#ifdef MOVING_BOUNDARIES"+R(
)+R(kernel void update_moving_boundaries(const global float* u, global uchar* flags) { // mark/unmark nodes next to TYPE_S nodes with velocity!=0 with TYPE_MS
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute update_moving_boundaries() on halo
	const uchar flagsn = flags[n];
	const uchar flagsn_bo = flagsn&TYPE_BO; // extract boundary flags
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	uchar flagsj[def_velocity_set]; // cache neighbor flags for multiple readings
	for(uint i=1u; i<def_velocity_set; i++) flagsj[i] = flags[j[i]];
	if(flagsn_bo!=TYPE_S&&flagsn_bo!=TYPE_E&&!(flagsn&TYPE_T)) { // local lattice point is not solid and not equilibrium boundary and not temperature boundary
		bool next_to_moving_boundary = false;
		for(uint i=1u; i<def_velocity_set; i++) {
			next_to_moving_boundary = next_to_moving_boundary||((u[j[i]]!=0.0f||u[def_N+(ulong)j[i]]!=0.0f||u[2ul*def_N+(ulong)j[i]]!=0.0f)&&(flagsj[i]&TYPE_BO)==TYPE_S);
		}
		flags[n] = next_to_moving_boundary ? flagsn|TYPE_MS : flagsn&~TYPE_MS; // mark/unmark nodes next to TYPE_S nodes with velocity!=0 with TYPE_MS
	}
} // update_moving_boundaries()
)+"#endif"+R( // MOVING_BOUNDARIES



)+R(kernel void stream_collide)+"("+R(global fpxx* fi, global float* rho, global float* u, global uchar* flags, const ulong t, const float fx, const float fy, const float fz // ) { // main LBM kernel
)+"#ifdef FORCE_FIELD"+R(
	, const global float* F // argument order is important
)+"#endif"+R( // FORCE_FIELD
)+"#ifdef SURFACE"+R(
	, const global float* mass // argument order is important
)+"#endif"+R( // SURFACE
)+"#ifdef TEMPERATURE"+R(
	, global fpxx* gi, global float* T // argument order is important
)+"#endif"+R( // TEMPERATURE
)+") {"+R( // stream_collide()
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute stream_collide() on halo
	const uchar flagsn = flags[n]; // cache flags[n] for multiple readings
	const uchar flagsn_bo=flagsn&TYPE_BO, flagsn_su=flagsn&TYPE_SU; // extract boundary and surface flags
	if(flagsn_bo==TYPE_S||flagsn_su==TYPE_G) return; // if node is solid boundary or gas, just return

	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices

	float fhn[def_velocity_set]; // local DDFs
	load_f(n, fhn, fi, j, t); // perform streaming (part 2)

)+"#ifdef MOVING_BOUNDARIES"+R(
	if(flagsn_bo==TYPE_MS) apply_moving_boundaries(fhn, j, u, flags); // apply Dirichlet velocity boundaries if necessary (reads velocities of only neighboring boundary nodes, which do not change during simulation)
)+"#endif"+R( // MOVING_BOUNDARIES

	float rhon, uxn, uyn, uzn; // calculate local density and velocity for collision
)+"#ifndef EQUILIBRIUM_BOUNDARIES"+R(
	calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fi
)+"#else"+R( // EQUILIBRIUM_BOUNDARIES
	if(flagsn_bo==TYPE_E) {
		rhon = rho[               n]; // apply preset velocity/density
		uxn  = u[                 n];
		uyn  = u[    def_N+(ulong)n];
		uzn  = u[2ul*def_N+(ulong)n];
	} else {
		calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fi
	}
)+"#endif"+R( // EQUILIBRIUM_BOUNDARIES
	float fxn=fx, fyn=fy, fzn=fz; // force starts as constant volume force, can be modified before call of calculate_forcing_terms(...)
	float Fin[def_velocity_set]; // forcing terms

)+"#ifdef FORCE_FIELD"+R(
	{ // separate block to avoid variable name conflicts
		fxn += F[                 n]; // apply force field
		fyn += F[    def_N+(ulong)n];
		fzn += F[2ul*def_N+(ulong)n];
	}
)+"#endif"+R( // FORCE_FIELD

)+"#ifdef SURFACE"+R(
	if(flagsn_su==TYPE_I) { // node was interface, eventually initiate flag change
		bool TYPE_NO_F=true, TYPE_NO_G=true; // temporary flags for no fluid or gas neighbors
		for(uint i=1u; i<def_velocity_set; i++) {
			const uchar flagsji_su = flags[j[i]]&TYPE_SU; // extract SURFACE flags
			TYPE_NO_F = TYPE_NO_F&&flagsji_su!=TYPE_F;
			TYPE_NO_G = TYPE_NO_G&&flagsji_su!=TYPE_G;
		}
		const float massn = mass[n]; // load mass
		     if(massn>rhon || TYPE_NO_G) flags[n] = (flagsn&~TYPE_SU)|TYPE_IF; // set flag interface->fluid
		else if(massn<0.0f || TYPE_NO_F) flags[n] = (flagsn&~TYPE_SU)|TYPE_IG; // set flag interface->gas
	}
)+"#endif"+R( // SURFACE

)+"#ifdef TEMPERATURE"+R(
	{ // separate block to avoid variable name conflicts
		uint j7[7]; // neighbors of D3Q7 subset
		neighbors_temperature(n, j7);
		float ghn[7]; // read from gA and stream to gh (D3Q7 subset, periodic boundary conditions)
		load_g(n, ghn, gi, j7, t); // perform streaming (part 2)
		float Tn;
		if(flagsn&TYPE_T) {
			Tn = T[n]; // apply preset temperature
		} else {
			Tn = 0.0f;
			for(uint i=0u; i<7u; i++) Tn += ghn[i]; // calculate temperature from g
			Tn += 1.0f; // add 1.0f last to avoid digit extinction effects when summing up gi (perturbation method / DDF-shifting)
		}
		float geq[7]; // cache f_equilibrium[n]
		calculate_g_eq(Tn, uxn, uyn, uzn, geq); // calculate equilibrium DDFs
		if(flagsn&TYPE_T) {
			for(uint i=0u; i<7u; i++) ghn[i] = geq[i]; // just write geq to ghn (no collision)
		} else {
)+"#ifdef UPDATE_FIELDS"+R(
			T[n] = Tn; // update temperature field
)+"#endif"+R( // UPDATE_FIELDS
			for(uint i=0u; i<7u; i++) ghn[i] = fma(1.0f-def_w_T, ghn[i], def_w_T*geq[i]); // perform collision
		}
		store_g(n, ghn, gi, j7, t); // perform streaming (part 1)
		fxn -= fx*def_beta*(Tn-def_T_avg);
		fyn -= fy*def_beta*(Tn-def_T_avg);
		fzn -= fz*def_beta*(Tn-def_T_avg);
	}
)+"#endif"+R( // TEMPERATURE

	{ // separate block to avoid variable name conflicts
)+"#ifdef VOLUME_FORCE"+R( // apply force and collision operator, write to fi in video memory
		const float rho2 = 0.5f/rhon; // apply external volume force (Guo forcing, Krueger p.233f)
		uxn = clamp(fma(fxn, rho2, uxn), -def_c, def_c); // limit velocity (for stability purposes)
		uyn = clamp(fma(fyn, rho2, uyn), -def_c, def_c); // force term: F*dt/(2*rho)
		uzn = clamp(fma(fzn, rho2, uzn), -def_c, def_c);                                                                                             
		calculate_forcing_terms(uxn, uyn, uzn, fxn, fyn, fzn, Fin); // calculate volume force terms Fin from velocity field (Guo forcing, Krueger p.233f)
)+"#else"+R( // VOLUME_FORCE
		uxn = clamp(uxn, -def_c, def_c); // limit velocity (for stability purposes)
		uyn = clamp(uyn, -def_c, def_c); // force term: F*dt/(2*rho)
		uzn = clamp(uzn, -def_c, def_c);
		for(uint i=0u; i<def_velocity_set; i++) Fin[i] = 0.0f;
)+"#endif"+R( // VOLUME_FORCE
	}

)+"#ifndef EQUILIBRIUM_BOUNDARIES"+R(
)+"#ifdef UPDATE_FIELDS"+R(
	rho[               n] = rhon; // update density field
	u[                 n] = uxn; // update velocity field
	u[    def_N+(ulong)n] = uyn;
	u[2ul*def_N+(ulong)n] = uzn;
)+"#endif"+R( // UPDATE_FIELDS
)+"#else"+R( // EQUILIBRIUM_BOUNDARIES
)+"#ifdef UPDATE_FIELDS"+R(
	if(flagsn_bo!=TYPE_E) { // only update fields for non-TYPE_E nodes
		rho[               n] = rhon; // update density field
		u[                 n] = uxn; // update velocity field
		u[    def_N+(ulong)n] = uyn;
		u[2ul*def_N+(ulong)n] = uzn;
	}
)+"#endif"+R( // UPDATE_FIELDS
)+"#endif"+R( // EQUILIBRIUM_BOUNDARIES

	float feq[def_velocity_set]; // equilibrium DDFs
	calculate_f_eq(rhon, uxn, uyn, uzn, feq); // calculate equilibrium DDFs
	float w = def_w; // LBM relaxation rate w = dt/tau = dt/(nu/c^2+dt/2) = 1/(3*nu+1/2)

)+"#ifdef SUBGRID"+R(
	{ // Smagorinsky-Lilly subgrid turbulence model, source: https://arxiv.org/pdf/comp-gas/9401004.pdf, in the eq. below (26), it is "tau_0" not "nu_0", and "sqrt(2)/rho" (they call "rho" "n") is missing
		const float tau0 = 1.0f/w; // source 2: https://youtu.be/V8ydRrdCzl0
		float Hxx=0.0f, Hyy=0.0f, Hzz=0.0f, Hxy=0.0f, Hxz=0.0f, Hyz=0.0f; // non-equilibrium stress tensor
		for(uint i=1u; i<def_velocity_set; i++) {
			const float fneqi = fhn[i]-feq[i];
			const float cxi=c(i), cyi=c(def_velocity_set+i), czi=c(2u*def_velocity_set+i);
			Hxx += cxi*cxi*fneqi; //Hyx += cyi*cxi*fneqi; Hzx += czi*cxi*fneqi; // symmetric tensor
			Hxy += cxi*cyi*fneqi; Hyy += cyi*cyi*fneqi; //Hzy += czi*cyi*fneqi;
			Hxz += cxi*czi*fneqi; Hyz += cyi*czi*fneqi; Hzz += czi*czi*fneqi;
		}
		const float Q = sq(Hxx)+sq(Hyy)+sq(Hzz)+2.0f*(sq(Hxy)+sq(Hxz)+sq(Hyz)); // Q = H*H, turbulent eddy viscosity nut = (C*Delta)^2*|S|, intensity of local strain rate tensor |S|=sqrt(2*S*S)
		w = 2.0f/(tau0+sqrt(sq(tau0)+0.76421222f*sqrt(Q)/rhon)); // 0.76421222 = 18*sqrt(2)*(C*Delta)^2, C = 1/pi*(2/(3*CK))^(3/4) = Smagorinsky-Lilly constant, CK = 3/2 = Kolmogorov constant, Delta = 1 = lattice constant
	} // modity LBM relaxation rate by increasing effective viscosity in regions of high strain rate (add turbulent eddy viscosity), nu_eff = nu_0+nu_t
)+"#endif"+R( // SUBGRID

)+"#if defined(SRT)"+R(
)+"#ifdef VOLUME_FORCE"+R(
	const float c_tau = fma(w, -0.5f, 1.0f);
	for(uint i=0u; i<def_velocity_set; i++) Fin[i] *= c_tau;
)+"#endif"+R( // VOLUME_FORCE
)+"#ifndef EQUILIBRIUM_BOUNDARIES"+R(
	for(uint i=0u; i<def_velocity_set; i++) fhn[i] = fma(1.0f-w, fhn[i], fma(w, feq[i], Fin[i])); // perform collision (SRT)
)+"#else"+R( // EQUILIBRIUM_BOUNDARIES
	for(uint i=0u; i<def_velocity_set; i++) fhn[i] = flagsn_bo==TYPE_E ? feq[i] : fma(1.0f-w, fhn[i], fma(w, feq[i], Fin[i])); // perform collision (SRT)
)+"#endif"+R( // EQUILIBRIUM_BOUNDARIES

)+"#elif defined(TRT)"+R(
	const float wp = w; // TRT: inverse of "+" relaxation time
	const float wm = 1.0f/(0.25f/(1.0f/w-0.5f)+0.5f); // TRT: inverse of "-" relaxation time wm = 1.0f/(0.1875f/(3.0f*nu)+0.5f), nu = (1.0f/w-0.5f)/3.0f;
)+"#ifdef VOLUME_FORCE"+R(
	const float c_taup=fma(wp, -0.25f, 0.5f), c_taum=fma(wm, -0.25f, 0.5f); // source: https://arxiv.org/pdf/1901.08766.pdf
	float Fib[def_velocity_set]; // F_bar
	Fib[0] = Fin[0];
	for(uint i=1u; i<def_velocity_set; i+=2u) {
		Fib[i   ] = Fin[i+1u];
		Fib[i+1u] = Fin[i   ];
	}
	for(uint i=0u; i<def_velocity_set; i++) Fin[i] = fma(c_taup, Fin[i]+Fib[i], c_taum*(Fin[i]-Fib[i]));
)+"#endif"+R( // VOLUME_FORCE
	float fhb[def_velocity_set]; // fhn in inverse directions
	float feb[def_velocity_set]; // feq in inverse directions
	fhb[0] = fhn[0];
	feb[0] = feq[0];
	for(uint i=1u; i<def_velocity_set; i+=2u) {
		fhb[i   ] = fhn[i+1u];
		fhb[i+1u] = fhn[i   ];
		feb[i   ] = feq[i+1u];
		feb[i+1u] = feq[i   ];
	}
)+"#ifndef EQUILIBRIUM_BOUNDARIES"+R(
	for(uint i=0u; i<def_velocity_set; i++) fhn[i] = fma(0.5f*wp, feq[i]-fhn[i]+feb[i]-fhb[i], fma(0.5f*wm, feq[i]-feb[i]-fhn[i]+fhb[i], fhn[i]+Fin[i])); // perform collision (TRT)
)+"#else"+R( // EQUILIBRIUM_BOUNDARIES
	for(uint i=0u; i<def_velocity_set; i++) fhn[i] = flagsn_bo==TYPE_E ? feq[i] : fma(0.5f*wp, feq[i]-fhn[i]+feb[i]-fhb[i], fma(0.5f*wm, feq[i]-feb[i]-fhn[i]+fhb[i], fhn[i]+Fin[i])); // perform collision (TRT)
)+"#endif"+R( // EQUILIBRIUM_BOUNDARIES
)+"#endif"+R( // TRT

	store_f(n, fhn, fi, j, t); // perform streaming (part 1)
} // stream_collide()

)+"#ifdef SURFACE"+R(
)+R(kernel void surface_0(global fpxx* fi, const global float* rho, const global float* u, const global uchar* flags, global float* mass, const global float* massex, const global float* phi, const ulong t, const float fx, const float fy, const float fz) { // capture outgoing DDFs before streaming
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute surface_0() on halo
	const uchar flagsn = flags[n]; // cache flags[n] for multiple readings
	const uchar flagsn_bo=flagsn&TYPE_BO, flagsn_su=flagsn&TYPE_SU; // extract boundary and surface flags
	if(flagsn_bo==TYPE_S||flagsn_su==TYPE_G) return; // node processed here is fluid or interface

	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	float fhn[def_velocity_set]; // incoming DDFs
	load_f(n, fhn, fi, j, t); // load incoming DDFs
	float fon[def_velocity_set]; // outgoing DDFs
	fon[0] = fhn[0]; // fon[0] is already loaded in fhn[0]
	load_f_outgoing(n, fon, fi, j, t); // load outgoing DDFs

	float massn = mass[n];
	for(uint i=1u; i<def_velocity_set; i++) {
		massn += massex[j[i]]; // distribute excess mass from last step which is stored in neighbors
	}
	if(flagsn_su==TYPE_F) { // node is fluid
		for(uint i=1u; i<def_velocity_set; i++) massn += fhn[i]-fon[i]; // neighbor is fluid or interface node
	} else if(flagsn_su==TYPE_I) { // node is interface
		float phij[def_velocity_set]; // cache fill level of neighbor lattice points
		for(uint i=1u; i<def_velocity_set; i++) phij[i] = phi[j[i]]; // cache fill level of neighbor lattice points
		float rhon, uxn, uyn, uzn, rho_laplace=0.0f; // no surface tension if rho_laplace is not overwritten later
)+"#ifndef EQUILIBRIUM_BOUNDARIES"+R(
		calculate_rho_u(fon, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fon (not fhn)
)+"#else"+R( // EQUILIBRIUM_BOUNDARIES
		if(flagsn_bo==TYPE_E) {
			rhon = rho[               n]; // apply preset velocity/density
			uxn  = u[                 n];
			uyn  = u[    def_N+(ulong)n];
			uzn  = u[2ul*def_N+(ulong)n];
		} else {
			calculate_rho_u(fon, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fon (not fhn)
		}
)+"#endif"+R( // EQUILIBRIUM_BOUNDARIES
		uxn = clamp(uxn, -def_c, def_c); // limit velocity (for stability purposes)
		uyn = clamp(uyn, -def_c, def_c);
		uzn = clamp(uzn, -def_c, def_c);
		phij[0] = calculate_phi(rhon, massn, flagsn); // don't load phi[n] from memory, instead recalculate it with mass corrected by excess mass
		rho_laplace = def_6_sigma==0.0f ? 0.0f : def_6_sigma*curvature_calculation(n, phij, phi); // surface tension least squares fit (PLIC, most accurate)
		float feg[def_velocity_set]; // reconstruct f from neighbor gas lattice points
		const float rho2tmp = 0.5f/rhon; // apply external volume force (Guo forcing, Krueger p.233f)
		const float uxntmp = clamp(fma(fx, rho2tmp, uxn), -def_c, def_c); // limit velocity (for stability purposes)
		const float uyntmp = clamp(fma(fy, rho2tmp, uyn), -def_c, def_c); // force term: F*dt/(2*rho)
		const float uzntmp = clamp(fma(fz, rho2tmp, uzn), -def_c, def_c);
		calculate_f_eq(1.0f-rho_laplace, uxntmp, uyntmp, uzntmp, feg); // calculate gas equilibrium DDFs with constant ambient pressure
		uchar flagsj_su[def_velocity_set]; // cache neighbor flags for multiple readings
		for(uint i=1u; i<def_velocity_set; i++) flagsj_su[i] = flags[j[i]]&TYPE_SU;
		for(uint i=1u; i<def_velocity_set; i+=2u) { // calculate mass exchange between current node and fluid/interface nodes
			massn += flagsj_su[i   ]&(TYPE_F|TYPE_I) ? flagsj_su[i   ]==TYPE_F ? fhn[i+1]-fon[i   ] : 0.5f*(phij[i   ]+phij[0])*(fhn[i+1 ]-fon[i   ]) : 0.0f; // neighbor is fluid or interface node
			massn += flagsj_su[i+1u]&(TYPE_F|TYPE_I) ? flagsj_su[i+1u]==TYPE_F ? fhn[i  ]-fon[i+1u] : 0.5f*(phij[i+1u]+phij[0])*(fhn[i   ]-fon[i+1u]) : 0.0f; // fluid : interface : gas
		}
		for(uint i=1u; i<def_velocity_set; i+=2u) { // calculate reconstructed gas DDFs
			fhn[i   ] = feg[i+1u]-fon[i+1u]+feg[i   ];
			fhn[i+1u] = feg[i   ]-fon[i   ]+feg[i+1u];
		}
		store_f_reconstructed(n, fhn, fi, j, t, flagsj_su); // store reconstructed gas DDFs that are streamed in during the following stream_collide()
	}
	mass[n] = massn;
}
)+R(kernel void surface_1(global uchar* flags) { // prevent neighbors from interface->fluid nodes to become/be gas nodes
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N) return; // execute surface_1() also on halo
	const uchar flagsn_sus = flags[n]&(TYPE_SU|TYPE_S); // extract SURFACE flags
	if(flagsn_sus==TYPE_IF) { // flag interface->fluid is set
		uint j[def_velocity_set]; // neighbor indices
		neighbors(n, j); // calculate neighbor indices
		for(uint i=1u; i<def_velocity_set; i++) {
			const uchar flagsji = flags[j[i]];
			const uchar flagsji_su = flagsji&(TYPE_SU|TYPE_S); // extract SURFACE flags
			const uchar flagsji_r = flagsji&~TYPE_SU; // extract all non-SURFACE flags
			if(flagsji_su==TYPE_IG) flags[j[i]] = flagsji_r|TYPE_I; // prevent interface neighbor nodes from becoming gas
			else if(flagsji_su==TYPE_G) flags[j[i]] = flagsji_r|TYPE_GI; // neighbor node was gas and must change to interface
		}
	}
} // possible types at the end of surface_1(): TYPE_F / TYPE_I / TYPE_G / TYPE_IF / TYPE_IG / TYPE_GI
)+R(kernel void surface_2(global fpxx* fi, const global float* rho, const global float* u, global uchar* flags, const ulong t) {  // apply flag changes and calculate excess mass
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N) return; // execute surface_2() also on halo
	const uchar flagsn_sus = flags[n]&(TYPE_SU|TYPE_S); // extract SURFACE flags
	if(flagsn_sus==TYPE_GI) { // initialize the fi of gas nodes that should become interface
		float rhon, uxn, uyn, uzn; // average over all fluid/interface neighbors
		average_neighbors_non_gas(n, rho, u, flags, &rhon, &uxn, &uyn, &uzn); // get average rho/u from all fluid/interface neighbors
		float feq[def_velocity_set];
		calculate_f_eq(rhon, uxn, uyn, uzn, feq); // calculate equilibrium DDFs
		uint j[def_velocity_set];
		neighbors(n, j);
		store_f(n, feq, fi, j, t); // write feq to fi in video memory
	} else if(flagsn_sus==TYPE_IG) { // flag interface->gas is set
		uint j[def_velocity_set]; // neighbor indices
		neighbors(n, j); // calculate neighbor indices
		for(uint i=1u; i<def_velocity_set; i++) {
			const uchar flagsji = flags[j[i]];
			const uchar flagsji_su = flagsji&(TYPE_SU|TYPE_S); // extract SURFACE flags
			const uchar flagsji_r = flagsji&~TYPE_SU; // extract all non-SURFACE flags
			if(flagsji_su==TYPE_F||flagsji_su==TYPE_IF) {
				flags[j[i]] = flagsji_r|TYPE_I; // prevent fluid or interface neighbors that turn to fluid from being/becoming fluid
			}
		}
	}
} // possible types at the end of surface_2(): TYPE_F / TYPE_I / TYPE_G / TYPE_IF / TYPE_IG / TYPE_GI
)+R(kernel void surface_3(const global float* rho, global uchar* flags, global float* mass, global float* massex, global float* phi) { // apply flag changes and calculate excess mass
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute surface_3() on halo
	const uchar flagsn_sus = flags[n]&(TYPE_SU|TYPE_S); // extract SURFACE flags
	if(flagsn_sus&TYPE_S) return;
	const float rhon = rho[n]; // density of node n
	float massn = mass[n]; // mass of node n
	float massexn = 0.0f; // excess mass of node n
	float phin = 0.0f;
	if(flagsn_sus==TYPE_F) { // regular fluid node
		massexn = massn-rhon; // dump mass-rho difference into excess mass
		massn = rhon; // fluid node mass has to equal rho
		phin = 1.0f;
	} else if(flagsn_sus==TYPE_I) { // regular interface node
		massexn = massn>rhon ? massn-rhon : massn<0.0f ? massn : 0.0f; // allow interface nodes with mass>rho or mass<0
		massn = clamp(massn, 0.0f, rhon);
		phin = calculate_phi(rhon, massn, TYPE_I); // calculate fill level for next step (only necessary for interface nodes)
	} else if(flagsn_sus==TYPE_G) { // regular gas node
		massexn = massn; // dump remaining mass into excess mass
		massn = 0.0f;
		phin = 0.0f;
	} else if(flagsn_sus==TYPE_IF) { // flag interface->fluid is set
		flags[n] = (flags[n]&~TYPE_SU)|TYPE_F; // node becomes fluid
		massexn = massn-rhon; // dump mass-rho difference into excess mass
		massn = rhon; // fluid node mass has to equal rho
		phin = 1.0f; // set phi[n] to 1.0f for fluid nodes
	} else if(flagsn_sus==TYPE_IG) { // flag interface->gas is set
		flags[n] = (flags[n]&~TYPE_SU)|TYPE_G; // node becomes gas
		massexn = massn; // dump remaining mass into excess mass
		massn = 0.0f; // gas mass has to be zero
		phin = 0.0f; // set phi[n] to 0.0f for gas nodes
	} else if(flagsn_sus==TYPE_GI) { // flag gas->interface is set
		flags[n] = (flags[n]&~TYPE_SU)|TYPE_I; // node becomes interface
		massexn = massn>rhon ? massn-rhon : massn<0.0f ? massn : 0.0f; // allow interface nodes with mass>rho or mass<0
		massn = clamp(massn, 0.0f, rhon);
		phin = calculate_phi(rhon, massn, TYPE_I); // calculate fill level for next step (only necessary for interface nodes)
	}
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	uint counter = 0u; // count (fluid|interface) neighbors
	for(uint i=1u; i<def_velocity_set; i++) { // simple model: distribute excess mass equally to all interface and fluid neighbors
		const uchar flagsji_su = flags[j[i]]&(TYPE_SU|TYPE_S); // extract SURFACE flags
		counter += (uint)(flagsji_su==TYPE_F||flagsji_su==TYPE_I||flagsji_su==TYPE_IF||flagsji_su==TYPE_GI); // avoid branching
	}
	massn += counter>0u ? 0.0f : massexn; // if excess mass can't be distributed to neighboring interface or fluid nodes, add it to local mass (ensure mass conservation)
	massexn = counter>0u ? massexn/(float)counter : 0.0f; // divide excess mass up for all interface or fluid neighbors
	mass[n] = massn; // update mass
	massex[n] = massexn; // update excess mass
	phi[n] = phin; // update phi
} // possible types at the end of surface_3(): TYPE_F / TYPE_I / TYPE_G
)+"#endif"+R( // SURFACE

)+R(kernel void update_fields)+"("+R(const global fpxx* fi, global float* rho, global float* u, const global uchar* flags, const ulong t, const float fx, const float fy, const float fz // ) { // calculate fields from DDFs
)+"#ifdef FORCE_FIELD"+R(
	, const global float* F // argument order is important
)+"#endif"+R( // FORCE_FIELD
)+"#ifdef TEMPERATURE"+R(
	, const global fpxx* gi, global float* T // argument order is important
)+"#endif"+R( // TEMPERATURE
)+") {"+R( // update_fields()
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute update_fields() on halo
	const uchar flagsn = flags[n];
	const uchar flagsn_bo=flagsn&TYPE_BO, flagsn_su=flagsn&TYPE_SU; // extract boundary and surface flags
	if(flagsn_bo==TYPE_S||flagsn_su==TYPE_G) return; // don't update fields for boundary or gas lattice points

	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	float fhn[def_velocity_set]; // local DDFs
	load_f(n, fhn, fi, j, t); // perform streaming (part 2)

)+"#ifdef MOVING_BOUNDARIES"+R(
	if(flagsn_bo==TYPE_MS) apply_moving_boundaries(fhn, j, u, flags); // apply Dirichlet velocity boundaries if necessary (reads velocities of only neighboring boundary nodes, which do not change during simulation)
)+"#endif"+R( // MOVING_BOUNDARIES

	float rhon, uxn, uyn, uzn; // calculate local density and velocity for collision
	calculate_rho_u(fhn, &rhon, &uxn, &uyn, &uzn); // calculate density and velocity fields from fi
	float fxn=fx, fyn=fy, fzn=fz; // force starts as constant volume force, can be modified before call of calculate_forcing_terms(...)

)+"#ifdef FORCE_FIELD"+R(
	{ // separate block to avoid variable name conflicts
		fxn += F[                 n]; // apply force field
		fyn += F[    def_N+(ulong)n];
		fzn += F[2ul*def_N+(ulong)n];
	}
)+"#endif"+R( // FORCE_FIELD

)+"#ifdef TEMPERATURE"+R(
	{ // separate block to avoid variable name conflicts
		uint j7[7]; // neighbors of D3Q7 subset
		neighbors_temperature(n, j7);
		float ghn[7]; // read from gA and stream to gh (D3Q7 subset, periodic boundary conditions)
		load_g(n, ghn, gi, j7, t); // perform streaming (part 2)
		float Tn;
		if(flagsn&TYPE_T) {
			Tn = T[n]; // apply preset temperature
		} else {
			Tn = 0.0f;
			for(uint i=0u; i<7u; i++) Tn += ghn[i]; // calculate temperature from g
			Tn += 1.0f; // add 1.0f last to avoid digit extinction effects when summing up gi (perturbation method / DDF-shifting)
			T[n] = Tn; // update temperature field
		}
		fxn -= fx*def_beta*(Tn-def_T_avg);
		fyn -= fy*def_beta*(Tn-def_T_avg);
		fzn -= fz*def_beta*(Tn-def_T_avg);
	}
)+"#endif"+R( // TEMPERATURE

	{ // separate block to avoid variable name conflicts
)+"#ifdef VOLUME_FORCE"+R( // apply force and collision operator, write to fi in video memory
		const float rho2 = 0.5f/rhon; // apply external volume force (Guo forcing, Krueger p.233f)
		uxn = clamp(fma(fxn, rho2, uxn), -def_c, def_c); // limit velocity (for stability purposes)
		uyn = clamp(fma(fyn, rho2, uyn), -def_c, def_c); // force term: F*dt/(2*rho)
		uzn = clamp(fma(fzn, rho2, uzn), -def_c, def_c);
)+"#else"+R( // VOLUME_FORCE
		uxn = clamp(uxn, -def_c, def_c); // limit velocity (for stability purposes)
		uyn = clamp(uyn, -def_c, def_c); // force term: F*dt/(2*rho)
		uzn = clamp(uzn, -def_c, def_c);
)+"#endif"+R( // VOLUME_FORCE
	}

)+"#ifndef EQUILIBRIUM_BOUNDARIES"+R(
	rho[               n] = rhon; // update density field
	u[                 n] = uxn; // update velocity field
	u[    def_N+(ulong)n] = uyn;
	u[2ul*def_N+(ulong)n] = uzn;
)+"#else"+R( // EQUILIBRIUM_BOUNDARIES
	if(flagsn_bo!=TYPE_E) { // only update fields for non-TYPE_E nodes
		rho[               n] = rhon; // update density field
		u[                 n] = uxn; // update velocity field
		u[    def_N+(ulong)n] = uyn;
		u[2ul*def_N+(ulong)n] = uzn;
	}
)+"#endif"+R( // EQUILIBRIUM_BOUNDARIES
} // update_fields()

)+"#ifdef FORCE_FIELD"+R(
)+R(kernel void calculate_force_on_boundaries(const global fpxx* fi, const global uchar* flags, const ulong t, global float* F) { // calculate force from the fluid on solid boundaries from fi directly
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N||is_halo(n)) return; // don't execute calculate_force_on_boundaries() on halo
	if((flags[n]&TYPE_BO)!=TYPE_S) return; // only continue for solid boundary nodes
	uint j[def_velocity_set]; // neighbor indices
	neighbors(n, j); // calculate neighbor indices
	float fhn[def_velocity_set]; // local DDFs
	load_f(n, fhn, fi, j, t); // perform streaming (part 2)
	float Fb=1.0f, fx=0.0f, fy=0.0f, fz=0.0f;
	calculate_rho_u(fhn, &Fb, &fx, &fy, &fz); // abuse calculate_rho_u() method for calculating force
	F[                 n] = 2.0f*fx*Fb; // 2 times because fi are reflected on solid boundary nodes (bounced-back)
	F[    def_N+(ulong)n] = 2.0f*fy*Fb;
	F[2ul*def_N+(ulong)n] = 2.0f*fz*Fb;
} // calculate_force_on_boundaries()
)+R(kernel void reset_force_field(global float* F) { // reset force field
	const uint n = get_global_id(0); // n = x+(y+z*Ny)*Nx
	if(n>=(uint)def_N) return; // execute reset_force_field() also on halo
	F[                 n] = 0.0f;
	F[    def_N+(ulong)n] = 0.0f;
	F[2ul*def_N+(ulong)n] = 0.0f;
} // reset_force_field()
)+R(void spread_force(volatile global float* F, const float3 p, const float3 Fn) {
	const float xa=p.x-0.5f+1.5f*def_Nx, ya=p.y-0.5f+1.5f*def_Ny, za=p.z-0.5f+1.5f*def_Nz; // subtract lattice offsets
	const uint xb=(uint)xa, yb=(uint)ya, zb=(uint)za; // integer casting to find bottom left corner
	const float x1=xa-(float)xb, y1=ya-(float)yb, z1=za-(float)zb; // calculate interpolation factors
	for(uint c=0u; c<8u; c++) { // count over eight corner points
		const uint i=(c&0x04u)>>2, j=(c&0x02u)>>1, k=c&0x01u; // disassemble c into corner indices ijk
		const uint x=(xb+i)%def_Nx, y=(yb+j)%def_Ny, z=(zb+k)%def_Nz; // calculate corner lattice positions
		const uint n = x+(y+z*def_Ny)*def_Nx; // calculate lattice linear index
		const float d = (1.0f-fabs(x1-(float)i))*(1.0f-fabs(y1-(float)j))*(1.0f-fabs(z1-(float)k)); // force spreading
		atomic_add_f(&F[                 n], Fn.x*d); // F[                 n] += Fn.x*d;
		atomic_add_f(&F[    def_N+(ulong)n], Fn.y*d); // F[    def_N+(ulong)n] += Fn.y*d;
		atomic_add_f(&F[2ul*def_N+(ulong)n], Fn.z*d); // F[2ul*def_N+(ulong)n] += Fn.z*d;
	}
} // spread_force()
)+"#endif"+R( // FORCE_FIELD

)+"#ifdef PARTICLES"+R(
)+R(kernel void integrate_particles)+"("+R(global float* particles, const global float* u // ) {
)+"#ifdef FORCE_FIELD"+R(
	, volatile global float* F, const float fx, const float fy, const float fz
)+"#endif"+R( // FORCE_FIELD
)+") {"+R( // integrate_particles()
	const uint n = get_global_id(0); // index of membrane points
	if(n>=(uint)def_particles_N) return;
	const float3 p0 = (float3)(particles[n], particles[def_particles_N+(ulong)n], particles[2ul*def_particles_N+(ulong)n]); // cache particle position
)+"#ifdef FORCE_FIELD"+R(
	if(def_particles_rho!=1.0f) {
		const float drho = def_particles_rho-1.0f; // density difference leads to particle buoyancy
		float3 Fn = (float3)(fx*drho, fy*drho, fz*drho); // F = F_p+F_f = (m_p-m_f)*g = (rho_p-rho_f)*g*V
		spread_force(F, p0, Fn); // do force spreading
	}
)+"#endif"+R( // FORCE_FIELD
	const float3 un = interpolate_u(mirror_position(p0), u); // trilinear interpolation of velocity at point p
	const float3 p = mirror_position(p0+un); // advect particles
	particles[                           n] = p.x;
	particles[    def_particles_N+(ulong)n] = p.y;
	particles[2ul*def_particles_N+(ulong)n] = p.z;
} // integrate_particles()
)+"#endif"+R( // PARTICLES


);} // ############################################################### end of OpenCL C code #####################################################################