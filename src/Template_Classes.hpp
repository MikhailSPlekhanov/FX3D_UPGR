#pragma once

#include "utilities.hpp"
#include "defines.hpp"
#include "opencl.hpp"
#include "graphics.hpp"
#include "units.hpp"
#include "info.hpp"


class DomainBase {
private:
	uint Nx = 1u, Ny = 1u, Nz = 1u; // (local) lattice dimensions
	uint Dx = 1u, Dy = 1u, Dz = 1u; // lattice domains
	int Ox = 0, Oy = 0, Oz = 0; // lattice domain offset
	ulong t = 0ull; // discrete time step in LBM units

	float nu = 1.0f / 6.0f; // kinematic shear viscosity
	float fx = 0.0f, fy = 0.0f, fz = 0.0f; // global force per volume
	float sigma = 0.0f; // surface tension coefficient
	float alpha = 1.0f, beta = 1.0f, T_avg = 1.0f; // alpha = thermal diffusion coefficient, beta = thermal expansion coefficient, T_avg = 1 = average temperature
	uint particles_N = 0u;
	float particles_rho = 1.0f;

	Device device; // OpenCL device associated with this LBM domain
	Kernel kernel_initialize; // initialization kernel
	Kernel kernel_stream_collide; // main LBM kernel
	Kernel kernel_update_fields; // reads DDFs and updates (rho, u, T) in device memory
	Memory<fpxx> fi; // LBM density distribution functions (DDFs); only exist in device memory
	ulong t_last_update_fields = 0ull; // optimization to not call kernel_update_fields multiple times if (rho, u, T) are already up-to-date

	void allocate(Device& device); // allocate all memory for data fields on host and device and set up kernels
	string device_defines() const; // returns preprocessor constants for embedding in OpenCL C code
public:
	Memory<float> rho; // density of every node
	Memory<float> u; // velocity of every node
	Memory<uchar> flags; // flags of every node

	Memory<char> transfer_buffer_p, transfer_buffer_m; // transfer buffers for multi-device domain communication, only allocate one set of transfer buffers in plus/minus directions, for all x/y/z transfers
	Kernel kernel_transfer[enum_transfer_field::enum_transfer_field_length][2]; // for each field one extract and one insert kernel
	void allocate_transfer(Device& device); // allocate all memory for multi-device transfer
	ulong get_area(const uint direction);
	void enqueue_transfer_extract_field(Kernel& kernel_transfer_extract_field, const uint direction, const uint bytes_per_cell);
	void enqueue_transfer_insert_field(Kernel& kernel_transfer_insert_field, const uint direction, const uint bytes_per_cell);

	DomainBase(const Device_Info& device_info, const uint Nx, const uint Ny, const uint Nz, const uint Dx, const uint Dy, const uint Dz, const int Ox, const int Oy, const int Oz, const float nu, const float fx, const float fy, const float fz, const float sigma, const float alpha, const float beta, const uint particles_N, const float particles_rho); // compiles OpenCL C code and allocates memory

	void enqueue_initialize(); // write all data fields to device and call kernel_initialize
	void enqueue_stream_collide(); // call kernel_stream_collide to perform one LBM time step
	void enqueue_update_fields(); // update fields (rho, u, T) manually
	void increment_time_step(); // increment time step
	void reset_time_step(); // reset time step
	void finish_queue();

	const Device& get_device() const { return device; }
	uint get_Nx() const { return Nx; } // get (local) lattice dimensions in x-direction
	uint get_Ny() const { return Ny; } // get (local) lattice dimensions in y-direction
	uint get_Nz() const { return Nz; } // get (local) lattice dimensions in z-direction
	ulong get_N() const { return (ulong)Nx * (ulong)Ny * (ulong)Nz; } // get (local) number of lattice points
	uint3 get_Dim() const { return uint3(Nx, Ny, Nz); }
	uint get_Dx() const { return Dx; } // get lattice domains in x-direction
	uint get_Dy() const { return Dy; } // get lattice domains in y-direction
	uint get_Dz() const { return Dz; } // get lattice domains in z-direction
	uint get_D() const { return Dx * Dy * Dz; } // get number of lattice domains
	float get_nu() const { return nu; } // get kinematic shear viscosity
	float get_tau() const { return 3.0f * get_nu() + 0.5f; } // get LBM relaxation time
	float get_fx() const { return fx; } // get global froce per volume
	float get_fy() const { return fy; } // get global froce per volume
	float get_fz() const { return fz; } // get global froce per volume
	float get_sigma() const { return sigma; } // get surface tension coefficient
	float get_alpha() const { return alpha; } // get thermal diffusion coefficient
	float get_beta() const { return beta; } // get thermal expansion coefficient
	ulong get_t() const { return t; } // get discrete time step in LBM units
	uint get_velocity_set() const; // get LBM velocity set
	void set_fx(const float fx) { this->fx = fx; } // set global froce per volume
	void set_fy(const float fy) { this->fy = fy; } // set global froce per volume
	void set_fz(const float fz) { this->fz = fz; } // set global froce per volume
	void set_f(const float fx, const float fy, const float fz) { set_fx(fx); set_fy(fy); set_fz(fz); } // set global froce per volume
	void set_nu(const float nu) { this->nu = nu; }

	void voxelize_mesh_on_device(const Mesh* mesh, const uchar flag = TYPE_S, const float3& rotation_center = float3(0.0f), const float3& linear_velocity = float3(0.0f), const float3& rotational_velocity = float3(0.0f)); // voxelize mesh
	void enqueue_unvoxelize_mesh_on_device(const Mesh* mesh, const uchar flag = TYPE_S); // remove voxelized triangle mesh from LBM grid
};

class DomainForceFieldUnit {
private:
	Kernel kernel_calculate_force_on_boundaries; // calculate forces from fluid on TYPE_S nodes
	Kernel kernel_reset_force_field; // reset force field (also on TYPE_S nodes)
public:
	Memory<float> F; // individual force for every node
	void enqueue_calculate_force_on_boundaries(); // calculate forces from fluid on TYPE_S nodes
};

class DomainSurfaceUnit {
private:
	Kernel kernel_surface_0; // additional kernel for computing mass conservation and mass flux computation
	Kernel kernel_surface_1; // additional kernel for flag handling
	Kernel kernel_surface_2; // additional kernel for flag handling
	Kernel kernel_surface_3; // additional kernel for flag handling and mass conservation
	Memory<float> mass; // fluid mass; phi=mass/rho
	Memory<float> massex; // excess mass; used for mass conservation
public:
	Memory<float> phi; // fill level of every node
	void enqueue_surface_0();
	void enqueue_surface_1();
	void enqueue_surface_2();
	void enqueue_surface_3();
};

class DomainTemperatureUnit {
private:
	Memory<fpxx> gi; // thermal DDFs
public:
	Memory<float> T; // temperature of every node
};

class DomainParticlesUnit {
private:
	Kernel kernel_integrate_particles; // intgegrates particles forward in time and couples particles to fluid
public:
	Memory<float> particles; // particle positions
	void enqueue_integrate_particles(); // intgegrates particles forward in time and couples particles to fluid
};

class DomainGraphicsBase {
private:
	Kernel kernel_clear; // reset bitmap and zbuffer
	Memory<int> bitmap; // bitmap for rendering
	Memory<int> zbuffer; // z-buffer for rendering
	Memory<float> camera_parameters; // contains camera position, rotation, field of view etc.

	DomainBase* lbm = nullptr;
	Kernel kernel_graphics_flags; // render flag lattice
	Kernel kernel_graphics_field; // render a colored velocity vector for each node
	Kernel kernel_graphics_streamline; // render streamlines
	Kernel kernel_graphics_q; // render vorticity (Q-criterion)

	ulong t_last_frame = 0ull; // optimization to not call draw_frame() multiple times if camera_parameters and LBM time step are unchanged
	bool update_camera(); // update camera_parameters and return true if they are changed from their previous state
public:
	DomainGraphicsBase() {} // default constructor
	DomainGraphicsBase(DomainBase* lbm) {
		this->lbm = lbm;
		if constexpr (SURFACE) { ///
			skybox_image = read_png(path_skybox); //!!!!!!!!!! image существует в блоке surface
		}
	}
	DomainGraphicsBase& operator=(const DomainGraphicsBase& graphics) { // copy assignment
		lbm = graphics.lbm;
		if constexpr (SURFACE) {
			skybox_image = graphics.get_skybox_image();
		}
		return *this;
	}
	void allocate(Device& device); // allocate memory for bitmap and zbuffer
	void enqueue_draw_frame(); // main rendering function, calls rendering kernels
	int* get_bitmap(); // returns pointer to bitmap
	int* get_zbuffer(); // returns pointer to zbuffer
	string device_defines() const; // returns preprocessor constants for embedding in OpenCL C code
};

class DomainGraphicsSurface {
private:

public:

};

class DomainGraphicsParticles {
private:
	Kernel kernel_graphics_particles;
};

class LBMBase {
private:
	uint Nx = 1u, Ny = 1u, Nz = 1u; // (global) lattice dimensions
	uint Dx = 1u, Dy = 1u, Dz = 1u; // lattice domains
	bool initialized = false; // becomes true after LBM::initialize() has been called

	void sanity_checks_constructor(const vector<Device_Info>& device_infos, const uint Nx, const uint Ny, const uint Nz, const uint Dx, const uint Dy, const uint Dz, const float nu, const float fx, const float fy, const float fz, const float sigma, const float alpha, const float beta, const uint particles_N, const float particles_rho); // sanity checks on grid resolution and extension support
	void sanity_checks_initialization(); // sanity checks during initialization on used extensions based on used flags
	void initialize(); // write all data fields to device and call kernel_initialize
	void do_time_step(); // call kernel_stream_collide to perform one LBM time step

	void communicate_field(const enum_transfer_field field, const uint bytes_per_cell);

	void communicate_fi();
	void communicate_rho_u_flags();

public:
	/* MemoryContainer - no problem */
};


template<typename T,  typename U, typename... Units> //  
class LBM_variadic : T { // если наследовать, то потом непонятно как стучаться
	T member;
	// тут всякие функции-члены T, которые используют U в качестве аргумента? а если используют не только U?? 
	// а надо ли вообще так делать, если использование через constexpr определится?
	// а как продолжать рекурсию??? желательно чтобы все Units были членами LBM_variadic
	
};

// или лучше так?
template<typename... Units>
class LBM_includeall : Units ...{

};


// возможный сценарий использования:
// вызываем конструктор типа LBM_variadic<SurfaceUnit, TemperatureUnit> (constexpr LBMParameterList list);
// 
//