#ifndef _smoke_solver3D_
#define _smoke_solver3D_
#include "array.h"
#include "tbb/tbb.h"
//#include "Multigrid3D.h"
#include "fluid_buffer3D.h"
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstdio>
#include "fluid_particle.h"
#include "vec.h"
#include "pcg_solver.h"
#include "array3.h"

#include "GeometricLevelGen.h"

#include "AlgebraicMultigrid.h"
//using namespace gf;
using namespace std;
class SmokeSolver3D
{
public:
	SmokeSolver3D(){}
	~SmokeSolver3D(){Finalize();}
	Vec3f emitter_pos;
	float emitter_r;
	uint  emitter_n;
	vector<Vec4f> tracers;
	float frand(float a, float b)
	{
		return a + (b-a)*((float)(rand()%RAND_MAX)/(float)RAND_MAX);
	}
	void setEmitter(const Vec3f & pos, const float &r, uint n)
	{
		emitter_pos = pos;
		emitter_r = r;
		emitter_n = n;
	}
	void emit_tracers()
	{
		//cout<<"emitting tracers:"<<endl;
		vector<Vec4f> tracers_temp;
		tracers_temp.resize(0);
		for (uint i=0;i<tracers.size();i++)
		{
			if (tracers[i][0]>2*_hx && 
				tracers[i][1]>2*_hx &&
				tracers[i][2]>2*_hx &&
				tracers[i][0]<_lx-2*_hx &&
				tracers[i][1]<0.86*_ly-_hx &&
				tracers[i][2]<_lz-2*_hx &&
				tracers[i][3>0.01])
			{
				tracers_temp.push_back(tracers[i]);
			}
		}
		tracers.swap(tracers_temp);
		uint num = 0;
		while(num<emitter_n)
		{
			float r = emitter_r;
			float x = frand(-r-_hx, r+_hx);
			float y = frand(-r-_hx, r+_hx);
			float z = frand(-r-_hx, r+_hx);

			if (x*x + y*y + z*z <= r*r)
			{
				tracers.push_back(Vec4f(emitter_pos[0]+x,
									    emitter_pos[1]+y,
										emitter_pos[2]+z,
										1.0)
								 );
				num++;
			}
		}
		//cout<<"emitting tracers done:"<<tracers.size()<<" tracers"<<endl;
	}
	void advect_tracers(float dt)
	{
		tbb::parallel_for((size_t)0,
			              (size_t)tracers.size(),
						  (size_t)1,
						  [&](size_t i)
		{
			Vec3f pos = Vec3f(tracers[i][0],tracers[i][1],tracers[i][2]);
			pos = traceRK3(pos,dt);
			tracers[i] = Vec4f(pos[0],pos[1],pos[2], tracers[i][3]/(1.0+0.1*dt));
		});
	}
	void write_tracers(char * file_path, int frame)
	{
		char file_name[256];
		sprintf(file_name,"%s/Particle_data%04d.bin", file_path,frame);
		float* data;
		data = new float[4*tracers.size()];
		
		tbb::parallel_for((size_t)0,
			(size_t)tracers.size(),
			(size_t)1,
			[&](size_t i)
		{
			data[i*4+0] = tracers[i][0]/_lx - 0.5;
			data[i*4+1] = tracers[i][1]/_lx;
			data[i*4+2] = tracers[i][2]/_lx - 0.5;
			data[i*4+3] = tracers[i][3];
		});

		FILE *data_file = fopen(file_name,"wb");
		fwrite(data,sizeof(float)*4,tracers.size(),data_file);
		fclose(data_file);
		delete []data;
	}
	uint _nx, _ny, _nz, _n;
	float _hx, _hy, _hz;
	float _lx, _ly, _lz;
	float _temp_decay;
	float _alpha, _beta;
	float _smoke_heat, _smoke_dens, _smoke_fuel;
	//buffers:
	buffer3Df _un, _vn, _wn, _utemp,_vtemp,_wtemp,_unp1,_vnp1,_wnp1;
	Array3f u_extrap,v_extrap,w_extrap;
	Array3c u_valid, v_valid, w_valid;
	buffer3Df _wxn, _wyn, _wzn, _wxnp1,_wynp1,_wznp1;
	buffer3Df _wxstr, _wystr, _wzstr;
	buffer3Df _Psix,_Psiy,_Psiz;
	buffer3Df _Tbf, _rho, _fuel, _Ttemp,_Ttempnp1;
	buffer3Dc _b_desc, _h_desc, _boundary_marker;
	//MGSolverf _ppe_solver;
	buffer3Df _p;
	buffer3Df _div, _burn_div;
	float _cfl;
	float _vort_confine_coef;
	//fluid_particle
	fluid_particle FLIP_Solver;
	
	buffer3Df vort_confine_x,vort_confine_y,vort_confine_z;


	//solver for the vortex part
	levelGen<double> amg_levelGen;
	FixedSparseMatrix<double> fixed_matrix;
	vector<FixedSparseMatrix<double> *> A_L;
	vector<FixedSparseMatrix<double> *> R_L;
	vector<FixedSparseMatrix<double> *> P_L;
	vector<Vec3i>                  S_L;
	int total_level;


	//Solver data
	PCGSolver<double> solver;
	SparseMatrixd matrix;
	std::vector<double> rhs;
	std::vector<double> pressure;

	void mark_boundary();
	void clearBoundary();
	Vec3f get_velocity(Vec3f & pos);
	Vec3f traceRK3(Vec3f & pos, float dt);
	Vec3f trace(float dt, Vec3f & pos);
	void getCFL(float dt );
	void extrapolate(Array3f & grid, buffer3Df & u, Array3c & valid);
	void set_vort_confine(float str) { _vort_confine_coef = str; }


	void Finalize()
	{
		_un.free();
		_utemp.free();
		_unp1.free();
		_vn.free();
		_vtemp.free();
		_vnp1.free();
		_wn.free();
		_wtemp.free();
		_wnp1.free();

		_Tbf.free();
		_rho.free();
		_fuel.free();
		_Ttemp.free();
		_Ttempnp1.free();

		_b_desc.free();
		_h_desc.free();
		//_ppe_solver.m_FinalMemoryEachLevel();
		_p.free();
		_div.free();
		_burn_div.free();

		_wxn.free();
		_wxnp1.free();
		_wyn.free();
		_wynp1.free();
		_wzn.free();
		_wznp1.free();
		_wxstr.free();
		_wystr.free();
		_wzstr.free();

		_Psix.free();
		_Psiy.free();
		_Psiz.free();
		vort_confine_x.free();
		vort_confine_y.free();
		vort_confine_z.free();
		_boundary_marker.free();
		                          
	}
	//init
	void init(uint nx, uint ny, uint nz, float L)
	{
		_nx=nx;
		_ny=ny;
		_nz=nz;
		_hx = L/(float)_nx;
		_hz = _hy = _hx;
		_lx = L; _ly = _hy*(float)ny; _lz = _hz*(float)nz;
		_temp_decay = 0;
		_alpha = _beta = 0;
		_smoke_dens = _smoke_heat = 0;

		FLIP_Solver.init(_nx,_ny,_nz,_hx);

		_un.init(_nx+1,_ny,_nz,_hx,0.5,0,0);
		_utemp.init(_nx+1,_ny,_nz,_hx,0.5,0,0);
		_unp1.init(_nx+1,_ny,_nz,_hx,0.5,0,0);
		_vn.init(_nx,_ny+1,_nz,_hx,0.0,0.5,0.0);
		_vtemp.init(_nx,_ny+1,_nz,_hx,0.0,0.5,0.0);
		_vnp1.init(_nx,_ny+1,_nz,_hx,0.0,0.5,0.0);
		_wn.init(_nx,_ny,_nz+1,_hx,0.0,0.0,0.5);
		_wtemp.init(_nx,_ny,_nz+1,_hx,0.0,0.0,0.5);
		_wnp1.init(_nx,_ny,_nz+1,_hx,0.0,0.0,0.5);
		u_valid.resize(_nx+1,_ny,_nz);
		v_valid.resize(_nx,_ny+1,_nz);
		w_valid.resize(_nx,_ny,_nz+1);
		u_extrap.resize(_nx+1,_ny,_nz);
		v_extrap.resize(_nx,_ny+1,_nz);
		w_extrap.resize(_nx,_ny,_nz+1);

		_Tbf.init(_nx,_ny,_nz,_hx,0,0,0);
		_rho.init(_nx,_ny,_nz,_hx,0,0,0);
		_fuel.init(_nx,_ny,_nz,_hx,0,0,0);
		_Ttemp.init(_nx,_ny,_nz,_hx,0,0,0);
		_Ttempnp1.init(_nx,_ny,_nz,_hx,0,0,0);
		vort_confine_x.init(_nx,_ny,_nz,_hx,0,0,0);
		vort_confine_y.init(_nx,_ny,_nz,_hx,0,0,0);
		vort_confine_z.init(_nx,_ny,_nz,_hx,0,0,0);

		_b_desc.init(_nx,_ny,_nz);
		_boundary_marker.init(_nx,_ny,_nz);
		_h_desc.init(_nx,_ny,_nz);
		//_ppe_solver.m_InitialSystem(_nx,_ny,_nz);
		_p.init(_nx,_ny,_nz);
		_div.init(_nx,_ny,_nz);
		_burn_div.init(_nx,_ny,_nz,_hx,0,0,0);


		_wxn.init(_nx,_ny,_nz,_hx,0,0.5,0.5);
		_wyn.init(_nx,_ny,_nz,_hx,0.5,0,0.5);
		_wzn.init(_nx,_ny,_nz,_hx,0.5,0.5,0);
		_wxnp1.init(_nx,_ny,_nz,_hx,0,0.5,0.5);
		_wynp1.init(_nx,_ny,_nz,_hx,0.5,0,0.5);
		_wznp1.init(_nx,_ny,_nz,_hx,0.5,0.5,0);
		_Psix.init(_nx,_ny,_nz);
		_Psiy.init(_nx,_ny,_nz);
		_Psiz.init(_nx,_ny,_nz);

		_wxstr.init(_nx,_ny,_nz,_hx,0,0.5,0.5);
		_wystr.init(_nx,_ny,_nz,_hx,0.5,0,0.5);
		_wzstr.init(_nx,_ny,_nz,_hx,0.5,0.5,0);



		matrix.resize(_nx*_ny*_nz);
		matrix.zero();
		int ni = _nx, nj = _ny, nk = _nz;
		int compute_num = ni*nj*nk;
		int slice = ni*nj;
		tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/ni;
			int i = thread_idx%ni;
			if(i>=0 && i<=ni-1 && j>=0 && j<=nj-1 && k>=0 && k<=nk-1)
			{
				int index = i + ni*j + ni*nj*k;

				

				if( _b_desc(i,j,k)==0 )//a fluid cell 
				{

					//right neighbour
					if( i+1<ni ) {//a fluid cell
						matrix.add_to_element(index, index, 1.0);
						matrix.add_to_element(index, index + 1, -1.0);
					}
					else 
					{
						matrix.add_to_element(index, index, 1.0);
					}
					

					//left neighbour
					if( i-1>=0 ) {
						matrix.add_to_element(index, index, 1.0);
						matrix.add_to_element(index, index - 1, -1.0);
					}
					else
					{
						matrix.add_to_element(index, index, 1.0);
					}



					//top neighbour
					if( j+1<nj ) {//a fluid cell
						matrix.add_to_element(index, index, 1.0);
						matrix.add_to_element(index, index + ni, -1.0);
					}
					else 
					{
						matrix.add_to_element(index, index, 1.0);
					}

					//bottom neighbour
					if( j-1>=0 ) {
						matrix.add_to_element(index, index, 1.0);
						matrix.add_to_element(index, index - ni, -1.0);
					}
					else 
					{
						matrix.add_to_element(index, index, 1.0);
					}

					//back neighbour
					if( k+1<nk ) {//a fluid cell
						matrix.add_to_element(index, index, 1.0);
						matrix.add_to_element(index, index + ni*nj, -1.0);
					}
					else
					{
						matrix.add_to_element(index, index, 1.0);
					}

					//front neighbour
					if( k-1>=0 ) {
						matrix.add_to_element(index, index, 1.0);
						matrix.add_to_element(index, index - ni*nj, -1.0);
					}
					else
					{

						matrix.add_to_element(index, index, 1.0);
					}
				}
			}
		});
		fixed_matrix.construct_from_matrix(matrix);
		matrix.zero();
		amg_levelGen.generateLevelsGalerkinCoarsening(A_L, R_L, P_L,S_L,total_level,fixed_matrix,ni,nj,nk);
	}
	void solve_stream_Poisson(buffer3Df &psi, buffer3Df &rhs)
	{
		vector<double> x;
		x.resize(_nx*_ny*_nz);
		x.assign(_nx*_ny*_nz, 0);
		vector<double> b;
		b.resize(_nx*_ny*_nz);
		b.assign(_nx*_ny*_nz, 0);


		int ni = _nx, nj = _ny, nk = _nz;
		int compute_num = ni*nj*nk;
		int slice = ni*nj;
		tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/ni;
			int i = thread_idx%ni;
			if(i>=0 && i<=ni-1 && j>=0 && j<=nj-1 && k>=0 && k<=nk-1)
			{
				int index = i + ni*j + ni*nj*k;
				b[index] = -rhs(i,j,k);
			}
		});



		double tolerance;
		int iterations;
		bool success = AMGPCGSolve(fixed_matrix,b,x,A_L,R_L,P_L,S_L,total_level,1e-3,5,tolerance,iterations,_nx,_ny,_nz);
		printf("Solver took %d iterations and had residual %e\n", iterations, tolerance);
		if(!success) {
			printf("WARNING: Pressure solve failed!************************************************\n");
		}

		tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
		{
			int k = thread_idx/slice;
			int j = (thread_idx%slice)/ni;
			int i = thread_idx%ni;
			if(i>=0 && i<=ni-1 && j>=0 && j<=nj-1 && k>=0 && k<=nk-1)
			{
				int index = i + ni*j + ni*nj*k;
				psi(i,j,k) = x[index];
			}
		});

	}
	void vorticity_confinement(float coeff);
	void add_vort_confinement(buffer3Df & field, buffer3Df & vortf);
	void reaction(float dt);
	void decay_vortices(float dt, buffer3Df & wxnp1, buffer3Df &wynp1, buffer3Df &wznp1);
	void FLIP_advect(float dt);
	void FLIP_IVOCK(float dt);
	void clampSmoke();
	void clampExtrema(float dt, buffer3Df & f_n, buffer3Df & f_np1);
	void clampExtrema_order1(float dt, buffer3Df & f_n, buffer3Df & f_np1);
	void setSmoke(double temp_decay, double alpha,double beta,double smoke_heat,double smoke_dens){
		_temp_decay = temp_decay;
		_alpha = alpha;
		_beta = beta;
		_smoke_heat = smoke_heat;
		_smoke_dens = smoke_dens;
	}
	void BFECC(float dt);
	void BFECC_IVOCK(float dt);
	void Best_Combine(float dt);
	void BFECC_field(float dt, buffer3Df &field_n,buffer3Df &field_np1,buffer3Df &field_AUX);
	void BFECC_field_order1(float dt, buffer3Df &field_n,buffer3Df &field_np1,buffer3Df &field_AUX);

	void MacCormack(float dt);
	void MacCormack_IVOCK(float dt);
	void MacCormack_field(float dt, buffer3Df &field_n,buffer3Df &field_np1,buffer3Df &field_AUX);
	void MacCormack_field_order1(float dt, buffer3Df &field_n,buffer3Df &field_np1,buffer3Df &field_AUX);
	void getDivergence();

	void advect(float dt);
	void advect_boundary(float dt, buffer3Df & field, buffer3Df &field_new);
	void simple_advect_field(float dt, buffer3Df & field, buffer3Df &field_new);
	void simple_advect_field_order1(float dt, buffer3Df & field, buffer3Df &field_new);
	void advect_field(float dt, buffer3Df & field, buffer3Df &field_new);
	void advect_field_cubic(float dt, buffer3Df & field, buffer3Df &field_new);
	void advect_field_cubic_clamp(float dt, buffer3Df & field, buffer3Df &field_new);

	void iVICK(float dt);
	
	void addCurlPsi();
	void formRHS(float scale,float dt);
	void set_boundary(buffer3Dc & b_desc) {_b_desc.copy(b_desc);}
	void set_heat(buffer3Dc & h_desc)
	{
		_h_desc.copy(h_desc);
		
	}
	void compute_curl();
	void stretch(float dt);
	
	void time_step(float dt,int adv_type);
	void heat_decay(float dt);
	void gen_heat(float dt);
	void diffuse_heat(float dt);
	void diffuse_buffer(float c, buffer3Df & diffuse_field);
	void add_force(float dt);
	void projection();
	void pcg_projection();
	void set_bc();
	void compute_rhs(float scale);
	void apply_grad();
	void output(uint nx, uint ny, uint nz, int frame, char* file_path);
};





#endif