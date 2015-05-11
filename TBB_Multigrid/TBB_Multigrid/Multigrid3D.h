
#ifndef _MULTI_GRID_H_
#define _MULTI_GRID_H_


#include <vector>
#include "array.h"
#include "tbb/tbb.h"
#include "fluid_buffer3D.h"
using namespace std;

namespace gf{
struct uniform_grid_descriptor_3D
{
	int gridx, gridy, gridz, system_size;
	uniform_grid_descriptor_3D()
		:gridx(0), gridy(0), gridz(0), system_size(0)
	{}
	uniform_grid_descriptor_3D(int x, int y, int z)
		:gridx(x), gridy(y), gridz(z), system_size(x*y*z)
	{}
	~uniform_grid_descriptor_3D()
	{}
	void set_grid_information(int x, int y, int z)
	{
		gridx = x; gridy = y;  gridz = z;
		system_size = x*y*z;
	}
};
template<class T>
class MultiGrid
{
public:
	MultiGrid(void){ m_bInitialized = false;}

	~MultiGrid(void){m_FinalMemoryEachLevel();}

	void m_InitialSystem(int gridx, int gridy, int gridz);
	void m_SetGridInformation(int x_, int y_, int z_);
	void m_AssignMemoryEachLevel();
	void m_FinalMemoryEachLevel();
	void m_ComputeLevels();
	//void m_CreateIndexBuffers();
	bool m_bInitialized;
	uniform_grid_descriptor_3D m_system_descriptor;
	int m_max_level;
	std::vector<Buffer3D<T>*> xk;
	std::vector<Buffer3D<T>*> xk_new;
	std::vector<Buffer3D<T>*> bk;
	std::vector<Buffer3D<T>*> rk;
	std::vector<uniform_grid_descriptor_3D> systemk;
};

template<class T>
void MultiGrid<T>::m_FinalMemoryEachLevel()
{
	if(m_bInitialized)
	{
		for(int i=1; i<m_max_level; ++i)
		{
			xk[i]->free();
			xk_new[i]->free();
			bk[i]->free();
			rk[i]->free();
		}
		rk[0]->free();
		xk_new[0]->free();
		m_bInitialized = false;
	}
}
template<class T>
void MultiGrid<T>::m_AssignMemoryEachLevel()
{
	/*
		if grid has been changed
		as soon as we know what the finest level grid should be
		we can allocate memory for each level's grids
		first compute how many steps it is needed to get to the
		coarest level, a grid whose smallest dimension is 2; 
		suppose its level is 0;
		then, for level=0 to finest;
		memory[level]=alloc corresponding memory;
	*/
	
	xk.resize(m_max_level);
	xk_new.resize(m_max_level);
	bk.resize(m_max_level);
	rk.resize(m_max_level);

	for (int i=0; i<m_max_level; i++)
	{
		xk[i] = new Buffer3D<T>();
		xk_new[i] = new Buffer3D<T>();
		bk[i] = new Buffer3D<T>();
		rk[i] = new Buffer3D<T>();
	}
	systemk.resize(m_max_level);
	systemk[0].set_grid_information(m_system_descriptor.gridx,m_system_descriptor.gridy, m_system_descriptor.gridz);
	for(int i=1; i<m_max_level; ++i)
	{
		systemk[i].set_grid_information((systemk[i-1].gridx/2)==0?1:systemk[i-1].gridx/2,
										(systemk[i-1].gridy/2)==0?1:systemk[i-1].gridy/2, 
										(systemk[i-1].gridz/2)==0?1:systemk[i-1].gridz/2);
	}
	
	for(int i=1; i<m_max_level; ++i)
	{
		double h = 1.0/((float)(systemk[i].gridx));
		xk[i]->init(systemk[i].gridx, systemk[i].gridy, systemk[i].gridz,h,-0.5,-0.5,-0.5);
		xk_new[i]->init(systemk[i].gridx, systemk[i].gridy, systemk[i].gridz,h,-0.5,-0.5,-0.5);
		bk[i]->init(systemk[i].gridx, systemk[i].gridy, systemk[i].gridz,h,-0.5,-0.5,-0.5);
		rk[i]->init(systemk[i].gridx, systemk[i].gridy, systemk[i].gridz,h,-0.5,-0.5,-0.5);
	}
	xk_new[0]->init(systemk[0].gridx, systemk[0].gridy, systemk[0].gridz,1.0/((float)(systemk[0].gridx)),-0.5,-0.5,-0.5);
	rk[0]->init(systemk[0].gridx, systemk[0].gridy, systemk[0].gridz,1.0/((float)(systemk[0].gridx)),-0.5,-0.5,-0.5);
	m_bInitialized = true;
	
}
template<class T>
void MultiGrid<T>::m_SetGridInformation(int x_, int y_, int z_)
{
	m_system_descriptor.set_grid_information(x_, y_, z_);
}
template<class T>
void MultiGrid<T>::m_ComputeLevels()
{
	int level = 0;
	int x = m_system_descriptor.gridx, y = m_system_descriptor.gridy, z = m_system_descriptor.gridz;
	while (x*y*z>512)
	{
		level = level + 1;
		x = x/2; if(x==0) x=1;
		y = y/2; if(y==0) y=1;
		z = z/2; if(z==0) z=1;
	}
	m_max_level = level+1;
}
template<class T>
void MultiGrid<T>::m_InitialSystem(int gridx, int gridy, int gridz)
{
	m_FinalMemoryEachLevel();
	m_SetGridInformation(gridx, gridy, gridz);
	m_ComputeLevels();
	m_AssignMemoryEachLevel();
	//mg_InitIdxBuffer(systemk[m_max_level-1].gridx,
	//	systemk[m_max_level-1].gridy,
	//	systemk[m_max_level-1].gridz);
	//m_CreateIndexBuffers();
}


template<class T> 
class MultiGridSolver3D: public MultiGrid<T>
{
public:
	MultiGridSolver3D(){_b_Dirichlet = false; }
	~MultiGridSolver3D(){}
	bool _b_Dirichlet;
	void m_Vcycle(Buffer3D<T>* x, Buffer3D<T>* b, T tol, T &residual, int level);
	void m_FullMultiGrid(Buffer3D<T>* x, Buffer3D<T>* b, T tol, T &residual);
	

	
	void m_Red_Black_Gauss(uniform_grid_descriptor_3D &system, Buffer3D<T> & b, Buffer3D<T> & x, Buffer3D<T> & x_new, int iter_time, double omega);
	void m_Restrict(uniform_grid_descriptor_3D &next_level, uniform_grid_descriptor_3D &curr_level, Buffer3D<T> & src_level, Buffer3D<T> & dst_level);
	void m_MonopleToMonopole(uniform_grid_descriptor_3D &next_level, uniform_grid_descriptor_3D &curr_level, Buffer3D<T> & src_level, Buffer3D<T> & dst_level);
	void m_cumulate_poles(uniform_grid_descriptor_3D & curr_level,
		                  uniform_grid_descriptor_3D & upper_level,
						  Buffer3D<T> & curr_b,
						  Buffer3D<T> & upper_b,
						  Buffer3D<T> & curr_x,
						  T curr_h,
						  T upper_h);
	void m_Toplevel(Buffer3D<T> &x, Buffer3D<T> &b, T h);
	void m_bottomlevel(Buffer3D<T> &x, Buffer3D<T> &b, T h);
	void m_ExactSolve();
	void m_Prolongate(uniform_grid_descriptor_3D &next_level, uniform_grid_descriptor_3D &curr_level, Buffer3D<T> & src_level, Buffer3D<T> & dst_level);
	void m_ComputeResidual(uniform_grid_descriptor_3D &system, Buffer3D<T> & res, Buffer3D<T> & b, Buffer3D<T> & x);
	//void m_bottomup(Buffer3D<T>* x, Buffer3D<T>* b, T h, int level);
	void m_FastSummation(Buffer3D<T>* x, Buffer3D<T>* b, T h);
	T m_BarnesHutSummation(int level, int i, int j, int k, T curr_h, T h0, T wx, T wy, T wz);
	void m_applyOpenBoundaryCondition(Buffer3D<T> * mass, Buffer3D<T> * rhs, T h);
	
	
};
template<class T>
void MultiGridSolver3D<T>::m_applyOpenBoundaryCondition(Buffer3D<T> * mass, Buffer3D<T> * rhs, T h)
{
	xk[0] = rhs;
	bk[0] = mass;
	float *hk = new float[m_max_level];
	hk[0] = (float)h;
	for(int i=1; i<m_max_level; i++)
	{
		xk[i]->setZero();
		bk[i]->setZero();
		rk[i]->setZero();
	}
	rk[0]->setZero();

	//bottom-up pass
	for(int i=0; i<m_max_level-1; ++i)
	{
		(*(rk[i])).copy(*(bk[i]));

		m_MonopleToMonopole(systemk[i+1], 
			systemk[i],
			*(rk[i]),
			*(bk[i+1]));
		hk[i+1] = 2.0*hk[i];
	}
	int compute_element = systemk[0].gridz * systemk[0].gridy;
	
	//left face
	//for (int k=0;k<systemk[0].gridz;k++)for (int j=0; j<systemk[0].gridy;j++)
	tbb::parallel_for(0,compute_element,1,[&](int thread_idx)
	{
		int j = thread_idx%(systemk[0].gridy);
		int k = thread_idx/systemk[0].gridy;
		float wx = -0.5*h;
		float wy = ((float)j + 0.5)*h;
		float wz = ((float)k + 0.5)*h;
		T potential = 0.0;
		for (int kk=0;kk<systemk[m_max_level-1].gridz; kk++)
			for (int jj=0;jj<systemk[m_max_level-1].gridy; jj++)
				for (int ii=0;ii<systemk[m_max_level-1].gridx; ii++)
				{
					potential += m_BarnesHutSummation(m_max_level-1,ii,jj,kk,hk[m_max_level-1],h,wx,wy,wz);
				}
				(*(xk[0]))(0,j,k) -= potential;
	});


	//right face
	compute_element = systemk[0].gridz * systemk[0].gridy;
	//for (int k=0;k<systemk[0].gridz;k++)for (int j=0; j<systemk[0].gridy;j++)
	tbb::parallel_for(0,compute_element,1,[&](int thread_idx)
	{
		int j = thread_idx%(systemk[0].gridy);
		int k = thread_idx/systemk[0].gridy;
		float wx = ((float)(systemk[0].gridx)+0.5)*h;
		float wy = ((float)j + 0.5)*h;
		float wz = ((float)k + 0.5)*h;
		T potential = 0.0;
		for (int kk=0;kk<systemk[m_max_level-1].gridz; kk++)
			for (int jj=0;jj<systemk[m_max_level-1].gridy; jj++)
				for (int ii=0;ii<systemk[m_max_level-1].gridx; ii++)
				{
					potential += m_BarnesHutSummation(m_max_level-1,ii,jj,kk,hk[m_max_level-1],h,wx,wy,wz);
				}
				(*(xk[0]))(systemk[0].gridx-1,j,k) -= potential;
	});


	//top face
	compute_element = systemk[0].gridz * systemk[0].gridx;
	//for (int k=0;k<systemk[0].gridz;k++)for (int i=0; i<systemk[0].gridx;i++)
	tbb::parallel_for(0,compute_element,1,[&](int thread_idx)
	{
		int i = thread_idx%systemk[0].gridx;
		int k = thread_idx/systemk[0].gridx;
		float wx = ((float)i + 0.5)*h;
		float wy = ((float)(systemk[0].gridy)+0.5)*h;
		float wz = ((float)k + 0.5)*h;
		T potential = 0.0;
		for (int kk=0;kk<systemk[m_max_level-1].gridz; kk++)
			for (int jj=0;jj<systemk[m_max_level-1].gridy; jj++)
				for (int ii=0;ii<systemk[m_max_level-1].gridx; ii++)
				{
					potential += m_BarnesHutSummation(m_max_level-1,ii,jj,kk,hk[m_max_level-1],h,wx,wy,wz);
				}
				(*(xk[0]))(i,systemk[0].gridy-1,k) -= potential;
	});
	//bottom face
	compute_element = systemk[0].gridz * systemk[0].gridx;
	//for (int k=0;k<systemk[0].gridz;k++)for (int i=0; i<systemk[0].gridx;i++)
	tbb::parallel_for(0,compute_element,1,[&](int thread_idx)
	{
		int i = thread_idx%systemk[0].gridx;
		int k = thread_idx/systemk[0].gridx;
		float wx = ((float)i + 0.5)*h;
		float wy = (-0.5)*h;
		float wz = ((float)k + 0.5)*h;
		T potential = 0.0;
		for (int kk=0;kk<systemk[m_max_level-1].gridz; kk++)
			for (int jj=0;jj<systemk[m_max_level-1].gridy; jj++)
				for (int ii=0;ii<systemk[m_max_level-1].gridx; ii++)
				{
					potential += m_BarnesHutSummation(m_max_level-1,ii,jj,kk,hk[m_max_level-1],h,wx,wy,wz);
				}
				(*(xk[0]))(i,0,k) -= potential;
	});
	//front face
	compute_element = systemk[0].gridy * systemk[0].gridx;
	//for (int j=0;j<systemk[0].gridy;j++)for (int i=0; i<systemk[0].gridx;i++)
	tbb::parallel_for(0,compute_element,1,[&](int thread_idx)
	{
		int i = thread_idx%systemk[0].gridx;
		int j = thread_idx/systemk[0].gridx;
		float wx = ((float)i + 0.5)*h;
		float wy = ((float)j + 0.5)*h;
		float wz = (-0.5)*h;
		T potential = 0.0;
		for (int kk=0;kk<systemk[m_max_level-1].gridz; kk++)
			for (int jj=0;jj<systemk[m_max_level-1].gridy; jj++)
				for (int ii=0;ii<systemk[m_max_level-1].gridx; ii++)
				{
					potential += m_BarnesHutSummation(m_max_level-1,ii,jj,kk,hk[m_max_level-1],h,wx,wy,wz);
				}
				(*(xk[0]))(i,j,0) -= potential;
	});
	//back face
	compute_element = systemk[0].gridy * systemk[0].gridx;
	//for (int j=0;j<systemk[0].gridy;j++)for (int i=0; i<systemk[0].gridx;i++)
	tbb::parallel_for(0,compute_element,1,[&](int thread_idx)
	{
		int i = thread_idx%systemk[0].gridx;
		int j = thread_idx/systemk[0].gridx;
		float wx = ((float)i + 0.5)*h;
		float wy = ((float)j + 0.5)*h;
		float wz = ((float)(systemk[0].gridz) + 0.5)*h;
		T potential = 0.0;
		for (int kk=0;kk<systemk[m_max_level-1].gridz; kk++)
			for (int jj=0;jj<systemk[m_max_level-1].gridy; jj++)
				for (int ii=0;ii<systemk[m_max_level-1].gridx; ii++)
				{
					potential += m_BarnesHutSummation(m_max_level-1,ii,jj,kk,hk[m_max_level-1],h,wx,wy,wz);
				}
				(*(xk[0]))(i,j,systemk[0].gridz-1) -= potential;
	});
	

	for(int i=0+1; i<m_max_level; i++)
	{
		xk[i]->setZero();
		bk[i]->setZero();
		rk[i]->setZero();

	}
	delete [] hk;
}
template<class T>
T MultiGridSolver3D<T>::m_BarnesHutSummation(int level, int i, int j, int k, T curr_h, T h0, T wx, T wy, T wz)
{
	T px = ((float)i + 0.5)*curr_h;
	T py = ((float)j + 0.5)*curr_h;
	T pz = ((float)k + 0.5)*curr_h;

	T dx = fabs(px - wx);
	T dy = fabs(py - wy);
	T dz = fabs(pz - wz);
	T test_d = max(max(dx,dy),dz);
	//if level is 0 or far enough, return potential value at this point.
	if(level==0 || test_d > 1.001*(0.5*curr_h + 0.5*h0))
	{
		T r = sqrt(dx*dx + dy*dy + dz*dz);
		if(r>1e-7)
			return (*(bk[level]))(i,j,k)*0.07957747154/r;
		else
			return 0;
	}
	//else loop over 8 children, do BarnesHutSummation
	else
	{
		T sum = 0;
		for (int kk=0;kk<=1;kk++)
		for (int jj=0;jj<=1;jj++)
		for (int ii=0;ii<=1;ii++)
		{
			sum += m_BarnesHutSummation(level-1, i*2+ii, j*2+jj, k*2+kk, curr_h*0.5, h0,wx,wy,wz);
		}
		return sum;
	}

}



template<class T>
void MultiGridSolver3D<T>::m_Toplevel(Buffer3D<T> &x, Buffer3D<T> &b, T h)
{
	int computational_elements = x._blockx*x._blocky*x._blockz;
	uint slice = x._blockx*x._blocky;

	tbb::parallel_for(0, computational_elements, 1, [&](int thread_idx) 
	{

		int bbk = thread_idx/slice;
		int bbj = (thread_idx%slice)/x._blockx;
		int bbi = thread_idx%(x._blockx);

		for (int kk=0;kk<8;kk++)for(int jj=0;jj<8;jj++)for(int ii=0;ii<8;ii++)
		{
			int i=bbi*8+ii, j=bbj*8+jj, k=bbk*8+kk;
			if( i<x._nx
				&& j<x._ny 
				&& k<x._nz)
			{
				float px = ((float)i + 0.5)*h, py = ((float)j + 0.5)*h, pz = ((float)k + 0.5)*h;
				float sum = 0.0;
				for (int kkk=0;kkk<b._nz;kkk++)
				for (int jjj=0;jjj<b._ny;jjj++)
				for (int iii=0;iii<b._nx;iii++)
				{
					if(!(abs(iii-i)<=1&&abs(jjj-j)<=1&&abs(kkk-k)<=1))
					{
						float px2 = ((float)iii + 0.5)*h;
						float py2 = ((float)jjj + 0.5)*h;
						float pz2 = ((float)kkk + 0.5)*h;
						float dx = px2-px;
						float dy = py2-py;
						float dz = pz2-pz;
						float len = sqrt(dx*dx+dy*dy+dz*dz);

						sum += 0.07957747154/len*b(iii,jjj,kkk);
					}
				}
				x(i,j,k) = sum;
			}
		}
	});
}
template<class T>
void MultiGridSolver3D<T>::m_bottomlevel(Buffer3D<T> &x, Buffer3D<T> &b, T h)
{
	int computational_elements = x._blockx*x._blocky*x._blockz;
	uint slice = x._blockx*x._blocky;

	tbb::parallel_for(0, computational_elements, 1, [&](int thread_idx) 
	{

		int bbk = thread_idx/slice;
		int bbj = (thread_idx%slice)/x._blockx;
		int bbi = thread_idx%(x._blockx);

		for (int kk=0;kk<8;kk++)for(int jj=0;jj<8;jj++)for(int ii=0;ii<8;ii++)
		{
			int i=bbi*8+ii, j=bbj*8+jj, k=bbk*8+kk;
			if( i<x._nx
				&& j<x._ny 
				&& k<x._nz)
			{
				float px = ((float)i + 0.5)*h, py = ((float)j + 0.5)*h, pz = ((float)k + 0.5)*h;
				float sum = 0.0;
				for (int kkk=max(0,k-1);kkk<=min(b._nz-1,k+1);kkk++)
					for (int jjj=max(0,j-1);jjj<=min(b._ny-1,j+1);jjj++)
						for (int iii=max(0,i-1);iii<=min(b._nx-1,i+1);iii++)
						{
							if(!(iii==i && jjj==j && kkk==k))
							{

								float px2 = ((float)iii + 0.5)*h;
								float py2 = ((float)jjj + 0.5)*h;
								float pz2 = ((float)kkk + 0.5)*h;
								float dx = px2-px;
								float dy = py2-py;
								float dz = pz2-pz;
								float len = sqrt(dx*dx+dy*dy+dz*dz);

								sum += 0.07957747154/len*b(iii,jjj,kkk);
							}
							
						}
						x(i,j,k) += sum;
			}
		}
	});
}
template<class T>
void MultiGridSolver3D<T>::m_cumulate_poles(uniform_grid_descriptor_3D & curr_level, 
	uniform_grid_descriptor_3D & upper_level, 
	Buffer3D<T> & curr_b, 
	Buffer3D<T> & upper_b, 
	Buffer3D<T> & curr_x, 
	T curr_h, T upper_h)
{
	//for all x in current level
	//loop over all neighbor 1 cells of the upper level b
	//add their contribution,

	int computational_elements = curr_x._blockx*curr_x._blocky*curr_x._blockz;
	uint slice = curr_x._blockx*curr_x._blocky;

	tbb::parallel_for(0, computational_elements, 1, [&](int thread_idx) 
	{

		int bbk = thread_idx/slice;
		int bbj = (thread_idx%slice)/curr_x._blockx;
		int bbi = thread_idx%(curr_x._blockx);

		for (int kk=0;kk<8;kk++)for(int jj=0;jj<8;jj++)for(int ii=0;ii<8;ii++)
		{
			int i=bbi*8+ii, j=bbj*8+jj, k=bbk*8+kk;
			if( i<curr_x._nx
				&& j<curr_x._ny 
				&& k<curr_x._nz)
			{
				int upper_i = i/2;
				int upper_j = j/2;
				int upper_k = k/2;
				float px = ((float)i + 0.5)*curr_h;
				float py = ((float)j + 0.5)*curr_h;
				float pz = ((float)k + 0.5)*curr_h;
				float sum = 0.0;
				for (int kkk=max(0,upper_k-1);kkk<=min(upper_b._nz-1,upper_k+1);kkk++){

				for (int jjj=max(0,upper_j-1);jjj<=min(upper_b._ny-1,upper_j+1);jjj++){

				for (int iii=max(0,upper_i-1);iii<=min(upper_b._nx-1,upper_i+1);iii++)
				{
					if(!(iii==upper_i && jjj==upper_j && kkk==upper_k))
					{


						float px2 = ((float)iii + 0.5)*upper_h;
						float py2 = ((float)jjj + 0.5)*upper_h;
						float pz2 = ((float)kkk + 0.5)*upper_h;
						float dx = fabs(px2-px);
						float dy = fabs(py2-py);
						float dz = fabs(pz2-pz);
						float test_d = max(dx,max(dy,dz));
						//here we have two cases,
						//a) the upper cell is not connected to current small cell
						if(test_d > 1.505 * curr_h){

							float len = sqrt(dx*dx+dy*dy+dz*dz);

							sum += 0.07957747154/len*upper_b(iii,jjj,kkk);
						}
						else //b) the upper cell is connected to the small cell
						{
							//consider the eight sub cell,
							for (int kkkk=0;kkkk<=1;kkkk++)
							for (int jjjj=0;jjjj<=1;jjjj++)
							for (int iiii=0;iiii<=1;iiii++)
							{
								int sub_i = iii*2 + iiii;
								int sub_j = jjj*2 + jjjj;
								int sub_k = kkk*2 + kkkk;
								int di = abs(sub_i-i);
								int dj = abs(sub_j-j);
								int dk = abs(sub_k-k);
								if(!(di<=1&&dj<=1&&dk<=1))
								{
									float px3 = ((float)sub_i + 0.5)*curr_h;
									float py3 = ((float)sub_j + 0.5)*curr_h;
									float pz3 = ((float)sub_k + 0.5)*curr_h;
									float dx2 = fabs(px3-px);
									float dy2 = fabs(py3-py);
									float dz2 = fabs(pz3-pz);
									float len = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);

									sum += 0.07957747154/len*curr_b(sub_i,sub_j,sub_k);
								}
							}
						}
					}

				}}}//end for kkk,jjj,iii
				curr_x(i,j,k) += sum;
			}
		}
	});



}
template<class T>
void MultiGridSolver3D<T>::m_FastSummation(Buffer3D<T>* x, Buffer3D<T>* b, T h)
{
	xk[0] = x;
	bk[0] = b;
	float *hk = new float[m_max_level];
	hk[0] = (float)h;
	for(int i=1; i<m_max_level; i++)
	{
		xk[i]->setZero();
		bk[i]->setZero();
		rk[i]->setZero();
	}
	rk[0]->setZero();

	//bottom-up pass
	for(int i=0; i<m_max_level-1; ++i)
	{
		(*(rk[i])).copy(*(bk[i]));

		m_MonopleToMonopole(systemk[i+1], 
			systemk[i],
			*(rk[i]),
			*(bk[i+1]));
		hk[i+1] = 2.0*hk[i];
	}
	m_Toplevel(*(xk[m_max_level-1]), *(bk[m_max_level-1]),hk[m_max_level-1]);
	//top-down pass
	for (int i=m_max_level-2; i>=0; i--)
	{
		m_Prolongate(systemk[i], 
			systemk[i+1],
			*(xk[i+1]),
			*(xk[i]));
		m_cumulate_poles(systemk[i],
			systemk[i+1],
			*(bk[i]),
			*(bk[i+1]),
			*(xk[i]),hk[i],hk[i+1]);
	}
	m_bottomlevel(*(xk[0]), *(bk[0]),hk[0]);

	for(int i=0+1; i<m_max_level; i++)
	{
		xk[i]->setZero();
		bk[i]->setZero();
		rk[i]->setZero();

	}
	delete [] hk;
}

//
//template<class T>
//void MultiGridSolver3D<T>::form_openBC(Buffer3D<T> & b0, Buffer3D<T> & bTop, T h0, T hTop)
//{
//	int computational_elements = bTop._blockx*bTop._blocky*bTop._blockz;
//	uint slice = bTop._blockx*bTop._blocky;
//
//	tbb::parallel_for(0, computational_elements, 1, [&](int thread_idx) 
//	{
//
//		int bbk = thread_idx/slice;
//		uint bbj = (thread_idx%slice)/bTop._blockx;
//		uint bbi = thread_idx%(bTop._blockx);
//
//		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
//		{
//			uint i=bbi*8+ii, j=bbj*8+jj, k=bbk*8+kk;
//			if( i<systemk[m_max_level-1].gridx 
//				&& j<systemk[m_max_level-1].gridy 
//				&& k<systemk[m_max_level-1].gridz)
//			{
//				bTop(i,j,k) = -bTop(i,j,k)*hTop;//now this is the mass
//			}
//		}
//	});
//
//	computational_elements = b0._blockx*b0._blocky*b0._blockz;
//	slice = b0._blockx*b0._blocky;
//
//	tbb::parallel_for(0, computational_elements, 1, [&](int thread_idx) 
//	{
//
//		int bbk = thread_idx/slice;
//		uint bbj = (thread_idx%slice)/b0._blockx;
//		uint bbi = thread_idx%(b0._blockx);
//
//		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
//		{
//			uint i=bbi*8+ii, j=bbj*8+jj, k=bbk*8+kk;
//			if( i<b0._nx 
//				&& j<b0._ny 
//				&& k<b0._nz)
//			{
//				if (i==0 )
//				{
//					float px = ((float)i - 0.5)*h0, py = ((float)j + 0.5)*h0, pz = ((float)k + 0.5)*h0;
//					float sum = 0;
//					for(int kkk=0;kkk<bTop._nz;kkk++)
//						for(int jjj = 0; jjj<bTop._ny;jjj++)
//							for (int iii=0;iii<bTop._nz;iii++)
//							{
//								float px2 = ((float)iii + 0.5)*hTop;
//								float py2 = ((float)jjj + 0.5)*hTop;
//								float pz2 = ((float)kkk + 0.5)*hTop;
//
//								float dx = px2-px;
//								float dy = py2-py;
//								float dz = pz2-pz;
//								float len = sqrt(dx*dx+dy*dy+dz*dz);
//
//								sum += 0.07957747154/len*bTop(iii,jjj,kkk);
//							}
//							b0(i,j,k) -= sum;
//				}
//				if (i==b0._nx-1)
//				{
//					float px = ((float)i + 1.5)*h0, py = ((float)j + 0.5)*h0, pz = ((float)k + 0.5)*h0;
//					float sum = 0;
//					for(int kkk=0;kkk<bTop._nz;kkk++)
//						for(int jjj = 0; jjj<bTop._ny;jjj++)
//							for (int iii=0;iii<bTop._nz;iii++)
//							{
//								float px2 = ((float)iii + 0.5)*hTop;
//								float py2 = ((float)jjj + 0.5)*hTop;
//								float pz2 = ((float)kkk + 0.5)*hTop;
//
//								float dx = px2-px;
//								float dy = py2-py;
//								float dz = pz2-pz;
//								float len = sqrt(dx*dx+dy*dy+dz*dz);
//
//								sum += 0.07957747154/len*bTop(iii,jjj,kkk);
//							}
//							b0(i,j,k) -= sum;
//				}
//				if (j==0)
//				{
//					float px = ((float)i + 0.5)*h0, py = ((float)j - 0.5)*h0, pz = ((float)k + 0.5)*h0;
//					float sum = 0;
//					for(int kkk=0;kkk<bTop._nz;kkk++)
//						for(int jjj = 0; jjj<bTop._ny;jjj++)
//							for (int iii=0;iii<bTop._nz;iii++)
//							{
//								float px2 = ((float)iii + 0.5)*hTop;
//								float py2 = ((float)jjj + 0.5)*hTop;
//								float pz2 = ((float)kkk + 0.5)*hTop;
//
//								float dx = px2-px;
//								float dy = py2-py;
//								float dz = pz2-pz;
//								float len = sqrt(dx*dx+dy*dy+dz*dz);
//
//								sum += 0.07957747154/len*bTop(iii,jjj,kkk);
//							}
//							b0(i,j,k) -= sum;
//				}
//				if (j==b0._ny-1)
//				{
//					float px = ((float)i + 0.5)*h0, py = ((float)j + 1.5)*h0, pz = ((float)k + 0.5)*h0;
//					float sum = 0;
//					for(int kkk=0;kkk<bTop._nz;kkk++)
//						for(int jjj = 0; jjj<bTop._ny;jjj++)
//							for (int iii=0;iii<bTop._nz;iii++)
//							{
//								float px2 = ((float)iii + 0.5)*hTop;
//								float py2 = ((float)jjj + 0.5)*hTop;
//								float pz2 = ((float)kkk + 0.5)*hTop;
//
//								float dx = px2-px;
//								float dy = py2-py;
//								float dz = pz2-pz;
//								float len = sqrt(dx*dx+dy*dy+dz*dz);
//
//								sum += 0.07957747154/len*bTop(iii,jjj,kkk);
//							}
//							b0(i,j,k) -= sum;
//				}
//
//				if (k==0)
//				{
//					float px = ((float)i + 0.5)*h0, py = ((float)j - 0.5)*h0, pz = ((float)k - 0.5)*h0;
//					float sum = 0;
//					for(int kkk=0;kkk<bTop._nz;kkk++)
//						for(int jjj = 0; jjj<bTop._ny;jjj++)
//							for (int iii=0;iii<bTop._nz;iii++)
//							{
//								float px2 = ((float)iii + 0.5)*hTop;
//								float py2 = ((float)jjj + 0.5)*hTop;
//								float pz2 = ((float)kkk + 0.5)*hTop;
//
//								float dx = px2-px;
//								float dy = py2-py;
//								float dz = pz2-pz;
//								float len = sqrt(dx*dx+dy*dy+dz*dz);
//
//								sum += 0.07957747154/len*bTop(iii,jjj,kkk);
//							}
//							b0(i,j,k) -= sum;
//				}
//				if (k==b0._nz-1)
//				{
//					float px = ((float)i + 0.5)*h0, py = ((float)j + 1.5)*h0, pz = ((float)k + 1.5)*h0;
//					float sum = 0;
//					for(int kkk=0;kkk<bTop._nz;kkk++)
//						for(int jjj = 0; jjj<bTop._ny;jjj++)
//							for (int iii=0;iii<bTop._nz;iii++)
//							{
//								float px2 = ((float)iii + 0.5)*hTop;
//								float py2 = ((float)jjj + 0.5)*hTop;
//								float pz2 = ((float)kkk + 0.5)*hTop;
//
//								float dx = px2-px;
//								float dy = py2-py;
//								float dz = pz2-pz;
//								float len = sqrt(dx*dx+dy*dy+dz*dz);
//
//								sum += 0.07957747154/len*bTop(iii,jjj,kkk);
//							}
//							b0(i,j,k) -= sum;
//				}
//
//
//
//			}
//		}
//	});
//}

//template<class T>
//void MultiGridSolver3D<T>::m_bottomup(Buffer3D<T>* x, Buffer3D<T>* b, T h, int level)
//{
//	xk[level] = x;
//	bk[level] = b;
//
//	for(int i=level+1; i<m_max_level; i++)
//	{
//		xk[i]->setZero();
//		bk[i]->setZero();
//		rk[i]->setZero();
//	}
//	rk[level]->setZero();
//
//	float grid_h = (float)h;
//	for(int i=level; i<1; ++i)
//	{
//		m_Red_Black_Gauss(systemk[i],
//			*(bk[i]),
//			*(xk[i]),
//			*(xk_new[i]),
//			0,
//			1.0);
//
//		m_ComputeResidual(systemk[i],
//			*(rk[i]),
//			*(bk[i]),
//			*(xk[i]));
//
//		m_Restrict(systemk[i+1], 
//			systemk[i],
//			*(rk[i]),
//			*(bk[i+1]));
//		grid_h = grid_h*2.0;
//	}
//	
//	form_openBC(*(bk[0]),
//		        *(bk[1]),
//				h,
//				grid_h);
//	
//}
template<class T>
void MultiGridSolver3D<T>::m_Red_Black_Gauss(uniform_grid_descriptor_3D &system,
	Buffer3D<T> & b, Buffer3D<T> & x, Buffer3D<T> & x_new, int iter_time, double omega)
{
	int computational_elements = x._blockx*x._blocky*x._blockz;
	uint slice = x._blockx*x._blocky;
	//x_new.copy(x);
	for(int iteration=0; iteration<iter_time; iteration++)
	{

		//red sweep
		//tbb::parallel_for(0, computational_elements, 1, [&](int thread_idx) {

		//	uint bk = thread_idx/slice;
		//	uint bj = (thread_idx%slice)/x._blockx;
		//	uint bi = thread_idx%(x._blockx);

		//	for (uint kk=0;kk<x._blockN;kk++)for(uint jj=0;jj<x._blockN;jj++)for(uint ii=0;ii<x._blockN;ii++)
		//	{
		//		uint i=bi*x._blockN+ii, j=bj*x._blockN+jj, k=bk*x._blockN+kk;
		//		if((j+k+i)%2==0 && i<system.gridx && j<system.gridy && k<system.gridz)
		//		{

		//			//T rhs = b(i,j,k);
		//			uint li = (i==0)? i:i-1;
		//			uint ri = (i==system.gridx-1)? i : i+1;
		//			uint ti = (j==system.gridy-1)? j : j+1;
		//			uint di = (j==0)? j:j-1;
		//			uint fi = (k==0)? k:k-1;
		//			uint bi = (k==system.gridz-1)?k:k+1;


		//			T lv = 0;
		//			if(ii>0) lv = x(li,j,k);
		//			T rv = 0;
		//			if(ii<x._blockN-1) rv = x(ri,j,k);
		//			T tv = 0;
		//			if(jj<x._blockN-1) tv = x(i,ti,k);
		//			T dv = 0;
		//			if(jj>0) dv = x(i,di,k);
		//			T fv = 0;
		//			if(kk>0) fv = x(i,j,fi);
		//			T bv = 0;
		//			if(kk<x._blockN-1) bv = x(i,j,bi);
		//			if(_b_Dirichlet){

		//				if(i==0) lv = 0;
		//				if(i==system.gridx-1) rv = 0;
		//				if(j==system.gridy-1) tv = 0;
		//				if(j==0) dv = 0;
		//				if(k==0) fv = 0;
		//				if(k==system.gridz-1) bv = 0;
		//			}

		//			T	rr =  lv + rv;
		//			rr += tv + dv;
		//			rr += fv + bv;
		//			//rr-=rhs;

		//			x_new(i,j,k) = rr;

		//		}
		//	}
		//	for (uint kk=0;kk<x._blockN;kk++)for(uint jj=0;jj<x._blockN;jj++) //ii=0
		//	{
		//		uint i=bi*x._blockN, j=bj*x._blockN+jj, k=bk*x._blockN+kk;
		//		if((j+k+i)%2==0 && i<system.gridx && j<system.gridy && k<system.gridz)
		//		{

		//			//T rhs = b(i,j,k);
		//			uint li = (i==0)? i:i-1;
		//			T lv = x(li,j,k);
		//	
		//			if(_b_Dirichlet){
		//				if(i==0) lv = 0;
		//			}

		//			x_new(i,j,k) += lv;

		//		}
		//	}
		//	for (uint kk=0;kk<x._blockN;kk++)for(uint jj=0;jj<x._blockN;jj++) //ii=x._blockN-1
		//	{
		//		uint i=bi*x._blockN+x._blockN-1, j=bj*x._blockN+jj, k=bk*x._blockN+kk;
		//		if((j+k+i)%2==0 && i<system.gridx && j<system.gridy && k<system.gridz)
		//		{

		//			//T rhs = b(i,j,k);
		//			uint ri = (i==system.gridx-1)? i : i+1;
		//			T rv = x(ri,j,k);

		//			if(_b_Dirichlet){
		//				if(i==system.gridx-1) rv = 0;
		//			}

		//			x_new(i,j,k) += rv;

		//		}
		//	}


		//	for (uint kk=0;kk<x._blockN;kk++)for(uint ii=0;ii<x._blockN;ii++) //jj=0
		//	{
		//		uint i=bi*x._blockN+ii, j=bj*x._blockN+0, k=bk*x._blockN+kk;
		//		if((j+k+i)%2==0 && i<system.gridx && j<system.gridy && k<system.gridz)
		//		{

		//			//T rhs = b(i,j,k);
		//			uint di = (j==0)? j:j-1;
		//			T dv = x(i,di,k);

		//			if(_b_Dirichlet){
		//				if(j==0) dv = 0;
		//			}

		//			x_new(i,j,k) += dv;

		//		}
		//	}
		//	for (uint kk=0;kk<x._blockN;kk++)for(uint ii=0;ii<x._blockN;ii++) //jj=x._blockN-1
		//	{
		//		uint i=bi*x._blockN+ii, j=bj*x._blockN+x._blockN-1, k=bk*x._blockN+kk;
		//		if((j+k+i)%2==0 && i<system.gridx && j<system.gridy && k<system.gridz)
		//		{

		//			//T rhs = b(i,j,k);
		//			uint ti = (j==system.gridy-1)? j : j+1;
		//			T tv = x(i,ti,k);

		//			if(_b_Dirichlet){
		//				if(j==system.gridy-1) tv = 0;
		//			}

		//			x_new(i,j,k) += tv;

		//		}
		//	}
		//	for (uint jj=0;jj<x._blockN;jj++)for(uint ii=0;ii<x._blockN;ii++) //kk=0
		//	{
		//		uint i=bi*x._blockN+ii, j=bj*x._blockN+jj, k=bk*x._blockN+0;
		//		if((j+k+i)%2==0 && i<system.gridx && j<system.gridy && k<system.gridz)
		//		{

		//			//T rhs = b(i,j,k);
		//			uint fi = (k==0)? k:k-1;
		//			
		//			T fv = x(i,j,fi);

		//			if(_b_Dirichlet){
		//				if(k==0) fv = 0;
		//			}

		//			x_new(i,j,k) += fv;

		//		}
		//	}
		//	for (uint jj=0;jj<x._blockN;jj++)for(uint ii=0;ii<x._blockN;ii++) //kk=x._blockN-1
		//	{
		//		uint i=bi*x._blockN+ii, j=bj*x._blockN+jj, k=bk*x._blockN+x._blockN-1;
		//		if((j+k+i)%2==0 && i<system.gridx && j<system.gridy && k<system.gridz)
		//		{


		//			uint bi = (k==system.gridz-1)?k:k+1;
		//			T bv = x(i,j,bi);

		//			if(_b_Dirichlet){
		//				if(k==system.gridz-1) bv = 0;
		//			}

		//			x_new(i,j,k) += bv;

		//		}
		//	}
		//	for (uint kk=0;kk<x._blockN;kk++)for(uint jj=0;jj<x._blockN;jj++)for(uint ii=0;ii<x._blockN;ii++)
		//	{
		//		uint i=bi*x._blockN+ii, j=bj*x._blockN+jj, k=bk*x._blockN+kk;
		//		if((j+k+i)%2==0 && i<system.gridx && j<system.gridy && k<system.gridz)
		//		{

		//			T rhs = b(i,j,k);
		//			x(i,j,k) = (1-omega)*x(i,j,k) + omega * (x_new(i,j,k)-rhs)/6.0;
		//		}
		//	}
		//	
		//});
		//black sweep
		tbb::parallel_for(0, computational_elements, 1, [&](int thread_idx) 
		{

			int bk = thread_idx/slice;
			uint bj = (thread_idx%slice)/x._blockx;
			uint bi = thread_idx%(x._blockx);

			for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
			{
				uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
				if((j+k+i)%2==1 && i<system.gridx && j<system.gridy && k<system.gridz)
				{

					T rhs = b(i,j,k);
					uint li = (i==0)? i:i-1;
					uint ri = (i==system.gridx-1)? i : i+1;
					uint ti = (j==system.gridy-1)? j : j+1;
					uint di = (j==0)? j:j-1;
					uint fi = (k==0)? k:k-1;
					uint bi = (k==system.gridz-1)?k:k+1;


					T lv = x(li,j,k);
					T rv = x(ri,j,k);
					T tv = x(i,ti,k);
					T dv = x(i,di,k);
					T fv = x(i,j,fi);
					T bv = x(i,j,bi);
					if(_b_Dirichlet){

						if(i==0) lv = 0;
						if(i==system.gridx-1) rv = 0;
						if(j==system.gridy-1) tv = 0;
						if(j==0) dv = 0;
						if(k==0) fv = 0;
						if(k==system.gridz-1) bv = 0;
					}

					T	rr =  lv + rv;
					rr += tv + dv;
					rr += fv + bv;
					rr-=rhs;

					x(i,j,k) = (1-omega)*x(i,j,k) + omega * rr/6.0;

				}
			}

		});
		tbb::parallel_for(0, computational_elements, 1, [&](int thread_idx) 
		{

			int bk = thread_idx/slice;
			uint bj = (thread_idx%slice)/x._blockx;
			uint bi = thread_idx%(x._blockx);

			for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
			{
				uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
				if((j+k+i)%2==0 && i<system.gridx && j<system.gridy && k<system.gridz)
				{

					T rhs = b(i,j,k);
					uint li = (i==0)? i:i-1;
					uint ri = (i==system.gridx-1)? i : i+1;
					uint ti = (j==system.gridy-1)? j : j+1;
					uint di = (j==0)? j:j-1;
					uint fi = (k==0)? k:k-1;
					uint bi = (k==system.gridz-1)?k:k+1;


					T lv = x(li,j,k);
					T rv = x(ri,j,k);
					T tv = x(i,ti,k);
					T dv = x(i,di,k);
					T fv = x(i,j,fi);
					T bv = x(i,j,bi);
					if(_b_Dirichlet){

						if(i==0) lv = 0;
						if(i==system.gridx-1) rv = 0;
						if(j==system.gridy-1) tv = 0;
						if(j==0) dv = 0;
						if(k==0) fv = 0;
						if(k==system.gridz-1) bv = 0;
					}

					T	rr =  lv + rv;
					rr += tv + dv;
					rr += fv + bv;
					rr-=rhs;

					x(i,j,k) = (1-omega)*x(i,j,k) + omega * rr/6.0;

				}
			}

		});
		
	
	}

}
template<class T>
void MultiGridSolver3D<T>::m_ComputeResidual(uniform_grid_descriptor_3D &system, Buffer3D<T> & res, Buffer3D<T> & b, Buffer3D<T> & x)
{
	int computational_elements = x._blockx*x._blocky*x._blockz;
	uint slice = x._blockx*x._blocky;

	
		//red sweep
	tbb::parallel_for(0, computational_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/x._blockx;
		uint bi = thread_idx%(x._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<system.gridx && j<system.gridy && k<system.gridz)
			{
			
				T rhs = b(i,j,k);
				T center = x(i,j,k);
				uint li = (i==0)? i:i-1;
				uint ri = (i==system.gridx-1)? i : i+1;
				uint ti = (j==system.gridy-1)? j : j+1;
				uint di = (j==0)? j:j-1;
				uint fi = (k==0)? k:k-1;
				uint bi = (k==system.gridz-1)?k:k+1;


				T lv = x(li,j,k);
				T rv = x(ri,j,k);
				T tv = x(i,ti,k);
				T dv = x(i,di,k);
				T fv = x(i,j,fi);
				T bv = x(i,j,bi);
				if(_b_Dirichlet){

					if(i==0) lv = 0;
					if(i==system.gridx-1) rv = 0;
					if(j==system.gridy-1) tv = 0;
					if(j==0) dv = 0;
					if(k==0) fv = 0;
					if(k==system.gridz-1) bv = 0;
				}

				T	rr =  lv + rv;
				rr += tv + dv;
				rr += fv + bv;
				rr -= 6*center;

				res(i,j,k) = rhs - rr;
				
			}
		}
	});
}
template<class T>
void MultiGridSolver3D<T>::m_ExactSolve()
{
	m_Red_Black_Gauss(systemk[m_max_level-1], *(bk[m_max_level-1]), *(xk[m_max_level-1]),*(xk_new[m_max_level-1]), 100,1.0);
}

template<class T>
void MultiGridSolver3D<T>::m_Restrict(uniform_grid_descriptor_3D &next_level, uniform_grid_descriptor_3D &curr_level, Buffer3D<T> & src_level, Buffer3D<T> & dst_level)
{
	//int computational_elements = next_level.system_size;
	//
	//uint slice = next_level.gridx*next_level.gridy;
	//uint slice_fine = curr_level.gridx*curr_level.gridy;
	//uint ystride = next_level.gridx;
	//uint ystride_fine = curr_level.gridx;

	int computational_elements = dst_level._blockx*dst_level._blocky*dst_level._blockz;
	uint slice = dst_level._blockx*dst_level._blocky;


	//red sweep
	tbb::parallel_for(0, computational_elements, 1, [&dst_level,&slice,&src_level,&next_level](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/dst_level._blockx;
		uint bi = thread_idx%(dst_level._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<next_level.gridx && j<next_level.gridy && k<next_level.gridz)
			{

				uint i_fine = i*2, j_fine = j*2, k_fine = k*2;


				//uint tid2i2j2k = k_fine*slice_fine + j_fine*ystride_fine + i_fine;
				//uint tid2i2j_12k = tid2i2j2k + ystride_fine;
				//uint tid2i2j2k_1 = tid2i2j2k + slice_fine;
				//uint tid2i2j_12k_1 = tid2i2j2k_1 + ystride_fine;

				T r = 0;
				r += src_level(i_fine, j_fine, k_fine);
				r += src_level(i_fine+1, j_fine, k_fine);
				r += src_level(i_fine, j_fine+1, k_fine);
				r += src_level(i_fine+1, j_fine+1, k_fine);
				r += src_level(i_fine, j_fine, k_fine+1);
				r += src_level(i_fine+1, j_fine, k_fine+1);
				r += src_level(i_fine, j_fine+1, k_fine+1);
				r += src_level(i_fine+1, j_fine+1, k_fine+1);

				r*=0.5;

				dst_level(i,j,k) = r;
			}
		}

	});
}
template<class T>
void MultiGridSolver3D<T>::m_MonopleToMonopole(uniform_grid_descriptor_3D &next_level, uniform_grid_descriptor_3D &curr_level, Buffer3D<T> & src_level, Buffer3D<T> & dst_level)
{
	//int computational_elements = next_level.system_size;
	//
	//uint slice = next_level.gridx*next_level.gridy;
	//uint slice_fine = curr_level.gridx*curr_level.gridy;
	//uint ystride = next_level.gridx;
	//uint ystride_fine = curr_level.gridx;

	int computational_elements = dst_level._blockx*dst_level._blocky*dst_level._blockz;
	uint slice = dst_level._blockx*dst_level._blocky;


	//red sweep
	tbb::parallel_for(0, computational_elements, 1, [&dst_level,&slice,&src_level,&next_level](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/dst_level._blockx;
		uint bi = thread_idx%(dst_level._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<next_level.gridx && j<next_level.gridy && k<next_level.gridz)
			{

				uint i_fine = i*2, j_fine = j*2, k_fine = k*2;


				//uint tid2i2j2k = k_fine*slice_fine + j_fine*ystride_fine + i_fine;
				//uint tid2i2j_12k = tid2i2j2k + ystride_fine;
				//uint tid2i2j2k_1 = tid2i2j2k + slice_fine;
				//uint tid2i2j_12k_1 = tid2i2j2k_1 + ystride_fine;

				T r = 0;
				r += src_level(i_fine, j_fine, k_fine);
				r += src_level(i_fine+1, j_fine, k_fine);
				r += src_level(i_fine, j_fine+1, k_fine);
				r += src_level(i_fine+1, j_fine+1, k_fine);
				r += src_level(i_fine, j_fine, k_fine+1);
				r += src_level(i_fine+1, j_fine, k_fine+1);
				r += src_level(i_fine, j_fine+1, k_fine+1);
				r += src_level(i_fine+1, j_fine+1, k_fine+1);

				dst_level(i,j,k) = r;
			}
		}

	});
}

template<class T>
void MultiGridSolver3D<T>::m_Prolongate(uniform_grid_descriptor_3D &next_level, uniform_grid_descriptor_3D &curr_level, Buffer3D<T> & src_level, Buffer3D<T> & dst_level)
{
	int computational_elements = dst_level._blockx*dst_level._blocky*dst_level._blockz;
	uint slice = dst_level._blockx*dst_level._blocky;


	//red sweep
	tbb::parallel_for(0, computational_elements, 1, [&dst_level,&slice,&next_level,&curr_level,&src_level](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/dst_level._blockx;
		uint bi = thread_idx%(dst_level._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<next_level.gridx && j<next_level.gridy && k<next_level.gridz)
			{
				//double h = 1.0/((float)(next_level.gridx));
				//float wx = ((float)i+0.5) * h,wy = ((float)j+0.5) * h,wz = ((float)k+0.5) * h;
				uint i_coarse = i/2, j_coarse = j/2, k_coarse = k/2;
				T r = src_level(i_coarse,j_coarse,k_coarse);
				//T r = src_level.sample_cubic(wx,wy,wz);
				dst_level(i,j,k) += r;
			}
		}
	});
}

template<class T>
void MultiGridSolver3D<T>::m_Vcycle(Buffer3D<T> * x, Buffer3D<T> * b, T tol, T &residual, int level)
{
	xk[level] = x;
	bk[level] = b;

	for(int i=level+1; i<m_max_level; i++)
	{
		xk[i]->setZero();
		bk[i]->setZero();
		rk[i]->setZero();
	}
	rk[level]->setZero();
	

	for(int i=level; i<m_max_level-1; ++i)
	{
		m_Red_Black_Gauss(systemk[i],
			*(bk[i]),
			*(xk[i]),
			*(xk_new[i]),
			5,
			1.0);

		m_ComputeResidual(systemk[i],
			*(rk[i]),
			*(bk[i]),
			*(xk[i]));

		m_Restrict(systemk[i+1], 
			systemk[i],
			*(rk[i]),
			*(bk[i+1]));
	}
	m_ExactSolve();
	for(int i=m_max_level-2; i>=level; --i)
	{
		m_Prolongate(systemk[i], 
			systemk[i+1],
			*(xk[i+1]),
			*(xk[i]));
		
		m_Red_Black_Gauss(systemk[i],
			*(bk[i]),
			*(xk[i]),
			*(xk_new[i]),
			5,
			1.0);
		
		
	}
	/*m_Red_Black_Gauss(systemk[0],
		*(bk[0]),
		*(xk[0]),
		*(xk_new[0]),
		2,
		1.0);*/
	for(int i=level+1; i<m_max_level; i++)
	{
		xk[i]->setZero();
		bk[i]->setZero();
		rk[i]->setZero();

	}
}

template<class T>
void MultiGridSolver3D<T>::m_FullMultiGrid(Buffer3D<T> * x, Buffer3D<T> * b, T tol, T &residual)
{
	xk[0] = x;
	bk[0] = b;
	T res;
	for(int i=1; i<m_max_level; i++)
	{
		xk[i]->setZero();
		bk[i]->setZero();
		rk[i]->setZero();
	}
	rk[0]->setZero();

	for(int i=0; i<m_max_level-1; ++i)
	{


		m_ComputeResidual(systemk[i],
			*(rk[i]),
			*(bk[i]),
			*(xk[i]));

		m_Restrict(systemk[i+1], 
			systemk[i],
			*(rk[i]),
			*(bk[i+1]));
	}
	m_ExactSolve();
	for(int i=m_max_level-2; i>=0; --i)
	{
		m_Prolongate(systemk[i], 
			systemk[i+1],
			*(xk[i+1]),
			*(xk[i]));

		m_Vcycle(xk[i], bk[i], 1e-10, res, i);
	}
	for(int i=0; i<1; i++)
	{
		m_Vcycle(xk[0],b,1e-10, res,0);
	}
}

typedef MultiGridSolver3D<float> MGSolverf;
typedef MultiGridSolver3D<double> MGSolverd;

}

#endif