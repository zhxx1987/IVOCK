#ifndef _fluid_particle_h_
#define _fluid_particle_h_
#include "array.h"
#include "fluid_buffer3D.h"
#include <vector>
#include <random>
#include "tbb/tbb.h"
using namespace std;




struct particle 
{
	float pos[3];
	float vel[3];
	float dens,temp;//other quantities
	int hash_idx;
	int sub_hash_idx;
	particle(){
		pos[0] = 0;
		pos[1] = 0;
		pos[2] = 0;
		vel[0] = 0;
		vel[1] = 0;
		vel[2] = 0;
		dens=0;
		temp=0;
		hash_idx = 0;
		sub_hash_idx=0;
	}
	~particle(){}
	
	particle(float px, float py, float pz, 
		float vx, float vy, float vz,
		float rho, float heat)
	{
		pos[0] = px;
		pos[1] = py;
		pos[2] = pz;
		
		vel[0] = vx;
		vel[1] = vy;
		vel[2] = vz;
		dens=rho;
		temp=heat;
		hash_idx = 0;
		sub_hash_idx = 0;
	}
	
	particle(const particle & p)
	{
		for (int i=0;i<3;i++)
		{
			pos[i]=p.pos[i];
			vel[i] = p.vel[i];

		}
		dens=p.dens;
		temp=p.temp;
		hash_idx = p.hash_idx;
		sub_hash_idx = p.sub_hash_idx;
	}
};

class fluid_particle
{
public:
	int _nx,_ny,_nz;
	float _h;
	vector<particle> particles;
	buffer3Df coef;
	vector<vector<particle>> particle_hash;
	fluid_particle(){}
	~fluid_particle(){}
	float frand(float a, float b){
		float w = (float)(rand()%RAND_MAX)/(float)RAND_MAX;
		return a + (b-a)*w;
	}
	void init(int nx, int ny, int nz, float h)
	{
		_nx = nx;
		_ny = ny;
		_nz = nz;
		coef.init(nx+1,ny+1,nz+1,h,0,0,0);
		_h = h;
		particles.resize(0);
		particle_hash.resize(_nx*_ny*_nz);
		for (int k=0;k<nz;k++)for(int j=0;j<ny;j++)for(int i=0;i<nx;i++)
		{
			float x = (float)i*_h;
			float y = (float)j*_h;
			float z = (float)k*_h;
			for (int p=0;p<6;p++)
			{
				particles.push_back(particle(x+frand(-0.49*_h, 0.49*_h), 
					y+frand(-0.49*_h,0.49*_h), 
					z+frand(-0.49*_h,0.49*_h),
					0,0,0,0,0));
			}
			//for (int pk = -1; pk<=1;pk=pk+2)
			//	for (int pj = -1; pj<=1;pj=pj+2)
			//		for (int pi = -1; pi<=1;pi=pi+2)
			//{
			//	float px = x + (float)pi * 0.25 *_h;
			//	float py = y + (float)pj * 0.25 *_h;
			//	float pz = z + (float)pk * 0.25 *_h;
			//	particles.push_back(particle(px,py,pz,0,0,0,0,0));

			//}
		}
		sort();
	}

	void move_particles(float dt, buffer3Df &un, buffer3Df &vn, buffer3Df &wn)
	{
		float c1 = 2.0/9.0*dt, c2 = 3.0/9.0 * dt, c3 = 4.0/9.0 * dt;
		int num = particles.size();
		tbb::parallel_for(0, num, 1, [&](int thread_idx) 
		{
			float px = particles[thread_idx].pos[0], 
				py = particles[thread_idx].pos[1], 
				pz = particles[thread_idx].pos[2];
			float u1=un.sample_linear(px,py,pz);
			float v1=vn.sample_linear(px,py,pz);
			float w1=wn.sample_linear(px,py,pz);


			float midx = px + 0.5*dt*u1, 
				midy = py + 0.5*dt*v1, 
				midz = pz + 0.5*dt*w1;
			float u2 = un.sample_linear(midx,midy, midz);
			float v2 = vn.sample_linear(midx,midy, midz);
			float w2 = wn.sample_linear(midx,midy, midz);
			float midx2 = px + 0.75*dt*u2;
			float midy2 = py + 0.75*dt*v2;
			float midz2 = pz + 0.75*dt*w2;
			float u3 = un.sample_linear(midx2,midy2,midz2);
			float v3 = vn.sample_linear(midx2,midy2,midz2);
			float w3 = wn.sample_linear(midx2,midy2,midz2);

			px = px + c1 * u1 + c2 * u2 + c3*u3;
			py = py + c1 * v1 + c2 * v2 + c3*v3;
			pz = pz + c1 * w1 + c2 * w2 + c3*w3;
			px = min(max(0,px),(float)(_nx-1)*_h);
			py = min(max(0,py),(float)(_ny-1)*_h);
			pz = min(max(0,pz),(float)(_nz-1)*_h);
			particles[thread_idx].pos[0] = px; 
			particles[thread_idx].pos[1] = py; 
			particles[thread_idx].pos[2] = pz;
		});
	}

	float H(float r)
	{
		float res = 0;
		if(r>=-1 && r<0) res = 1+r;
		if(r>=0 && r<1) res = 1-r;
		return res;
	}
	float compute_weight(float gx,float gy,float gz,
		                 float px,float py,float pz)
	{
		//k(x,y,z) = H(dx/hx)H(dy/hx)H(dz/hx)
		//H(r) = 1-r 0<=r<1  1+r -1<=r<0 0 else;
		float dx = px - gx;
		float dy = py - gy;
		float dz = pz - gz;
		return H(dx/_h)*H(dy/_h)*H(dz/_h);
	}
	
	void particle_to_grid(float dt, buffer3Df & field, buffer3Df & un,buffer3Df & vn,buffer3Df & wn,int component)
	{
		coef.setZero();
		int compute_elements = field._blockx*field._blocky*field._blockz;

		int slice = field._blockx*field._blocky;

		tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

			uint bk = thread_idx/slice;
			uint bj = (thread_idx%slice)/field._blockx;
			uint bi = thread_idx%(field._blockx);

			for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
			{
				uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
				if(i<field._nx&&j<field._ny&&k<field._nz)
				{
					float world_x = ((float)i-field._ox)*_h;
					float world_y = ((float)j-field._oy)*_h;
					float world_z = ((float)k-field._oz)*_h;
					float total_w = 1e-4;
					float total_val = 0;
					for (int hkk=-1; hkk<=1;hkk++)
						for(int hjj=-1;hjj<=1;hjj++)
							for(int hii=-1;hii<=1;hii++)
							{
								int iii = (int)i+hii, jjj = (int)j+hjj, kkk = (int)k+hkk;
								if(iii>=0&&iii<_nx &&jjj>=0 && jjj<_ny && kkk>=0 && kkk<_nz)
								{

									int hash_idx = kkk*_nx*_ny + jjj*_nx + iii;

									for (int p=0; p<particle_hash[hash_idx].size();p++)
									{
										float weight = compute_weight(world_x,world_y,world_z,
											particle_hash[hash_idx][p].pos[0],
											particle_hash[hash_idx][p].pos[1],
											particle_hash[hash_idx][p].pos[2]);
										total_w += weight;
										if(component<3)
											total_val += weight*(particle_hash[hash_idx][p].vel[component]);
										if(component==3)
											total_val += weight*(particle_hash[hash_idx][p].dens);
										if(component==4)
											total_val += weight*(particle_hash[hash_idx][p].temp);

									}
								}
							}
							coef(i,j,k) = total_w;
							field(i,j,k) = total_val/total_w;
				}

			}
		});



		tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

			uint bk = thread_idx/slice;
			uint bj = (thread_idx%slice)/field._blockx;
			uint bi = thread_idx%(field._blockx);

			for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
			{
				uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
				if(i>0&&i<field._nx-1&&j>0&&j<field._ny-1&&k>0&&k<field._nz-1)
				{
					if (coef(i,j,k)<1e-3)
					{
						field(i,j,k) = (field(i-1,j,k) + field(i+1,j,k)
							+field(i,j-1,k) + field(i,j+1,k)
							+field(i,j,k-1) + field(i,j,k+1))/6.0;

					}
				}

			}
		});



	}
	void write_to_grid(float dt, buffer3Df &unp1, buffer3Df &vnp1, buffer3Df &wnp1, 
		buffer3Df &rho, buffer3Df &heat,buffer3Df &un, buffer3Df &vn, buffer3Df &wn)
	{
		unp1.setZero();
		vnp1.setZero();
		wnp1.setZero();
		rho.setZero();
		heat.setZero();
		particle_to_grid(dt,unp1,un,vn,wn,0);
		particle_to_grid(dt,vnp1,un,vn,wn,1);
		particle_to_grid(dt,wnp1,un,vn,wn,2);
		particle_to_grid(dt,rho, un,vn,wn,3);
		particle_to_grid(dt,heat,un,vn,wn,4);
		
		
	}
	void getDelta(buffer3Df &field_old, buffer3Df &field_new, buffer3Df &fieldtemp)
	{
		int compute_elements = fieldtemp._blockx*fieldtemp._blocky*fieldtemp._blockz;

		int slice = fieldtemp._blockx*fieldtemp._blocky;

		tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

			uint bk = thread_idx/slice;
			uint bj = (thread_idx%slice)/fieldtemp._blockx;
			uint bi = thread_idx%(fieldtemp._blockx);

			for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
			{
				uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
				if( i<fieldtemp._nx && j<fieldtemp._ny && k<fieldtemp._nz )
				{
					fieldtemp(i,j,k) = field_new(i,j,k) - field_old(i,j,k);
				}
			}
		});
	}

	void interpolate(buffer3Df &uold, buffer3Df &vold, buffer3Df &wold, 
		buffer3Df &unew, buffer3Df &vnew, buffer3Df &wnew,
		buffer3Df &utemp, buffer3Df &vtemp, buffer3Df &wtemp)
	{
		utemp.setZero();
		vtemp.setZero();
		wtemp.setZero();
		getDelta(uold,unew,utemp);
		getDelta(vold,vnew,vtemp);
		getDelta(wold,wnew,wtemp);

		int num = particles.size();
		tbb::parallel_for(0, num, 1, [&](int thread_idx) 
		{
			float px = particles[thread_idx].pos[0], 
				py = particles[thread_idx].pos[1], 
				pz = particles[thread_idx].pos[2];

			particles[thread_idx].vel[0] = (
				0.03 * unew.sample_linear(px,py,pz) 
								+ 0.97 *(particles[thread_idx].vel[0] + utemp.sample_linear(px,py,pz))
				);

			particles[thread_idx].vel[1] = (
				0.03 * vnew.sample_linear(px,py,pz) 
								+ 0.97 *(particles[thread_idx].vel[1] + vtemp.sample_linear(px,py,pz))
				);

			particles[thread_idx].vel[2] = (
				0.03 * wnew.sample_linear(px,py,pz) 
								+ 0.97 *(particles[thread_idx].vel[2] + wtemp.sample_linear(px,py,pz))
				);

			//particles[thread_idx].dens = rho.sample_linear(px,py,pz);
			//particles[thread_idx].temp = heat.sample_linear(px,py,pz);

		});


	}
	void FLIP(float dt, buffer3Df &u, buffer3Df &v, buffer3Df &w, 
						buffer3Df &uold, buffer3Df &vold, buffer3Df &wold,
						buffer3Df &utemp, buffer3Df &vtemp, buffer3Df &wtemp,
						buffer3Df &rho, buffer3Df &heat)
	{
		resample(uold,vold,wold,rho,heat);
		interpolate(uold,vold,wold,u,v,w,utemp,vtemp,wtemp);
		move_particles(dt,u,v,w);//check
		sort();

		write_to_grid(dt,utemp,vtemp,wtemp, rho, heat,u,v,w);

	}
	void genHeat(float dt, buffer3Dc & h_desc, float smoke_dens, float smoke_heat)
	{
		int num = particles.size();
		tbb::parallel_for(0, num, 1, [&](int thread_idx) 
		{
			float px = particles[thread_idx].pos[0], 
				py = particles[thread_idx].pos[1], 
				pz = particles[thread_idx].pos[2];
			int grid_i = (int)(px/_h + 0.5);
			int grid_j = (int)(py/_h + 0.5);
			int grid_k = (int)(pz/_h + 0.5);
			if (h_desc(grid_i,grid_j,grid_k)==1)
			{
				particles[thread_idx].dens = smoke_dens;
				particles[thread_idx].temp = smoke_heat;
				int hash_idx = particles[thread_idx].hash_idx;
				int sub_idx  = particles[thread_idx].sub_hash_idx;
				particle_hash[hash_idx][sub_idx].dens = smoke_dens;
				particle_hash[hash_idx][sub_idx].temp = smoke_heat;
			}
		});
		//for (int i=0;i<particles.size();i++)
		//{
		//	if(particles[i].dens>0)
		//	{
		//		printf("%f ",particles[i].dens);
		//	}
		//}
		//printf("\nafter gen_heat\n");
		//sort();

	}
	void heat_decay(float dt, float temp_decay)
	{
		int num = particles.size();
		tbb::parallel_for(0, num, 1, [&](int thread_idx) 
		{
			
			particles[thread_idx].temp = particles[thread_idx].temp /(1.0 + temp_decay*dt);
			//particles[thread_idx].dens = particles[thread_idx].dens /(1.0 + 0.05*dt);
			
			if (particles[thread_idx].pos[1]/_h >= 0.9*(float)_ny)
			{
				particles[thread_idx].temp = 0;
				particles[thread_idx].dens = 0;
			}
			int hash_idx = particles[thread_idx].hash_idx;
			int sub_idx  = particles[thread_idx].sub_hash_idx;
			particle_hash[hash_idx][sub_idx].dens = particles[thread_idx].dens;
			particle_hash[hash_idx][sub_idx].temp = particles[thread_idx].temp;
		});
	}
	void sort()
	{
		//particle_hash.clear();
		//particle_hash.resize(_nx*_ny*_nz);
		for (int i=0;i<particle_hash.size();i++)
		{
			particle_hash[i].resize(0);
		}
		int num = particles.size();
		tbb::parallel_for(0, num, 1, [&](int thread_idx) 
		{
			float px = particles[thread_idx].pos[0], 
				py = particles[thread_idx].pos[1], 
				pz = particles[thread_idx].pos[2];
			int grid_i = (int)(px/_h + 0.5);
			int grid_j = (int)(py/_h + 0.5);
			int grid_k = (int)(pz/_h + 0.5);
			particles[thread_idx].hash_idx = (grid_k * _ny + grid_j)*_nx + grid_i;

		});
		for (int i=0;i<particles.size();i++)
		{
			particle_hash[particles[i].hash_idx].push_back(particles[i]);
			particles[i].sub_hash_idx = particle_hash[particles[i].hash_idx].size() - 1;
		}
	}
	void resample(buffer3Df &un, buffer3Df &vn, buffer3Df &wn, buffer3Df & rho, buffer3Df &heat)
	{
		
		//include reseeding, merge, delete
		vector<particle> new_particles;
		//copy particles
		for (int i=0; i<particle_hash.size();i++)
		{
			if(particle_hash[i].size()>2 && particle_hash[i].size()<8)
			{
				for (int p = 0 ; p<particle_hash[i].size(); p++)
				{
					new_particles.push_back(particle_hash[i][p]);
				}
			}
		}
		//reseed particles
		for (int k=0;k<_nz;k++)for(int j=0;j<_ny;j++)for(int i=0;i<_nx;i++)
		{
			if(particle_hash[ (k*_ny + j)*_nx + i ].size()<=2)//needs reseeding
			{
				
				float x = (float)i*_h;
				float y = (float)j*_h;
				float z = (float)k*_h;
				for (int p=0;p<5;p++)
				{
					float px = x + frand(-0.49*_h, 0.49*_h);
					float py = y + frand(-0.49*_h, 0.49*_h);
					float pz = z + frand(-0.49*_h, 0.49*_h);
					float u = un.sample_linear(px,py,pz);
					float v = vn.sample_linear(px,py,pz);
					float w = wn.sample_linear(px,py,pz);
					float dens = rho.sample_linear(px,py,pz);
					float temp = heat.sample_linear(px,py,pz);

					new_particles.push_back(particle(px, py, pz,u,v,w, dens,temp));
				}
				//for (int pk = -1; pk<=1;pk=pk+2)
				//	for(int pj = -1; pj<=1;pj=pj+2)
				//		for (int pi = -1; pi<=1;pi=pi+2)
				//{
				//	float px = x + (float)pi * 0.25 *_h;
				//	float py = y + (float)pj * 0.25 *_h;
				//	float pz = z + (float)pk * 0.25 *_h;
				//	float u = un.sample_linear(px,py,pz);
				//	float v = vn.sample_linear(px,py,pz);
				//	float w = wn.sample_linear(px,py,pz);
				//	float dens = rho.sample_linear(px,py,pz);
				//	float temp = heat.sample_linear(px,py,pz);
				//	new_particles.push_back(particle(px, py, pz,u,v,w, dens,temp));
				//}
			}
		}
		//merge particles
		for (int k=0;k<_nz;k++)for(int j=0;j<_ny;j++)for(int i=0;i<_nx;i++)
		{
			if(particle_hash[ (k*_ny + j)*_nx + i ].size()>=8)//needs merge
			{
				for (int p=0;p<6;p++)
				{
					new_particles.push_back(particle_hash[ (k*_ny + j)*_nx + i ][p]);
				}

				//float px = 0, py=0, pz = 0;
				//float pu = 0, pv=0, pw = 0;
				//float prho=0, ptemp=0;
				//int num_p = particle_hash[ (k*_ny + j)*_nx + i ].size()-8;
				//for (int p = 8; p <particle_hash[ (k*_ny + j)*_nx + i ].size() ; p++)
				//{
				//	px += particle_hash[ (k*_ny + j)*_nx + i ][p].pos[0];
				//	py += particle_hash[ (k*_ny + j)*_nx + i ][p].pos[1];
				//	pz += particle_hash[ (k*_ny + j)*_nx + i ][p].pos[2];

				//	pu += particle_hash[ (k*_ny + j)*_nx + i ][p].vel[0];
				//	pv += particle_hash[ (k*_ny + j)*_nx + i ][p].vel[1];
				//	pw += particle_hash[ (k*_ny + j)*_nx + i ][p].vel[2];

				//	prho += particle_hash[ (k*_ny + j)*_nx + i ][p].dens;
				//	ptemp += particle_hash[ (k*_ny + j)*_nx + i ][p].temp;



				//}
				//px=px/(float)num_p; py = py/(float)num_p; pz = pz/(float)num_p;
				//pu=pu/(float)num_p; pv = pv/(float)num_p; pw = pw/(float)num_p;
				//prho=prho/(float)num_p; ptemp=ptemp/(float)num_p;
				//new_particles.push_back(particle(px,py,pz,pu,pv,pw,prho,ptemp));

				
			}
		}
		particles.resize(0);
		for (int p=0;p<new_particles.size();p++)
		{
			particles.push_back(new_particles[p]);
		}
	}

};










#endif