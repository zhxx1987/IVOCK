#include "Smoke_solver3D.h"
#include <fstream> 

using namespace std;

void SmokeSolver3D::advect( float dt )
{
	_utemp.setZero();
	_vtemp.setZero();
	_wtemp.setZero();
	vort_confine_x.setZero();
	vort_confine_y.setZero();
	vort_confine_z.setZero();

	compute_curl();
	vorticity_confinement(dt*_vort_confine_coef);


	advect_field_cubic_clamp(dt,_un,_utemp);
	advect_field_cubic_clamp(dt,_vn,_vtemp);
	advect_field_cubic_clamp(dt,_wn,_wtemp);

	//tbb::parallel_for(0,3,1,[&](int i)
	//{
	//	paralle_advection(i, dt);
	//});



	_Ttemp.setZero();
	_Ttempnp1.setZero();
	BFECC_field(dt,_rho,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_rho,_Ttempnp1);
	_rho.copy(_Ttempnp1);


	_Ttemp.setZero();
	_Ttempnp1.setZero();
	BFECC_field(dt,_Tbf,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_Tbf,_Ttempnp1);
	_Tbf.copy(_Ttempnp1);

	_Ttemp.setZero();
	_Ttempnp1.setZero();
	BFECC_field(dt,_fuel,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_fuel,_Ttempnp1);
	_fuel.copy(_Ttempnp1);

	clampSmoke();


	_un.copy(_utemp);
	_vn.copy(_vtemp);
	_wn.copy(_wtemp);

	add_vort_confinement(_un, vort_confine_x);
	add_vort_confinement(_vn, vort_confine_y);
	add_vort_confinement(_wn, vort_confine_z);


}
void SmokeSolver3D::getCFL(float dt )
{
	float max_v = _hx;
	for (int k=0; k<_nz;k++)for (int j=0; j<_ny;j++)for (int i=0; i<_nx+1;i++)
	{
		if (fabs(_un(i,j,k))>max_v)
		{
			max_v = _un(i,j,k);
		}
	}
	for (int k=0; k<_nz;k++)for (int j=0; j<_ny+1;j++)for (int i=0; i<_nx;i++)
	{
		if (fabs(_vn(i,j,k))>max_v)
		{
			max_v = _vn(i,j,k);
		}
	}
	for (int k=0; k<_nz+1;k++)for (int j=0; j<_ny;j++)for (int i=0; i<_nx;i++)
	{
		if (fabs(_wn(i,j,k))>max_v)
		{
			max_v = _wn(i,j,k);
		}
	}
	_cfl = _hx/max_v;
	printf("cfl radio:%lf\n",dt/_cfl);

}
Vec3f SmokeSolver3D::get_velocity(Vec3f & pos)
{
	float u = _un.sample_linear(pos[0],pos[1],pos[2]);
	float v = _vn.sample_linear(pos[0],pos[1],pos[2]);
	float w = _wn.sample_linear(pos[0],pos[1],pos[2]);
	return Vec3f(u,v,w);
}
Vec3f SmokeSolver3D::traceRK3(Vec3f & pos, float dt)
{
	float c1 = 2.0/9.0*dt, c2 = 3.0/9.0 * dt, c3 = 4.0/9.0 * dt;
	Vec3f input = pos;
	Vec3f velocity1 = get_velocity(input);
	Vec3f midp1 = input + ((float)(0.5*dt))*velocity1;
	Vec3f velocity2 = get_velocity(midp1);
	Vec3f midp2 = input + ((float)(0.75*dt))*velocity2;
	Vec3f velocity3 = get_velocity(midp2);
	//velocity = get_velocity(input + 0.5f*dt*velocity);
	//input += dt*velocity;
	input = input + c1*velocity1 + c2*velocity2 + c3*velocity3;
	return input;
}
Vec3f SmokeSolver3D::trace(float dt, Vec3f & pos)
{
	if (dt>0)
	{
		float T=dt;
		Vec3f opos=pos;
		float t = 0;
		float substep = _cfl;
		while(t < T) {

			if(t + substep > T)
				substep = T - t;
			opos = traceRK3(opos,-substep);
			t+=substep;
		}
		//opos[0] = min(max(0.0f,opos[0]),(float)(_nx-1)*_hx );
		//opos[1] = min(max(0.0f,opos[1]),(float)(_ny-1)*_hx );
		//opos[2] = min(max(0.0f,opos[2]),(float)(_nz-1)*_hx );
		return opos;
	}
	else
	{
		float T = -dt;
		Vec3f opos=pos;
		float t = 0;
		float substep = _cfl;
		while(t < T) {

			if(t + substep > T)
				substep = T - t;
			opos = traceRK3(opos,substep);
			t+=substep;
		}
		//opos[0] = min(max(0.0f,opos[0]),(float)(_nx-1)*_hx );
		//opos[1] = min(max(0.0f,opos[1]),(float)(_ny-1)*_hx );
		//opos[2] = min(max(0.0f,opos[2]),(float)(_nz-1)*_hx );
		return opos;
	}

}
void SmokeSolver3D::time_step( float dt , int adv_type)
{
	getCFL(dt);
	float t = 0;
	printf("tracer step\n");
	if(_cfl<0.001) _cfl=0.001;
	while(t < dt) {
		float substep = 2.0*_cfl;   
		if(t + substep > dt)
			substep = dt - t;

		emit_tracers();
		advect_tracers(substep);
		t+=substep;
	}
	cout<<"emitting tracers done:"<<tracers.size()<<" tracers"<<endl;
	printf("tracer step done\n");
	switch (adv_type)
	{
	case 0:
		advect(dt);
		break;
	case 1:
		iVICK(dt);
		break;
	case 2:
		BFECC(dt);
		break;
	case 3:
		BFECC_IVOCK(dt);
		break;
	case 4:
		MacCormack(dt);
		break;
	case 5:
		MacCormack_IVOCK(dt);
		break;
	case 6:FLIP_advect(dt);
		break;
	case 7:FLIP_IVOCK(dt);
		break;
	case 8:Best_Combine(dt);
		break;
	default:
		break;
	}
	//clearBoundary();
	if(adv_type>=6&&adv_type<=7){


		FLIP_Solver.heat_decay(0.5*dt,_temp_decay);
		FLIP_Solver.genHeat(0.5*dt,_h_desc,_smoke_dens,_smoke_heat);

	}
	else{

		heat_decay(0.5*dt);
		//diffuse_heat(0.5*dt);
		gen_heat(0.5*dt);
		//reaction(dt);
	}
	int compute_num = _nx*_ny*_nz;
	int slice = _nx*_ny;
	tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/_nx;
		int i = thread_idx%_nx;
		if ( _b_desc(i,j,k)==1 )//empty
		{
			_un(i,j,k) = 0;
			_un(i+1,j,k) = 0;
			_vn(i,j,k) = 0;
			_vn(i,j+1,k) = 0;
			_wn(i,j,k) = 0;
			_wn(i,j,k+1) = 0;
		}
	});
	clearBoundary();
	add_force(dt);
	if(adv_type>=6&&adv_type<=7){


		FLIP_Solver.heat_decay(0.5*dt,_temp_decay);
		FLIP_Solver.genHeat(0.5*dt,_h_desc,_smoke_dens,_smoke_heat);

	}
	else{

		heat_decay(0.5*dt);
		//diffuse_heat(0.5*dt);
		gen_heat(0.5*dt);
	}
	//for (int i=0;i<1;i++)
	//{
	//	projection();
	//}
	pcg_projection();
	clearBoundary();
	//extrapolate(u_extrap,_Tbf,u_valid);
	//extrapolate(w_extrap,_fuel,w_valid);
	//set_bc();
}

void SmokeSolver3D::diffuse_heat(float dt)
{
	_Ttemp.setZero();
	_Ttemp.copy(_Tbf);

	int compute_elements = _rho._blockx*_rho._blocky*_rho._blockz;

	int slice = _rho._blockx*_rho._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		int bk = thread_idx/slice;
		int bj = (thread_idx%slice)/_rho._blockx;
		int bi = thread_idx%(_rho._blockx);

		for (int kk=0;kk<8;kk++)for(int jj=0;jj<8;jj++)for(int ii=0;ii<8;ii++)
		{
			int i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			_Ttemp(i,j,k) += 0.00001 * (
				_Tbf.at(i+1,j,k) + _Tbf.at(i-1,j,k) +
				_Tbf.at(i,j+1,k) + _Tbf.at(i,j-1,k) +
				_Tbf.at(i,j,k-1) + _Tbf.at(i,j,k+1) -
				6*_Tbf(i,j,k)
				);
		}
	});
	_Tbf.copy(_Ttemp);

}
void SmokeSolver3D::diffuse_buffer(float c, buffer3Df & diffuse_field)
{
	_Ttemp.setZero();
	int compute_elements = diffuse_field._blockx*diffuse_field._blocky*diffuse_field._blockz;

	int slice = diffuse_field._blockx*diffuse_field._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		int bk = thread_idx/slice;
		int bj = (thread_idx%slice)/diffuse_field._blockx;
		int bi = thread_idx%(diffuse_field._blockx);

		for (int kk=0;kk<8;kk++)for(int jj=0;jj<8;jj++)for(int ii=0;ii<8;ii++)
		{
			int i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			_Ttemp(i,j,k) = diffuse_field(i,j,k) + c * (
				diffuse_field.at(i+1,j,k) + diffuse_field.at(i-1,j,k) +
				diffuse_field.at(i,j+1,k) + diffuse_field.at(i,j-1,k) +
				diffuse_field.at(i,j,k-1) + diffuse_field.at(i,j,k+1) -
				6*diffuse_field(i,j,k)
				);
		}
	});
	diffuse_field.copy(_Ttemp);
	_Ttemp.setZero();
}
void SmokeSolver3D::heat_decay( float dt )
{
	int compute_elements = _rho._blockx*_rho._blocky*_rho._blockz;

	int slice = _rho._blockx*_rho._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_rho._blockx;
		uint bi = thread_idx%(_rho._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			float dT = _Tbf(i,j,k)/2000.0;
			float dT2 = dT*dT;
			float dT4 = dT2*dT2;
			//_Tbf(i,j,k) = _Tbf(i,j,k) - dt*_temp_decay*dT4;
			_Tbf(i,j,k) = _Tbf(i,j,k)*exp(-_temp_decay*dt);
			_rho(i,j,k) = _rho(i,j,k)*exp(-0.1*dt);
			//_fuel(i,j,k) = _fuel(i,j,k)/(1.0 + 1.0*dt);


		}
	});
}

void SmokeSolver3D::gen_heat( float dt )
{
	int compute_elements = _rho._blockx*_rho._blocky*_rho._blockz;

	int slice = _rho._blockx*_rho._blocky;
	int randam_force = rand();
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) 
		//for(int thread_idx=0; thread_idx<compute_elements; thread_idx++)
	{

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_rho._blockx;
		uint bi = thread_idx%(_rho._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(_h_desc(i,j,k)==1)//there is a source
			{
				_Tbf(i,j,k) = _smoke_heat;
				_rho(i,j,k) = _smoke_dens;
				//_fuel(i,j,k) = _smoke_dens;
				//_vn(i,j,k) = 3.0;
			}
			//if(_h_desc(i,j,k)==1&& (randam_force % 10) <=2 )//there is a source
			//{
			//	//_Tbf(i,j,k) = _smoke_heat;
			//	//_rho(i,j,k) = _smoke_dens;
			//	//_fuel(i,j,k) = _smoke_dens;
			//	float c = ((randam_force%10)>5)?1.0f:-1.0f;
			//	_vn(i,j,k) += dt*(randam_force % 10);

			//}
		}
	});
}

void SmokeSolver3D::reaction(float dt)
{
	//if fuel > 0 and temperature > ignitition:
	//r = rate*dt; // rate could depend on temperature-ignition or on fuel, or just a simple constant
	//burn = fuel*(1-1/(1+r)); // one backwards euler step to estimate how much fuel is burnt
	//fuel -= burn;
	//smoke += burn;
	//divergence = alpha*burn; // alpha controls the expansion due to combustion
	//temperature += beta*burn; // beta controls the heat release due to combustion
	//end
	float rate = 10.0;
	float T_ignition = 350.0f;
	float expansion_rate = 30.0;
	float combustion_rate = 200.0;

	int compute_elements = _rho._blockx*_rho._blocky*_rho._blockz;

	int slice = _rho._blockx*_rho._blocky;
	_burn_div.setZero();
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) 
		//for(int thread_idx=0; thread_idx<compute_elements; thread_idx++)
	{

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_rho._blockx;
		uint bi = thread_idx%(_rho._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i>0 && i<_nx-1 && j>0 &&j<_ny-1 && k>0&&k<_nz-1)
			{

				if(_fuel(i,j,k)>0 && _Tbf(i,j,k)>T_ignition)
				{
					//float r = rate * dt;
					//float burn = _fuel(i,j,k)*(1.0 - 1.0/(1.0+r));
					//burn = max(burn, 0.0f);
					float r = rate * dt;
					float burn = (_fuel(i,j,k) - r)>0?r:_fuel(i,j,k);
					_fuel(i,j,k) -= burn;
					_fuel(i,j,k) = max(_fuel(i,j,k),0.0f);
					_rho(i,j,k) += 2.5*burn;
					_burn_div(i,j,k) = expansion_rate*burn;
					if(_b_desc(i,j,k)!=0) _burn_div(i,j,k) = 0;
					_Tbf(i,j,k) += combustion_rate*burn;


				}
			}
		}
	});
}
void SmokeSolver3D::add_force( float dt )
{
	int compute_elements = _rho._blockx*_rho._blocky*_rho._blockz;

	int slice = _rho._blockx*_rho._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_rho._blockx;
		uint bi = thread_idx%(_rho._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx && j<_ny && k<_nz)
			{
				float density = _rho(i,j,k);
				float temperature = _Tbf(i,j,k);
				float f = -dt*_alpha*density + dt*_beta*temperature;

				_vn(i,j,k) += 0.5*f;
			}
		}
	});
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_rho._blockx;
		uint bi = thread_idx%(_rho._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx && j>0&&j<_ny && k<_nz)
			{
				float density = _rho(i,j,k);
				float temperature = _Tbf(i,j,k);
				float f = -dt*_alpha*density + dt*_beta*temperature;

				_vn(i,j+1,k) += 0.5*f;
			}
		}
	});
}

void SmokeSolver3D::projection()
{
//	//set_bc();
//	_Ttemp.setZero();
//	compute_rhs(_hx*_hx*_hx);
//	_Ttemp.copy(_div);
//	compute_rhs(-_hx*_hx);
//	float res;
//	//_ppe_solver.m_applyOpenBoundaryCondition(&_Ttemp,&_div,_hx);
//	_ppe_solver._b_Dirichlet = true;
//	_p.setZero();
//	for(int iter=0;iter<4;iter++){
//
//		_ppe_solver.m_FullMultiGrid(&_p,&_div,1e-5,res);
//	}
//	apply_grad();
//	//set_bc();
}

void SmokeSolver3D::set_bc()
{
	int compute_elements = _b_desc._blockx*_b_desc._blocky*_b_desc._blockz;
	float scale = 1.0/_hx * _hx * _hx;
	int slice = _b_desc._blockx*_b_desc._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_b_desc._blockx;
		uint bi = thread_idx%(_b_desc._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx && j<_ny && k<_nz)
			{

				if (_b_desc(i,j,k)==1)//solid boundary
				{
					_un(i,j,k) = 0;
					_un(i+1,j,k) = 0;
					_vn(i,j,k) = 0;
					_vn(i,j+1,k) = 0;
					_wn(i,j,k) = 0;
					_wn(i,j,k+1) = 0;

				}
				//other descriptions goes here
				if(_b_desc(i,j,k)==2)//open boundary in top j
				{
					//float val = 0.5*(_vn(i,j,k)+_vn(i,j+1,k));

					_vn(i,j+1,k) = _vn(i,j,k);
				}
				if(_b_desc(i,j,k)==3)//open boundary in bottom j
				{
					//float val = 0.5*(_vn(i,j,k)+_vn(i,j+1,k));

					_vn(i,j,k) = _vn(i,j+1,k);
				}
				if(_b_desc(i,j,k)==4)//open boundary in left i
				{
					//float val = 0.5*(_vn(i,j,k)+_vn(i,j+1,k));

					_un(i,j,k) = _un(i+1,j,k);
				}
				if(_b_desc(i,j,k)==5)//open boundary in front k
				{
					//float val = 0.5*(_vn(i,j,k)+_vn(i,j+1,k));

					_wn(i,j,k) = _wn(i,j,k+1);
				}
				if(_b_desc(i,j,k)==6)//open boundary in right i
				{
					//float val = 0.5*(_vn(i,j,k)+_vn(i,j+1,k));

					_un(i+1,j,k) = _un(i,j,k);
				}
				if(_b_desc(i,j,k)==7)//open boundary in back k
				{
					//float val = 0.5*(_vn(i,j,k)+_vn(i,j+1,k));

					_wn(i,j,k+1) = _wn(i,j,k);
				}
			}
		}
	});
}
void SmokeSolver3D::getDivergence()
{
	int compute_elements = _div._blockx*_div._blocky*_div._blockz;
	float scale = 1.0/_hx;
	int slice = _div._blockx*_div._blocky;
	_burn_div.setZero();
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_p._blockx;
		uint bi = thread_idx%(_p._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx && j<_ny && k<_nz)
			{
				_burn_div(i,j,k) = scale * (_un(i+1,j,k)-_un(i,j,k)+_vn(i,j+1,k)-_vn(i,j,k)+_wn(i,j,k+1)-_wn(i,j,k));
			}
		}
	});
}
void SmokeSolver3D::compute_rhs(float scale)
{
	int compute_elements = _div._blockx*_div._blocky*_div._blockz;
	//float scale = 1.0/_hx * _hx * _hx;
	int slice = _div._blockx*_div._blocky;
	_div.setZero();
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_p._blockx;
		uint bi = thread_idx%(_p._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx-1 && j<_ny-1 && k<_nz-1 && i>0 &&j>0 &&k>0)
			{

				_div(i,j,k) = -1.0/_hx*(_un(i+1,j,k)-_un(i,j,k)+_vn(i,j+1,k)-_vn(i,j,k)+_wn(i,j,k+1)-_wn(i,j,k));
				_div(i,j,k) += _burn_div(i,j,k);
				_div(i,j,k) = _div(i,j,k)*scale;
			}
		}
	});
}

void SmokeSolver3D::apply_grad()
{
	int compute_elements = _p._blockx*_p._blocky*_p._blockz;
	float scale = 1.0/_hx;
	int slice = _p._blockx*_p._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_p._blockx;
		uint bi = thread_idx%(_p._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx && j<_ny && k<_nz)
			{

				if(i>0)
					_un(i,j,k) -= scale*_p(i,j,k);
				if(j>0)
					_vn(i,j,k) -= scale*_p(i,j,k);
				if(k>0)
					_wn(i,j,k) -= scale*_p(i,j,k);
			}
		}
	});
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_p._blockx;
		uint bi = thread_idx%(_p._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx && j<_ny && k<_nz)
			{
				if(i<_nx-1)
					_un(i+1,j,k) += scale*_p(i,j,k);
				if(j<_ny-1)
					_vn(i,j+1,k) += scale*_p(i,j,k);
				if(k<_nz-1)
					_wn(i,j,k+1) += scale*_p(i,j,k);
			}
		}
	});
}

void SmokeSolver3D::advect_boundary(float dt, buffer3Df & field, buffer3Df &field_new)
{
	int compute_elements = field._blockx*field._blocky*field._blockz;

	int slice = field._blockx*field._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/field._blockx;
		uint bi = thread_idx%(field._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(!(i>10 && i<field._nx-10 && j>10 && j<field._ny-10 && k>10 && k<field._nz-10))
			{

				float world_x = ((float)i-field._ox)*_hx;
				float world_y = ((float)j-field._oy)*_hy;
				float world_z = ((float)k-field._oz)*_hz;
				Vec3f pos(world_x,world_y,world_z);
				//Vec3f trace_pos = trace(dt, pos);

				float u = _un.sample_linear(world_x,world_y,world_z);
				float v = _vn.sample_linear(world_x,world_y,world_z);
				float w = _wn.sample_linear(world_x,world_y,world_z);

				float px = world_x - 0.5*dt * u, py = world_y - 0.5*dt *v, pz = world_z - 0.5*dt*w;
				u = _un.sample_linear(px,py,pz);
				v = _vn.sample_linear(px,py,pz);
				w = _wn.sample_linear(px,py,pz);

				px = world_x - dt*u; py = world_y - dt*v; pz = world_z - dt*w;

				float SLv = field.sample_linear(px,py,pz);

				field_new(i,j,k) = SLv;

			}
		}
	});
}

void SmokeSolver3D::simple_advect_field(float dt, buffer3Df & field, buffer3Df &field_new)
{
	int compute_elements = field._blockx*field._blocky*field._blockz;

	int slice = field._blockx*field._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/field._blockx;
		uint bi = thread_idx%(field._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<field._nx && j<field._ny && k<field._nz)
			{

				float world_x = ((float)i-field._ox)*_hx;
				float world_y = ((float)j-field._oy)*_hy;
				float world_z = ((float)k-field._oz)*_hz;
				//Vec3f pos(world_x,world_y,world_z);
				//Vec3f trace_pos = trace(dt, pos);

				float u = _un.sample_linear(world_x,world_y,world_z);
				float v = _vn.sample_linear(world_x,world_y,world_z);
				float w = _wn.sample_linear(world_x,world_y,world_z);

				float px = world_x - 0.5*dt * u, 
					py = world_y - 0.5*dt *v, 
					pz = world_z - 0.5*dt*w;
				u = _un.sample_linear(px,py,pz);
				v = _vn.sample_linear(px,py,pz);
				w = _wn.sample_linear(px,py,pz);

				px = world_x - dt*u;
				py = world_y - dt*v;
				pz = world_z - dt*w;

				float SLv = field.sample_linear(px,py,pz);

				field_new(i,j,k) = SLv;

			}
		}
	});
}

void SmokeSolver3D::simple_advect_field_order1(float dt, buffer3Df & field, buffer3Df &field_new)
{
	int compute_elements = field._blockx*field._blocky*field._blockz;

	int slice = field._blockx*field._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/field._blockx;
		uint bi = thread_idx%(field._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<field._nx && j<field._ny && k<field._nz)
			{

				float world_x = ((float)i-field._ox)*_hx;
				float world_y = ((float)j-field._oy)*_hy;
				float world_z = ((float)k-field._oz)*_hz;
				//Vec3f pos(world_x,world_y,world_z);
				//Vec3f trace_pos = trace(dt, pos);

				float u = _un.sample_linear(world_x,world_y,world_z);
				float v = _vn.sample_linear(world_x,world_y,world_z);
				float w = _wn.sample_linear(world_x,world_y,world_z);



				float px = world_x - dt*u;
				float py = world_y - dt*v;
				float pz = world_z - dt*w;

				float SLv = field.sample_linear(px,py,pz);

				field_new(i,j,k) = SLv;

			}
		}
	});
}

void SmokeSolver3D::advect_field(float dt, buffer3Df & field, buffer3Df &field_new)
{
	int compute_elements = field._blockx*field._blocky*field._blockz;

	int slice = field._blockx*field._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/field._blockx;
		uint bi = thread_idx%(field._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<field._nx && j<field._ny && k<field._nz)
			{

				float world_x = ((float)i-field._ox)*_hx;
				float world_y = ((float)j-field._oy)*_hy;
				float world_z = ((float)k-field._oz)*_hz;
				Vec3f pos(world_x,world_y,world_z);
				//Vec3f trace_pos = trace(dt, pos);

				float u = _un.sample_linear(world_x,world_y,world_z);
				float v = _vn.sample_linear(world_x,world_y,world_z);
				float w = _wn.sample_linear(world_x,world_y,world_z);

				float px = world_x - 0.5*dt * u, py = world_y - 0.5*dt *v, pz = world_z - 0.5*dt*w;
				u = _un.sample_linear(px,py,pz);
				v = _vn.sample_linear(px,py,pz);
				w = _wn.sample_linear(px,py,pz);

				px = world_x - dt*u; py = world_y - dt*v; pz = world_z - dt*w;

				float SLv = field.sample_linear(px,
					py,
					pz);

				field_new(i,j,k) = SLv;

			}
		}
	});
}
void SmokeSolver3D::advect_field_cubic( float dt, buffer3Df & field, buffer3Df &field_new )
{
	int compute_elements = field._blockx*field._blocky*field._blockz;

	int slice = field._blockx*field._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/field._blockx;
		uint bi = thread_idx%(field._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<field._nx && j<field._ny && k<field._nz)
			{

				float world_x = ((float)i-field._ox)*_hx;
				float world_y = ((float)j-field._oy)*_hy;
				float world_z = ((float)k-field._oz)*_hz;
				Vec3f pos(world_x,world_y,world_z);
				Vec3f trace_pos = trace(dt, pos);

				//float u = _un.sample_linear(world_x,world_y,world_z);
				//float v = _vn.sample_linear(world_x,world_y,world_z);
				//float w = _wn.sample_linear(world_x,world_y,world_z);

				//float px = world_x - 0.5*dt * u, py = world_y - 0.5*dt *v, pz = world_z - 0.5*dt*w;
				//u = _un.sample_linear(px,py,pz);
				//v = _vn.sample_linear(px,py,pz);
				//w = _wn.sample_linear(px,py,pz);

				//px = world_x - dt*u; py = world_y - dt*v; pz = world_z - dt*w;

				float SLv = field.sample_cubic(trace_pos[0],trace_pos[1],trace_pos[2]);

				field_new(i,j,k) = SLv;



			}
		}
	});
}

void SmokeSolver3D::advect_field_cubic_clamp(float dt, buffer3Df & field, buffer3Df &field_new)
{
	int compute_elements = field._blockx*field._blocky*field._blockz;

	int slice = field._blockx*field._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/field._blockx;
		uint bi = thread_idx%(field._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<field._nx && j<field._ny && k<field._nz)
			{

				float world_x = ((float)i-field._ox)*_hx;
				float world_y = ((float)j-field._oy)*_hy;
				float world_z = ((float)k-field._oz)*_hz;
				Vec3f pos(world_x,world_y,world_z);
				//Vec3f trace_pos = trace(dt, pos);
				//
				float u = _un.sample_linear(world_x,world_y,world_z);
				float v = _vn.sample_linear(world_x,world_y,world_z);
				float w = _wn.sample_linear(world_x,world_y,world_z);

				float px = world_x - 0.5*dt * u, py = world_y - 0.5*dt *v, pz = world_z - 0.5*dt*w;
				u = _un.sample_linear(px,py,pz);
				v = _vn.sample_linear(px,py,pz);
				w = _wn.sample_linear(px,py,pz);

				px = world_x-dt*u, py = world_y-dt*v, pz = world_z-dt*w;
				float v0,v1,v2,v3,v4,v5,v6,v7;
				float SLvl = field.sample_cube_lerp(px,py,pz,
					v0,v1,v2,v3,v4,v5,v6,v7);

				float SLvc = field.sample_cubic(px,py,pz);

				float minv = min(min(min(min(min(min(min(v0,v1),v2),v3),v4),v5),v6),v7);
				float maxv = max(max(max(max(max(max(max(v0,v1),v2),v3),v4),v5),v6),v7);

				field_new(i,j,k) = SLvc;
				if( SLvc<=minv || SLvc >= maxv )
					field_new(i,j,k) = SLvl;

			}
		}
	});
}

//output to a pbrt scene....................
void SmokeSolver3D::output( uint nx, uint ny, uint nz, int frame, char* file_path)
{
	char file_name[256];
	sprintf(file_name,"%s/density_render.%04d.pbrt", file_path,frame);


	float *field = new float[nx*ny*nz];
	float *field_t=new float[nx*ny*nz];

	// normalize values
	//float max_dens = _rho(0,0,0);
	//for(int k=0;k<nz;k++)for(int j=0;j<ny;j++)for(int i=0;i<nx;i++) 
	//{
	//	if (_rho(0,0,0)>max_dens)
	//	{
	//		max_dens = _rho(i,j,k);
	//	}
	//}
	for(int k=0;k<nz;k++)for(int j=0;j<ny;j++)for(int i=0;i<nx;i++) 
	{
		if(_b_desc(i,j,k)==0){

			field[k*nx*ny + j*nx + i] = 0.5*max(_rho(i,j,k),0.0f)/10.0;
			field_t[k*nx*ny + j*nx + i] = _Tbf(i,j,k);
		}
		else
		{
			field[k*nx*ny + j*nx + i] =0;
			field_t[k*nx*ny + j*nx + i] = 0;

		}
	}


	std::fstream fout;
	fout.open(file_name, std::ios::out);

	int maxRes = (nx > ny) ? nx : ny;
	maxRes = (maxRes > nz) ? maxRes : nz;

	const float xSize = 1.0 / (float)maxRes * (float)nx;
	const float ySize = 1.0 / (float)maxRes * (float)ny;
	const float zSize = 1.0 / (float)maxRes * (float)nz;





	// dimensions
	fout<<"Volume \"volumegrid\" \n";
	fout<<" \"integer nx\""<<nx<<"\n";
	fout<<" \"integer ny\""<<ny<<"\n";
	fout<<" \"integer nz\""<<nz<<"\n";
	fout<<" \"point p0\" [ "<<-xSize/2.0<<" "<<0.0<<" "<<-zSize/2.0<<"] \"point p1\" "<<"[ "<<xSize/2.0<<" "<<ySize<<" "<<zSize/2.0<<" ]"<<"\n";
	fout<<" \"float density\" [ \n";
	for (int i = 0; i < nx*ny*nz; i++) 
		fout<<field[i]<<" ";
	fout<<"] \n";
	fout.close();

	//sprintf(file_name,"%s/density_render.%04d.bin", file_path,frame);
	//FILE *data_file = fopen(file_name,"wb");
	//fwrite(field,sizeof(float),nx*ny*nz,data_file);
	//fclose(data_file);

	//char file_name2[256];
	//sprintf(file_name2,"%s/density.%04d.bin", file_path,frame);
	//FILE *data_file = fopen(file_name2,"wb");
	//fwrite(field,sizeof(float),nx*ny*nz,data_file);
	//fclose(data_file);

	//char file_name3[256];
	//sprintf(file_name3,"%s/temperature.%04d.bin", file_path,frame);
	//data_file = fopen(file_name3,"wb");
	//fwrite(field_t,sizeof(float),nx*ny*nz,data_file);
	//fclose(data_file);




	delete[] field;
	delete[] field_t;
}
void SmokeSolver3D::mark_boundary()
{
	_boundary_marker.setZero();
	int compute_elements = _b_desc._blockx*_b_desc._blocky*_b_desc._blockz;

	int slice = _b_desc._blockx*_b_desc._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		int bk = thread_idx/slice;
		int bj = (thread_idx%slice)/_b_desc._blockx;
		int bi = thread_idx%(_b_desc._blockx);

		for (int kk=0;kk<8;kk++)for(int jj=0;jj<8;jj++)for(int ii=0;ii<8;ii++)
		{
			int i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx && j<_ny && k<_nz)
			{
				if(_b_desc(i,j,k)==2)//solid
				{
					for (int kkk =max(0,k-8); kkk<min((int)_nz-1,k+8); kkk++)
						for (int jjj =max(0,j-8); jjj<min((int)_ny-1,j+8); jjj++)
							for (int iii =max(0,i-8); iii<min((int)_nx-1,i+8); iii++)
							{
								_boundary_marker(iii,jjj,kkk) = 1;
							}
				}
			}
		}
	});
}
void SmokeSolver3D::vorticity_confinement(float coeff)
{
	//after we get vorticity;
	int compute_elements = vort_confine_z._blockx*vort_confine_z._blocky*vort_confine_z._blockz;
	float scale = 0.5/_hx;
	int slice = vort_confine_z._blockx*vort_confine_z._blocky;
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/vort_confine_z._blockx;
		uint bi = thread_idx%(vort_confine_z._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i>5 && i<_nx-5 && j>5&&j<_ny-5&&k>5&&k<_nz-5)
			{
				float world_x = ((float)i )*_hx;
				float world_y = ((float)j )*_hy;
				float world_z = ((float)k )*_hz;

				float Uwx = _wxn.sample_linear(world_x, world_y+_hx, world_z);
				float Uwy = _wyn.sample_linear(world_x, world_y+_hx, world_z);
				float Uwz = _wzn.sample_linear(world_x, world_y+_hx, world_z);

				float Uw = sqrt(Uwx*Uwx+Uwy*Uwy+Uwz*Uwz);

				float Dwx = _wxn.sample_linear(world_x, world_y-_hx, world_z);
				float Dwy = _wyn.sample_linear(world_x, world_y-_hx, world_z);
				float Dwz = _wzn.sample_linear(world_x, world_y-_hx, world_z);

				float Dw = sqrt(Dwx*Dwx+Dwy*Dwy+Dwz*Dwz);


				float Rwx = _wxn.sample_linear(world_x+_hx, world_y, world_z);
				float Rwy = _wyn.sample_linear(world_x+_hx, world_y, world_z);
				float Rwz = _wzn.sample_linear(world_x+_hx, world_y, world_z);

				float Rw = sqrt(Rwx*Rwx+Rwy*Rwy+Rwz*Rwz);

				float Lwx = _wxn.sample_linear(world_x-_hx, world_y, world_z);
				float Lwy = _wyn.sample_linear(world_x-_hx, world_y, world_z);
				float Lwz = _wzn.sample_linear(world_x-_hx, world_y, world_z);

				float Lw = sqrt(Lwx*Lwx+Lwy*Lwy+Lwz*Lwz);

				float Fwx = _wxn.sample_linear(world_x, world_y, world_z-_hx);
				float Fwy = _wyn.sample_linear(world_x, world_y, world_z-_hx);
				float Fwz = _wzn.sample_linear(world_x, world_y, world_z-_hx);

				float Fw = sqrt(Fwx*Fwx+Fwy*Fwy+Fwz*Fwz);

				float Bwx = _wxn.sample_linear(world_x, world_y, world_z+_hx);
				float Bwy = _wyn.sample_linear(world_x, world_y, world_z+_hx);
				float Bwz = _wzn.sample_linear(world_x, world_y, world_z+_hx);

				float Bw = sqrt(Bwx*Bwx+Bwy*Bwy+Bwz*Bwz);



				float Nx = scale*(Rw-Lw);
				float Ny = scale*(Uw-Dw);
				float Nz = scale*(Bw-Fw);

				float N = sqrt(Nx*Nx + Ny*Ny + Nz*Nz + 1e-5);
				Nx/=N;
				Ny/=N;
				Nz/=N;

				float Cwx = _wxn.sample_linear(world_x, world_y, world_z);
				float Cwy = _wyn.sample_linear(world_x, world_y, world_z);
				float Cwz = _wzn.sample_linear(world_x, world_y, world_z);

				vort_confine_x(i,j,k) = coeff*_hx*(Ny*Cwz-Nz*Cwy);
				vort_confine_y(i,j,k) = coeff*_hx*(Nz*Cwx-Nx*Cwz);
				vort_confine_z(i,j,k) = coeff*_hx*(Nx*Cwy-Ny*Cwx);

			}
		}
	});
}
void SmokeSolver3D::add_vort_confinement(buffer3Df & field, buffer3Df & vortf)
{
	int compute_elements = field._blockx*field._blocky*field._blockz;

	int slice = field._blockx*field._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/field._blockx;
		uint bi = thread_idx%(field._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<field._nx && j<field._ny && k<field._nz)
			{
				float world_x = ((float)i-field._ox)*_hx;
				float world_y = ((float)j-field._oy)*_hy;
				float world_z = ((float)k-field._oz)*_hz;

				float f = vortf.sample_linear(world_x,world_y,world_z);
				field(i,j,k) += f;
			}
		}
	});
}
void SmokeSolver3D::compute_curl()
{
	int compute_elements = _wxn._blockx*_wxn._blocky*_wxn._blockz;
	float scale = 1.0/_hx;
	int slice = _wxn._blockx*_wxn._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_wxn._blockx;
		uint bi = thread_idx%(_wxn._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			float y = ((float)j+0.5)/(float)_ny;
			if(i<_nx-1 && j<_ny-1 && k<_nz-1)
			{

				if(i>0 && j>0 && k>0)
				{
					_wxn(i,j,k) = scale * ( _wn(i,j,k)- _wn(i,j-1,k) - _vn(i,j,k) + _vn(i,j,k-1) );
					_wyn(i,j,k) = scale * ( _un(i,j,k) - _un(i,j,k-1) - _wn(i,j,k) + _wn(i-1,j,k));
					_wzn(i,j,k) = scale * ( _vn(i,j,k) - _vn(i-1,j,k) - _un(i,j,k) + _un(i,j-1,k));
				}
				else if(i==0)
				{
					_wxn(i,j,k) = _wxn(i+1,j,k);
					_wyn(i,j,k) = _wyn(i+1,j,k);
					_wzn(i,j,k) = _wzn(i+1,j,k);
				}
				else if(j==0)
				{
					_wxn(i,j,k) = _wxn(i,j+1,k);
					_wyn(i,j,k) = _wyn(i,j+1,k);
					_wzn(i,j,k) = _wzn(i,j+1,k);
				}
				else if(k==0)
				{
					_wxn(i,j,k) = _wxn(i,j,k+1);
					_wyn(i,j,k) = _wyn(i,j,k+1);
					_wzn(i,j,k) = _wzn(i,j,k+1);
				}
			}
			else if(i==_nx-1)
			{
				_wxn(i,j,k) = _wxn(i-1,j,k);
				_wyn(i,j,k) = _wyn(i-1,j,k);
				_wzn(i,j,k) = _wzn(i-1,j,k);
			}
			else if(j==_ny-1)
			{
				_wxn(i,j,k) = _wxn(i,j-1,k);
				_wyn(i,j,k) = _wyn(i,j-1,k);
				_wzn(i,j,k) = _wzn(i,j-1,k);
			}
			else if(k==_nz-1)
			{
				_wxn(i,j,k) = _wxn(i,j,k-1);
				_wyn(i,j,k) = _wyn(i,j,k-1);
				_wzn(i,j,k) = _wzn(i,j,k-1);
			}
		}
	});
}

void SmokeSolver3D::stretch( float dt )
{
	int compute_elements = _wxn._blockx*_wxn._blocky*_wxn._blockz;
	float l = sqrt(3.0)*0.5+ 0.01;
	l = 1.0;
	float scale = dt/(_hx*2.0*l);
	int slice = _wxn._blockx*_wxn._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_wxn._blockx;
		uint bi = thread_idx%(_wxn._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			float y = ((float)j+0.5)/(float)_ny;
			if(i<_nx-1 && j<_ny-1 && k<_nz-1)
			{

				if(i>1 && j>1 && k>1)
				{
					float world_x = ((float)i- _wxn._ox)*_hx;
					float world_y = ((float)j- _wxn._oy)*_hy;
					float world_z = ((float)k- _wxn._oz)*_hz;

					float wx = _wxn(i,j,k);
					float wy = _wyn.sample_linear(world_x,world_y,world_z);
					float wz = _wzn.sample_linear(world_x,world_y,world_z);

					float dx = 0;
					float dy = 0;
					float dz = 0;

					float L = sqrt(wx*wx + wy*wy + wz*wz + 1e-5);
					//if (L > 1e-5)
					{
						dx = wx/L;
						dy = wy/L;
						dz = wz/L;
					}

					float dwxdt = _un.sample_linear(world_x + l*_hx*dx, 
						world_y + l*_hx*dy, 
						world_z + l*_hx*dz) 
						- _un.sample_linear(world_x - l*_hx*dx, 
						world_y - l*_hx*dy, 
						world_z - l*_hx*dz);

					//float dudx = _un.sample_linear(world_x+0.5*_hx,world_y,world_z) - _un.sample_linear(world_x-0.5*_hx,world_y,world_z);
					//float dudy = _un.sample_linear(world_x,world_y+0.5*_hx,world_z) - _un.sample_linear(world_x,world_y-0.5*_hx,world_z);
					//float dudz = _un.sample_linear(world_x,world_y,world_z+0.5*_hx) - _un.sample_linear(world_x,world_y,world_z-0.5*_hx);

					_wxstr(i,j,k) = scale * dwxdt * L;

				}
			}
		}
	});

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_wyn._blockx;
		uint bi = thread_idx%(_wyn._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			float y = ((float)j+0.5)/(float)_ny;
			if(i<_nx-1 && j<_ny-1 && k<_nz-1)
			{

				if(i>1 && j>1 && k>1)
				{
					float world_x = ((float)i- _wyn._ox)*_hx;
					float world_y = ((float)j- _wyn._oy)*_hy;
					float world_z = ((float)k- _wyn._oz)*_hz;

					float wx = _wxn.sample_linear(world_x,world_y,world_z);
					float wy = _wyn.sample_linear(world_x,world_y,world_z);
					float wz = _wzn.sample_linear(world_x,world_y,world_z);

					float dx = 0;
					float dy = 0;
					float dz = 0;

					float L = sqrt(wx*wx + wy*wy + wz*wz+1e-5);
					//if (L > 1e-5)
					{
						dx = wx/L;
						dy = wy/L;
						dz = wz/L;
					}

					float dwydt = _vn.sample_linear(world_x + l*_hx*dx, 
						world_y + l*_hx*dy, 
						world_z + l*_hx*dz) 
						- _vn.sample_linear(world_x - l*_hx*dx, 
						world_y - l*_hx*dy, 
						world_z - l*_hx*dz);

					//float dvdx = _vn.sample_linear(world_x+0.5*_hx,world_y,world_z) - _vn.sample_linear(world_x-0.5*_hx,world_y,world_z);
					//float dvdy = _vn.sample_linear(world_x,world_y+0.5*_hx,world_z) - _vn.sample_linear(world_x,world_y-0.5*_hx,world_z);
					//float dvdz = _vn.sample_linear(world_x,world_y,world_z+0.5*_hx) - _vn.sample_linear(world_x,world_y,world_z-0.5*_hx);

					_wystr(i,j,k) = scale * dwydt * L;

				}
			}
		}
	});

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_wzn._blockx;
		uint bi = thread_idx%(_wzn._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			float y = ((float)j+0.5)/(float)_ny;
			if(i<_nx-1 && j<_ny-1 && k<_nz-1)
			{

				if(i>1 && j>1 && k>1)
				{
					float world_x = ((float)i- _wzn._ox)*_hx;
					float world_y = ((float)j- _wzn._oy)*_hy;
					float world_z = ((float)k- _wzn._oz)*_hz;

					float wx = _wxn.sample_linear(world_x,world_y,world_z);
					float wy = _wyn.sample_linear(world_x,world_y,world_z);
					float wz = _wzn.sample_linear(world_x,world_y,world_z);

					float dx = 0;
					float dy = 0;
					float dz = 0;

					float L = sqrt(wx*wx + wy*wy + wz*wz+1e-5);
					//if (L > 1e-5)
					{
						dx = wx/L;
						dy = wy/L;
						dz = wz/L;
					}

					float dwzdt = _wn.sample_linear(world_x + l*_hx*dx, 
						world_y + l*_hx*dy, 
						world_z + l*_hx*dz) 
						- _wn.sample_linear(world_x - l*_hx*dx, 
						world_y - l*_hx*dy, 
						world_z - l*_hx*dz);
					//float dwdx = _wn.sample_linear(world_x+0.5*_hx,world_y,world_z) - _wn.sample_linear(world_x-0.5*_hx,world_y,world_z);
					//float dwdy = _wn.sample_linear(world_x,world_y+0.5*_hx,world_z) - _wn.sample_linear(world_x,world_y-0.5*_hx,world_z);
					//float dwdz = _wn.sample_linear(world_x,world_y,world_z+0.5*_hx) - _wn.sample_linear(world_x,world_y,world_z-0.5*_hx);

					_wzstr(i,j,k) = scale * dwzdt * L;

				}
			}
		}
	});
	
	

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_wzn._blockx;
		uint bi = thread_idx%(_wzn._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;

			if(i<_nx-5 && j<_ny-5 && k<_nz-5 && i>5 && j>5 && k>5)
			{

				//if(i>5 && j>5 && k>5)
				//{
				_wxstr(i,j,k) += _wxn(i,j,k);
				_wystr(i,j,k) += _wyn(i,j,k);
				_wzstr(i,j,k) += _wzn(i,j,k);

				//}
			}
			else
			{
				_wxstr(i,j,k)= _wxn(i,j,k);
				_wystr(i,j,k)= _wyn(i,j,k);
				_wzstr(i,j,k)= _wzn(i,j,k);
			}
		}
	});
	
	for (int iter=0; iter<3;iter++)
	{
		float coef = 0.01*dt; 
		coef = min(0.1f,coef);
		diffuse_buffer(coef,_wxstr);
		diffuse_buffer(coef,_wystr);
		diffuse_buffer(coef,_wzstr);
	}

}


void SmokeSolver3D::decay_vortices( float dt, buffer3Df & wxnp1, buffer3Df &wynp1, buffer3Df &wznp1)
{
	int compute_elements = _wxn._blockx*_wxn._blocky*_wxn._blockz;
	float scale = 2.0*dt/_hx;
	int slice = _wxn._blockx*_wxn._blocky;
	u_valid.set_zero();
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_wxn._blockx;
		uint bi = thread_idx%(_wxn._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx-5 && j<_ny-5 && k<_nz-5)
			{

				if(i>5 && j>5 && k>5)
				{
					float world_x = ((float)i- _wxn._ox)*_hx;
					float world_y = ((float)j- _wxn._oy)*_hy;
					float world_z = ((float)k- _wxn._oz)*_hz;

					float divu = _burn_div.sample_linear(world_x,world_y,world_z);

					float coef = -dt*divu;
					coef = min(coef,0.0f);
					wxnp1(i,j,k) = wxnp1(i,j,k)*exp(coef);
					//u_valid(i,j,k) = 1;

				}
			}
		}
	});
	v_valid.set_zero();
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_wyn._blockx;
		uint bi = thread_idx%(_wyn._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx-5 && j<_ny-5 && k<_nz-5)
			{

				if(i>5 && j>5 && k>5)
				{
					float world_x = ((float)i- _wyn._ox)*_hx;
					float world_y = ((float)j- _wyn._oy)*_hy;
					float world_z = ((float)k- _wyn._oz)*_hz;

					float divu = _burn_div.sample_linear(world_x,world_y,world_z);

					float coef = -dt*divu;
					coef = min(coef,0.0f);
					wynp1(i,j,k) = wynp1(i,j,k)*exp(coef);
					//v_valid(i,j,k) = 1;
				}
			}
		}
	});
	w_valid.set_zero();
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_wzn._blockx;
		uint bi = thread_idx%(_wzn._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx-5 && j<_ny-5 && k<_nz-5)
			{

				if(i>5 && j>5 && k>5)
				{
					float world_x = ((float)i- _wzn._ox)*_hx;
					float world_y = ((float)j- _wzn._oy)*_hy;
					float world_z = ((float)k- _wzn._oz)*_hz;

					float divu = _burn_div.sample_linear(world_x,world_y,world_z);

					float coef = -dt*divu;
					coef = min(coef,0.0f);
					wznp1(i,j,k) = wznp1(i,j,k)*exp(coef);
					//w_valid(i,j,k) = 1;
				}
			}
		}
	});
}


void SmokeSolver3D::addCurlPsi()
{
	int compute_elements = _un._blockx*_un._blocky*_un._blockz;
	float scale = 1.0/_hx;
	int slice = _un._blockx*_un._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_un._blockx;
		uint bi = thread_idx%(_un._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx-4 && j<_ny-4 && k<_nz-4)
				//if(i<_vn._nx && i>0 && j<_vn._ny-1&&j>0 && k<_vn._nz-1&&k>0)
			{

				if(i>4 && j>4 && k>4)
				{
					_un(i,j,k) += scale * ( _Psiz.at(i,j+1,k)- _Psiz.at(i,j,k) - _Psiy.at(i,j,k+1) + _Psiy.at(i,j,k) );

				}
			}
		}
	});

	compute_elements = _vn._blockx*_vn._blocky*_vn._blockz;
	slice = _vn._blockx*_vn._blocky;
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_vn._blockx;
		uint bi = thread_idx%(_vn._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx-4 && j<_ny-4 && k<_nz-4)
				//if(i<_vn._nx && i>0 && j<_vn._ny-1&&j>0 && k<_vn._nz-1&&k>0)
			{

				if(i>4 && j>4 && k>4)
				{
					_vn(i,j,k) += scale * ( _Psix.at(i,j,k+1)- _Psix.at(i,j,k) - _Psiz.at(i+1,j,k) + _Psiz.at(i,j,k) );

				}
			}
		}
	});

	compute_elements = _wn._blockx*_wn._blocky*_wn._blockz;
	slice = _wn._blockx*_wn._blocky;
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_wn._blockx;
		uint bi = thread_idx%(_wn._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx-4 && j<_ny-4 && k<_nz-4)
				//if(i<_vn._nx && i>0 && j<_vn._ny-1&&j>0 && k<_vn._nz-1&&k>0)
			{

				if(i>4 && j>4 && k>4)
				{
					_wn(i,j,k) += scale * ( _Psiy.at(i+1,j,k)- _Psiy.at(i,j,k) - _Psix.at(i,j+1,k) + _Psix.at(i,j,k) );

				}
			}
		}
	});

}
void SmokeSolver3D::formRHS(float scale,float dt)
{
	_wxstr.setZero();
	_wystr.setZero();
	_wzstr.setZero();
	int compute_elements = _wxn._blockx*_wxn._blocky*_wxn._blockz;
	//float scale = _hx*_hx*_hx;
	int slice = _wxn._blockx*_wxn._blocky;
	int d = (int)(0.125*((float)_nx));
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx)
	{
		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_wxn._blockx;
		uint bi = thread_idx%(_wxn._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			float y = ((float)j+0.5)/(float)_ny;
			if(i<_nx-d && j<_ny-d && k<_nz-d && i>d &&j>d &&k>d&&/*y<=0.7&&*/_boundary_marker(i,j,k)==0)
			{
				_wxstr(i,j,k) = scale * (_wxnp1(i,j,k) - _wxn(i,j,k));
				_wystr(i,j,k) = scale * (_wynp1(i,j,k) - _wyn(i,j,k));
				_wzstr(i,j,k) = scale * (_wznp1(i,j,k) - _wzn(i,j,k));
			}
			else
			{
				_wxstr(i,j,k) = 0;
				_wystr(i,j,k) = 0;
				_wzstr(i,j,k) = 0;
			}
		}
	});

	

}

void SmokeSolver3D::iVICK( float dt )
{
	mark_boundary();

	_utemp.setZero();
	_vtemp.setZero();
	_wtemp.setZero();
	_unp1.setZero();
	_vnp1.setZero();
	_wnp1.setZero();
	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	_wxstr.setZero();
	_wystr.setZero();
	_wzstr.setZero();
	_wxnp1.setZero();
	_wynp1.setZero();
	_wznp1.setZero();


	vort_confine_x.setZero();
	vort_confine_y.setZero();
	vort_confine_z.setZero();

	compute_curl();
	vorticity_confinement(dt*_vort_confine_coef);

	advect_field_cubic_clamp(dt,_un,_unp1 );
	advect_field_cubic_clamp(dt,_vn,_vnp1 );
	advect_field_cubic_clamp(dt,_wn,_wnp1 );





	stretch(dt); //_wxstr,_wystr,_wzstr = stretch(w,u);
	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	getDivergence();
	decay_vortices(dt, _wxstr, _wystr, _wzstr);


	advect_field_cubic_clamp(dt,_wxstr,_wxnp1 );
	advect_field_cubic_clamp(dt,_wystr,_wynp1 );
	advect_field_cubic_clamp(dt,_wzstr,_wznp1 );
	


	_Ttemp.setZero();
	_Ttempnp1.setZero();
	advect_field_cubic_clamp(dt,_rho,_Ttempnp1);
	_rho.copy(_Ttempnp1);


	_Ttemp.setZero();
	_Ttempnp1.setZero();
	advect_field_cubic_clamp(dt,_Tbf,_Ttempnp1);
	_Tbf.copy(_Ttempnp1);

	_Ttemp.setZero();
	_Ttempnp1.setZero();
	advect_field_cubic_clamp(dt,_fuel,_Ttempnp1);
	_fuel.copy(_Ttempnp1);

	clampSmoke();

	_un.copy(_unp1);
	_vn.copy(_vnp1);
	_wn.copy(_wnp1);
	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	compute_curl();


	_Psix.setZero();
	_Psiy.setZero();
	_Psiz.setZero();
	formRHS(_hx*_hx*_hx,dt);
	_Psix.copy(_wxstr);
	_Psiy.copy(_wystr);
	_Psiz.copy(_wzstr);
	_wxstr.setZero();
	_wystr.setZero();
	_wzstr.setZero();
	formRHS(-_hx*_hx,dt);
//	float res;
//	_ppe_solver.m_applyOpenBoundaryCondition(&_Psix,&_wxstr,_hx);
//	_ppe_solver._b_Dirichlet = true;
//	_Psix.setZero();
//	_ppe_solver.m_Vcycle(&_Psix,&_wxstr, 1e-5,res,0);
//	_ppe_solver.m_applyOpenBoundaryCondition(&_Psiy,&_wystr,_hx);
//	_ppe_solver._b_Dirichlet = true;
//	_Psiy.setZero();
//	_ppe_solver.m_Vcycle(&_Psiy,&_wystr, 1e-5,res,0);
//	_ppe_solver.m_applyOpenBoundaryCondition(&_Psiz,&_wzstr,_hx);
//	_ppe_solver._b_Dirichlet = true;
//	_Psiz.setZero();
//	_ppe_solver.m_Vcycle(&_Psiz,&_wzstr, 1e-5,res,0);


	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psix,&_wxstr,_hx);
	solve_stream_Poisson(_Psix,_wxstr);
	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psiy,&_wystr,_hx);
	solve_stream_Poisson(_Psiy,_wystr);
	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psiz,&_wzstr,_hx);
	solve_stream_Poisson(_Psiz,_wzstr);




	addCurlPsi();
	add_vort_confinement(_un, vort_confine_x);
	add_vort_confinement(_vn, vort_confine_y);
	add_vort_confinement(_wn, vort_confine_z);
}

void SmokeSolver3D::BFECC_field( float dt, buffer3Df &field_n,buffer3Df &field_np1,buffer3Df &field_AUX )
{

	//advect field_n to field_AUX forward
	simple_advect_field(dt,field_n,field_AUX);

	//advect field_AUX to field_np1 backward
	simple_advect_field(-dt,field_AUX,field_np1);

	//field_AUX = field_n +0.5*(field_n - field_np1);
	int compute_elements = field_AUX._blockx*field_AUX._blocky*field_AUX._blockz;
	int slice = field_AUX._blockx*field_AUX._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx)
	{
		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/field_AUX._blockx;
		uint bi = thread_idx%(field_AUX._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<field_AUX._nx && j<field_AUX._ny && k<field_AUX._nz)
			{
				field_AUX(i,j,k) = field_n(i,j,k) + 0.5*(field_n(i,j,k)-field_np1(i,j,k));
			}
		}
	});
	field_np1.setZero();
	simple_advect_field(dt,field_AUX, field_np1);

}
void SmokeSolver3D::BFECC_field_order1( float dt, buffer3Df &field_n,buffer3Df &field_np1,buffer3Df &field_AUX )
{

	//advect field_n to field_AUX forward
	simple_advect_field_order1(dt,field_n,field_AUX);

	//advect field_AUX to field_np1 backward
	simple_advect_field_order1(-dt,field_AUX,field_np1);

	//field_AUX = field_n +0.5*(field_n - field_np1);
	int compute_elements = field_AUX._blockx*field_AUX._blocky*field_AUX._blockz;
	int slice = field_AUX._blockx*field_AUX._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx)
	{
		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/field_AUX._blockx;
		uint bi = thread_idx%(field_AUX._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<field_AUX._nx && j<field_AUX._ny && k<field_AUX._nz)
			{
				field_AUX(i,j,k) = field_n(i,j,k) + 0.5*(field_n(i,j,k)-field_np1(i,j,k));
			}
		}
	});
	field_np1.setZero();
	simple_advect_field_order1(dt,field_AUX, field_np1);

}
void SmokeSolver3D::BFECC( float dt )
{
	_utemp.setZero();
	_vtemp.setZero();
	_wtemp.setZero();
	_unp1.setZero();
	_vnp1.setZero();
	_wnp1.setZero();
	compute_curl();
	vorticity_confinement(dt*_vort_confine_coef);


	BFECC_field(dt,_un,_unp1,_utemp);
	BFECC_field(dt,_vn,_vnp1,_vtemp);
	BFECC_field(dt,_wn,_wnp1,_wtemp);
	clampExtrema(dt,_un,_unp1);
	clampExtrema(dt,_vn,_vnp1);
	clampExtrema(dt,_wn,_wnp1);





	_Ttemp.setZero();
	_Ttempnp1.setZero();
	BFECC_field(dt,_rho,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_rho,_Ttempnp1);
	_rho.copy(_Ttempnp1);


	_Ttemp.setZero();
	_Ttempnp1.setZero();
	BFECC_field(dt,_Tbf,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_Tbf,_Ttempnp1);
	_Tbf.copy(_Ttempnp1);

	_Ttemp.setZero();
	_Ttempnp1.setZero();
	BFECC_field(dt,_fuel,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_fuel,_Ttempnp1);
	_fuel.copy(_Ttempnp1);

	clampSmoke();


	_un.copy(_unp1);
	_vn.copy(_vnp1);
	_wn.copy(_wnp1);

	add_vort_confinement(_un, vort_confine_x);
	add_vort_confinement(_vn, vort_confine_y);
	add_vort_confinement(_wn, vort_confine_z);
}
void SmokeSolver3D::Best_Combine(float dt)
{
	mark_boundary();
	_utemp.setZero();
	_vtemp.setZero();
	_wtemp.setZero();
	_unp1.setZero();
	_vnp1.setZero();
	_wnp1.setZero();
	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	_wxstr.setZero();
	_wystr.setZero();
	_wzstr.setZero();
	_wxnp1.setZero();
	_wynp1.setZero();
	_wznp1.setZero();


	vort_confine_x.setZero();
	vort_confine_y.setZero();
	vort_confine_z.setZero();

	compute_curl();
	vorticity_confinement(dt*_vort_confine_coef);

	BFECC_field(dt,_un,_unp1,_utemp);
	BFECC_field(dt,_vn,_vnp1,_vtemp);
	BFECC_field(dt,_wn,_wnp1,_wtemp);
	clampExtrema(dt,_un,_unp1);
	clampExtrema(dt,_vn,_vnp1);
	clampExtrema(dt,_wn,_wnp1);



	stretch(dt); //_wxstr,_wystr,_wzstr = stretch(w,u);
	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	getDivergence();
	decay_vortices(dt, _wxstr, _wystr, _wzstr);
	advect_field_cubic_clamp(dt,_wxstr,_wxnp1);
	advect_field_cubic_clamp(dt,_wystr,_wynp1);
	advect_field_cubic_clamp(dt,_wzstr,_wznp1);
	for (int iter=0; iter<3;iter++)
	{
		float coef = 0.1*dt; 
		coef = min(0.1f,coef);
		diffuse_buffer(coef,_wxnp1);
		diffuse_buffer(coef,_wynp1);
		diffuse_buffer(coef,_wznp1);
	}
	//BFECC_field(dt,_wxstr,_wxnp1, _wxn);
	//BFECC_field(dt,_wystr,_wynp1, _wyn);
	//BFECC_field(dt,_wzstr,_wznp1, _wzn);
	//clampExtrema(dt, _wxstr,_wxnp1);
	//clampExtrema(dt, _wystr,_wynp1);
	//clampExtrema(dt, _wzstr,_wznp1);



	_Ttemp.setZero();
	_Ttempnp1.setZero();
	advect_field_cubic_clamp(dt,_rho,_Ttempnp1);
	_rho.copy(_Ttempnp1);


	_Ttemp.setZero();
	_Ttempnp1.setZero();
	advect_field_cubic_clamp(dt,_Tbf,_Ttempnp1);
	_Tbf.copy(_Ttempnp1);

	_Ttemp.setZero();
	_Ttempnp1.setZero();
	advect_field_cubic_clamp(dt,_fuel,_Ttempnp1);
	_fuel.copy(_Ttempnp1);

	clampSmoke();

	_un.copy(_unp1);
	_vn.copy(_vnp1);
	_wn.copy(_wnp1);
	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	compute_curl();


	_Psix.setZero();
	_Psiy.setZero();
	_Psiz.setZero();
	formRHS(0.95*_hx*_hx*_hx,dt);
	_Psix.copy(_wxstr);
	_Psiy.copy(_wystr);
	_Psiz.copy(_wzstr);
	_wxstr.setZero();
	_wystr.setZero();
	_wzstr.setZero();
	formRHS(-0.95*_hx*_hx,dt);
//	float res;
//	_ppe_solver.m_applyOpenBoundaryCondition(&_Psix,&_wxstr,_hx);
//	_ppe_solver._b_Dirichlet = true;
//	_Psix.setZero();
//	_ppe_solver.m_Vcycle(&_Psix,&_wxstr, 1e-5,res,0);
//
//	_ppe_solver.m_applyOpenBoundaryCondition(&_Psiy,&_wystr,_hx);
//	_ppe_solver._b_Dirichlet = true;
//	_Psiy.setZero();
//	_ppe_solver.m_Vcycle(&_Psiy,&_wystr, 1e-5,res,0);
//
//	_ppe_solver.m_applyOpenBoundaryCondition(&_Psiz,&_wzstr,_hx);
//	_ppe_solver._b_Dirichlet = true;
//	_Psiz.setZero();
//	_ppe_solver.m_Vcycle(&_Psiz,&_wzstr, 1e-5,res,0);

	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psix,&_wxstr,_hx);
	solve_stream_Poisson(_Psix,_wxstr);
	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psiy,&_wystr,_hx);
	solve_stream_Poisson(_Psiy,_wystr);
	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psiz,&_wzstr,_hx);
	solve_stream_Poisson(_Psiz,_wzstr);

	addCurlPsi();
	add_vort_confinement(_un, vort_confine_x);
	add_vort_confinement(_vn, vort_confine_y);
	add_vort_confinement(_wn, vort_confine_z);
}
void SmokeSolver3D::BFECC_IVOCK( float dt )
{
	mark_boundary();
	_utemp.setZero();
	_vtemp.setZero();
	_wtemp.setZero();
	_unp1.setZero();
	_vnp1.setZero();
	_wnp1.setZero();
	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	_wxstr.setZero();
	_wystr.setZero();
	_wzstr.setZero();
	_wxnp1.setZero();
	_wynp1.setZero();
	_wznp1.setZero();


	vort_confine_x.setZero();
	vort_confine_y.setZero();
	vort_confine_z.setZero();

	compute_curl();
	vorticity_confinement(dt*_vort_confine_coef);

	BFECC_field(dt,_un,_unp1,_utemp);
	BFECC_field(dt,_vn,_vnp1,_vtemp);
	BFECC_field(dt,_wn,_wnp1,_wtemp);
	clampExtrema(dt,_un,_unp1);
	clampExtrema(dt,_vn,_vnp1);
	clampExtrema(dt,_wn,_wnp1);



	stretch(dt); //_wxstr,_wystr,_wzstr = stretch(w,u);
	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	getDivergence();
	decay_vortices(dt, _wxstr, _wystr, _wzstr);


	BFECC_field(dt,_wxstr,_wxnp1,_wxn);
	BFECC_field(dt,_wystr,_wynp1,_wyn);
	BFECC_field(dt,_wzstr,_wznp1,_wzn);


	clampExtrema(dt, _wxstr,_wxnp1);
	clampExtrema(dt, _wystr,_wynp1);
	clampExtrema(dt, _wzstr,_wznp1);
	

	_Ttemp.setZero();
	_Ttempnp1.setZero();
	BFECC_field(dt,_rho,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_rho,_Ttempnp1);
	_rho.copy(_Ttempnp1);


	_Ttemp.setZero();
	_Ttempnp1.setZero();
	BFECC_field(dt,_Tbf,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_Tbf,_Ttempnp1);
	_Tbf.copy(_Ttempnp1);

	_Ttemp.setZero();
	_Ttempnp1.setZero();
	BFECC_field(dt,_fuel,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_fuel,_Ttempnp1);
	_fuel.copy(_Ttempnp1);

	//clampSmoke();

	_un.copy(_unp1);
	_vn.copy(_vnp1);
	_wn.copy(_wnp1);
	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	compute_curl();


	_Psix.setZero();
	_Psiy.setZero();
	_Psiz.setZero();
	formRHS(_hx*_hx*_hx*0.95,dt);
	_Psix.copy(_wxstr);
	_Psiy.copy(_wystr);
	_Psiz.copy(_wzstr);
	_wxstr.setZero();
	_wystr.setZero();
	_wzstr.setZero();
	formRHS(-_hx*_hx*0.95,dt);
	float res;

	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psix,&_wxstr,_hx);
	//_ppe_solver._b_Dirichlet = true;
	//_Psix.setZero();
	//_ppe_solver.m_Vcycle(&_Psix,&_wxstr, 1e-5,res,0);
	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psiy,&_wystr,_hx);
	//_ppe_solver._b_Dirichlet = true;
	//_Psiy.setZero();
	//_ppe_solver.m_Vcycle(&_Psiy,&_wystr, 1e-5,res,0);
	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psiz,&_wzstr,_hx);
	//_ppe_solver._b_Dirichlet = true;
	//_Psiz.setZero();
	//_ppe_solver.m_Vcycle(&_Psiz,&_wzstr, 1e-5,res,0);






	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psix,&_wxstr,_hx);
	solve_stream_Poisson(_Psix,_wxstr);
	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psiy,&_wystr,_hx);
	solve_stream_Poisson(_Psiy,_wystr);
	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psiz,&_wzstr,_hx);
	solve_stream_Poisson(_Psiz,_wzstr);





	addCurlPsi();
	add_vort_confinement(_un, vort_confine_x);
	add_vort_confinement(_vn, vort_confine_y);
	add_vort_confinement(_wn, vort_confine_z);
}

void SmokeSolver3D::clampSmoke()
{
	int compute_elements = _rho._blockx*_rho._blocky*_rho._blockz;

	int slice = _rho._blockx*_rho._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_rho._blockx;
		uint bi = thread_idx%(_rho._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<_nx && j<_ny && k<_nz)
			{
				//_rho(i,j,k) = max(0.0f,_rho(i,j,k));
				_Tbf(i,j,k) = min(max(_Tbf(i,j,k),0.0f),2000.0f);
				_fuel(i,j,k) = min(max(0.0f,_fuel(i,j,k)), 100.0f);
			}
		}
	});
}

void SmokeSolver3D::clampExtrema( float dt, buffer3Df & f_n, buffer3Df & f_np1 )
{
	int compute_elements = f_np1._blockx*f_np1._blocky*f_np1._blockz;

	int slice = f_np1._blockx*f_np1._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/f_np1._blockx;
		uint bi = thread_idx%(f_np1._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<f_np1._nx && j<f_np1._ny && k<f_np1._nz)
			{

				float world_x = ((float)i-f_np1._ox)*_hx;
				float world_y = ((float)j-f_np1._oy)*_hy;
				float world_z = ((float)k-f_np1._oz)*_hz;
				//Vec3f pos(world_x,world_y,world_z);
				//Vec3f trace_pos = trace(dt, pos);
				float u = _un.sample_linear(world_x,world_y,world_z);
				float v = _vn.sample_linear(world_x,world_y,world_z);
				float w = _wn.sample_linear(world_x,world_y,world_z);

				float px = world_x - 0.5*dt * u, py = world_y - 0.5*dt *v, pz = world_z - 0.5*dt*w;
				u = _un.sample_linear(px,py,pz);
				v = _vn.sample_linear(px,py,pz);
				w = _wn.sample_linear(px,py,pz);

				px = world_x - dt * u, py = world_y - dt *v, pz = world_z - dt*w;

				float v0,v1,v2,v3,v4,v5,v6,v7;
				//f_n.sample_cube(px,py,pz,v0,v1,v2,v3,v4,v5,v6,v7);
				float SLv = f_n.sample_cube_lerp(px,py,pz,
					v0,v1,v2,v3,v4,v5,v6,v7);

				float min_value = min(v0,min(v1,min(v2,min(v3,min(v4,min(v5,min(v6,v7)))))));
				float max_value = max(v0,max(v1,max(v2,max(v3,max(v4,max(v5,max(v6,v7)))))));

				if(f_np1(i,j,k)<min_value || f_np1(i,j,k)>max_value) 
				{

					f_np1(i,j,k) = SLv;
				}
				//f_np1(i,j,k) = max(min(max_value, f_np1(i,j,k)),min_value);

			}
		}
	});
}


void SmokeSolver3D::clampExtrema_order1( float dt, buffer3Df & f_n, buffer3Df & f_np1 )
{
	int compute_elements = f_np1._blockx*f_np1._blocky*f_np1._blockz;

	int slice = f_np1._blockx*f_np1._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/f_np1._blockx;
		uint bi = thread_idx%(f_np1._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<f_np1._nx && j<f_np1._ny && k<f_np1._nz)
			{

				float world_x = ((float)i-f_np1._ox)*_hx;
				float world_y = ((float)j-f_np1._oy)*_hy;
				float world_z = ((float)k-f_np1._oz)*_hz;
				//Vec3f pos(world_x,world_y,world_z);
				//Vec3f trace_pos = trace(dt, pos);
				float u = _un.sample_linear(world_x,world_y,world_z);
				float v = _vn.sample_linear(world_x,world_y,world_z);
				float w = _wn.sample_linear(world_x,world_y,world_z);



				float px = world_x - dt * u, py = world_y - dt *v, pz = world_z - dt*w;

				float v0,v1,v2,v3,v4,v5,v6,v7;
				//f_n.sample_cube(px,py,pz,v0,v1,v2,v3,v4,v5,v6,v7);
				float SLv = f_n.sample_cube_lerp(px,py,pz,
					v0,v1,v2,v3,v4,v5,v6,v7);

				float min_value = min(v0,min(v1,min(v2,min(v3,min(v4,min(v5,min(v6,v7)))))));
				float max_value = max(v0,max(v1,max(v2,max(v3,max(v4,max(v5,max(v6,v7)))))));


				if(f_np1(i,j,k)<min_value || f_np1(i,j,k)>max_value) 
				{

					f_np1(i,j,k) = SLv;
				}
				//f_np1(i,j,k) = max(min(max_value, f_np1(i,j,k)),min_value);

			}
		}
	});
}

void SmokeSolver3D::FLIP_advect( float dt )
{
	FLIP_Solver.FLIP(dt,_un,_vn,_wn,
		_unp1,_vnp1,_wnp1,
		_utemp,_vtemp,_wtemp,
		_rho,_Tbf);

	//_Ttemp.setZero();
	//_Ttempnp1.setZero();
	//BFECC_field(dt,_rho,_Ttempnp1,_Ttemp);
	//_rho.copy(_Ttempnp1);


	//_Ttemp.setZero();
	//_Ttempnp1.setZero();
	//BFECC_field(dt,_Tbf,_Ttempnp1,_Ttemp);
	//_Tbf.copy(_Ttempnp1);

	clampSmoke();
	_unp1.copy(_utemp);
	_vnp1.copy(_vtemp);
	_wnp1.copy(_wtemp);
	//instead of
	//_unp1.copy(_un);
	//_vnp1.copy(_vn);
	//_wnp1.copy(_wn);

	_un.copy(_utemp);
	_vn.copy(_vtemp);
	_wn.copy(_wtemp);



}
void SmokeSolver3D::FLIP_IVOCK(float dt)
{
	_utemp.setZero();
	_vtemp.setZero();
	_wtemp.setZero();


	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	_wxstr.setZero();
	_wystr.setZero();
	_wzstr.setZero();
	_wxnp1.setZero();
	_wynp1.setZero();
	_wznp1.setZero();


	vort_confine_x.setZero();
	vort_confine_y.setZero();
	vort_confine_z.setZero();

	compute_curl();
	vorticity_confinement(dt*_vort_confine_coef);

	FLIP_Solver.FLIP(dt,_un,_vn,_wn,
		_unp1,_vnp1,_wnp1,
		_utemp,_vtemp,_wtemp,
		_rho,_Tbf);

	printf("FLIP advection done\n");





	stretch(dt); //_wxstr,_wystr,_wzstr = stretch(w,u);
	_wxnp1.setZero();
	_wynp1.setZero();
	_wznp1.setZero();
	getDivergence();
	decay_vortices(dt, _wxstr, _wystr, _wzstr);
	advect_field_cubic_clamp(dt,_wxstr,_wxnp1);
	advect_field_cubic_clamp(dt,_wystr,_wynp1);
	advect_field_cubic_clamp(dt,_wzstr,_wznp1);
	

	_unp1.copy(_utemp);
	_vnp1.copy(_vtemp);
	_wnp1.copy(_wtemp);


	_un.copy(_utemp);
	_vn.copy(_vtemp);
	_wn.copy(_wtemp);


	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	compute_curl();
	_Psix.setZero();
	_Psiy.setZero();
	_Psiz.setZero();
	formRHS(0.95*_hx*_hx*_hx,dt);
	_Psix.copy(_wxstr);
	_Psiy.copy(_wystr);
	_Psiz.copy(_wzstr);
	_wxstr.setZero();
	_wystr.setZero();
	_wzstr.setZero();
	formRHS(-0.95*_hx*_hx,dt);
//	float res;
//	_ppe_solver.m_applyOpenBoundaryCondition(&_Psix,&_wxstr,_hx);
//	_ppe_solver._b_Dirichlet = true;
//	_Psix.setZero();
//	_ppe_solver.m_Vcycle(&_Psix,&_wxstr, 1e-5,res,0);
//	_ppe_solver.m_applyOpenBoundaryCondition(&_Psiy,&_wystr,_hx);
//	_ppe_solver._b_Dirichlet = true;
//	_Psiy.setZero();
//	_ppe_solver.m_Vcycle(&_Psiy,&_wystr, 1e-5,res,0);
//	_ppe_solver.m_applyOpenBoundaryCondition(&_Psiz,&_wzstr,_hx);
//	_ppe_solver._b_Dirichlet = true;
//	_Psiz.setZero();
//	_ppe_solver.m_Vcycle(&_Psiz,&_wzstr, 1e-5,res,0);


	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psix,&_wxstr,_hx);
	solve_stream_Poisson(_Psix,_wxstr);
	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psiy,&_wystr,_hx);
	solve_stream_Poisson(_Psiy,_wystr);
	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psiz,&_wzstr,_hx);
	solve_stream_Poisson(_Psiz,_wzstr);




	addCurlPsi();
	add_vort_confinement(_un, vort_confine_x);
	add_vort_confinement(_vn, vort_confine_y);
	add_vort_confinement(_wn, vort_confine_z);
	printf("IVOCK done\n");
}

void SmokeSolver3D::MacCormack( float dt )
{
	_utemp.setZero();
	_vtemp.setZero();
	_wtemp.setZero();
	_unp1.setZero();
	_vnp1.setZero();
	_wnp1.setZero();
	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	vort_confine_x.setZero();
	vort_confine_y.setZero();
	vort_confine_z.setZero();
	compute_curl();
	vorticity_confinement(dt*_vort_confine_coef);



	MacCormack_field(dt,_un,_unp1,_utemp);
	MacCormack_field(dt,_vn,_vnp1,_vtemp);
	MacCormack_field(dt,_wn,_wnp1,_wtemp);
	clampExtrema(dt,_un,_unp1);
	clampExtrema(dt,_vn,_vnp1);
	clampExtrema(dt,_wn,_wnp1);

	_Ttemp.setZero();
	_Ttempnp1.setZero();
	MacCormack_field(dt,_rho,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_rho,_Ttempnp1);
	_rho.copy(_Ttempnp1);


	_Ttemp.setZero();
	_Ttempnp1.setZero();
	MacCormack_field(dt,_Tbf,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_Tbf,_Ttempnp1);
	_Tbf.copy(_Ttempnp1);

	_Ttemp.setZero();
	_Ttempnp1.setZero();
	MacCormack_field(dt,_fuel,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_fuel,_Ttempnp1);
	_fuel.copy(_Ttempnp1);

	clampSmoke();


	_un.copy(_unp1);
	_vn.copy(_vnp1);
	_wn.copy(_wnp1);


	add_vort_confinement(_un, vort_confine_x);
	add_vort_confinement(_vn, vort_confine_y);
	add_vort_confinement(_wn, vort_confine_z);
}

void SmokeSolver3D::MacCormack_IVOCK( float dt )
{
	mark_boundary();

	_utemp.setZero();
	_vtemp.setZero();
	_wtemp.setZero();
	_unp1.setZero();
	_vnp1.setZero();
	_wnp1.setZero();
	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	_wxstr.setZero();
	_wystr.setZero();
	_wzstr.setZero();
	_wxnp1.setZero();
	_wynp1.setZero();
	_wznp1.setZero();

	vort_confine_x.setZero();
	vort_confine_y.setZero();
	vort_confine_z.setZero();

	compute_curl();
	vorticity_confinement(dt*_vort_confine_coef);


	MacCormack_field(dt,_un,_unp1,_utemp);
	MacCormack_field(dt,_vn,_vnp1,_vtemp);
	MacCormack_field(dt,_wn,_wnp1,_wtemp);
	clampExtrema(dt,_un,_unp1);
	clampExtrema(dt,_vn,_vnp1);
	clampExtrema(dt,_wn,_wnp1);



	stretch(dt); //_wxstr,_wystr,_wzstr = stretch(w,u);
	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	getDivergence();
	decay_vortices(dt, _wxstr, _wystr, _wzstr);


	MacCormack_field(dt,_wxstr,_wxnp1,_wxn);
	MacCormack_field(dt,_wystr,_wynp1,_wyn);
	MacCormack_field(dt,_wzstr,_wznp1,_wzn);
	clampExtrema(dt, _wxstr,_wxnp1);
	clampExtrema(dt, _wystr,_wynp1);
	clampExtrema(dt, _wzstr,_wznp1);
	

	_Ttemp.setZero();
	_Ttempnp1.setZero();
	MacCormack_field(dt,_rho,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_rho,_Ttempnp1);
	_rho.copy(_Ttempnp1);


	_Ttemp.setZero();
	_Ttempnp1.setZero();
	MacCormack_field(dt,_Tbf,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_Tbf,_Ttempnp1);
	_Tbf.copy(_Ttempnp1);

	_Ttemp.setZero();
	_Ttempnp1.setZero();
	MacCormack_field(dt,_fuel,_Ttempnp1,_Ttemp);
	clampExtrema(dt,_fuel,_Ttempnp1);
	_fuel.copy(_Ttempnp1);

	clampSmoke();

	_un.copy(_unp1);
	_vn.copy(_vnp1);
	_wn.copy(_wnp1);
	_wxn.setZero();
	_wyn.setZero();
	_wzn.setZero();
	compute_curl();


	_Psix.setZero();
	_Psiy.setZero();
	_Psiz.setZero();
	formRHS(_hx*_hx*_hx*0.95,dt);
	_Psix.copy(_wxstr);
	_Psiy.copy(_wystr);
	_Psiz.copy(_wzstr);
	_wxstr.setZero();
	_wystr.setZero();
	_wzstr.setZero();
	formRHS(-_hx*_hx*0.95,dt);
//	float res;
//	_ppe_solver.m_applyOpenBoundaryCondition(&_Psix,&_wxstr,_hx);
//	_ppe_solver._b_Dirichlet = true;
//	_Psix.setZero();
//	_ppe_solver.m_Vcycle(&_Psix,&_wxstr, 1e-5,res,0);
//	_ppe_solver.m_applyOpenBoundaryCondition(&_Psiy,&_wystr,_hx);
//	_ppe_solver._b_Dirichlet = true;
//	_Psiy.setZero();
//	_ppe_solver.m_Vcycle(&_Psiy,&_wystr, 1e-5,res,0);
//	_ppe_solver.m_applyOpenBoundaryCondition(&_Psiz,&_wzstr,_hx);
//	_ppe_solver._b_Dirichlet = true;
//	_Psiz.setZero();
//	_ppe_solver.m_Vcycle(&_Psiz,&_wzstr, 1e-5,res,0);


	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psix,&_wxstr,_hx);
	solve_stream_Poisson(_Psix,_wxstr);
	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psiy,&_wystr,_hx);
	solve_stream_Poisson(_Psiy,_wystr);
	//_ppe_solver.m_applyOpenBoundaryCondition(&_Psiz,&_wzstr,_hx);
	solve_stream_Poisson(_Psiz,_wzstr);





	addCurlPsi();
	add_vort_confinement(_un, vort_confine_x);
	add_vort_confinement(_vn, vort_confine_y);
	add_vort_confinement(_wn, vort_confine_z);
}

void SmokeSolver3D::MacCormack_field( float dt, buffer3Df &field_n,buffer3Df &field_np1,buffer3Df &field_AUX )
{
	simple_advect_field(dt,field_n,field_AUX);
	simple_advect_field(-dt,field_AUX,field_np1);
	int compute_elements = field_AUX._blockx*field_AUX._blocky*field_AUX._blockz;
	int slice = field_AUX._blockx*field_AUX._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx)
	{
		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/field_AUX._blockx;
		uint bi = thread_idx%(field_AUX._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<field_AUX._nx && j<field_AUX._ny && k<field_AUX._nz)
			{
				field_np1(i,j,k) = field_AUX(i,j,k) + 0.5*(field_n(i,j,k)-field_np1(i,j,k));
			}
		}
	});
}

void SmokeSolver3D::MacCormack_field_order1( float dt, buffer3Df &field_n,buffer3Df &field_np1,buffer3Df &field_AUX )
{
	simple_advect_field_order1(dt,field_n,field_AUX);
	simple_advect_field_order1(-dt,field_AUX,field_np1);
	int compute_elements = field_AUX._blockx*field_AUX._blocky*field_AUX._blockz;
	int slice = field_AUX._blockx*field_AUX._blocky;

	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx)
	{
		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/field_AUX._blockx;
		uint bi = thread_idx%(field_AUX._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			if(i<field_AUX._nx && j<field_AUX._ny && k<field_AUX._nz)
			{
				field_np1(i,j,k) = field_AUX(i,j,k) + 0.5*(field_n(i,j,k)-field_np1(i,j,k));
			}
		}
	});
}

void SmokeSolver3D::clearBoundary()
{
	int compute_elements = _rho._blockx*_rho._blocky*_rho._blockz;

	int slice = _rho._blockx*_rho._blocky;
	u_valid.set_zero();
	v_valid.set_zero();
	w_valid.set_zero();
	tbb::parallel_for(0, compute_elements, 1, [&](int thread_idx) {

		uint bk = thread_idx/slice;
		uint bj = (thread_idx%slice)/_rho._blockx;
		uint bi = thread_idx%(_rho._blockx);

		for (uint kk=0;kk<8;kk++)for(uint jj=0;jj<8;jj++)for(uint ii=0;ii<8;ii++)
		{
			uint i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;
			float y = (float)(j+0.5)/((float)_ny);
			//if(y>0.8 && y<=0.9)
			//{
			//	float coef = 1.8-y;

			//	_Tbf(i,j,k) *= 0;
			//	_rho(i,j,k) *= 0;
			//	_fuel(i,j,k) *= 0;
			//	_un(i,j,k) *= coef;
			//	_un(i+1,j,k) *= coef;

			//	_vn(i,j,k) *= coef;
			//	_vn(i,j+1,k) *= coef;

			//	_wn(i,j,k) *= coef ;
			//	_wn(i,j,k+1) *= coef;
			//}
			if(y>0.9 && y<=1)
			{
				float coef = (1-y)/0.1;

				_Tbf(i,j,k) *= 0;
				_rho(i,j,k) *= 0;
				_fuel(i,j,k) *= 0;
				//_un(i,j,k) *= coef;
				//_un(i+1,j,k) *= coef;

				//_vn(i,j,k) *= coef;
				//_vn(i,j+1,k) *= coef;

				//_wn(i,j,k) *= coef ;
				//_wn(i,j,k+1) *= coef;
			}

			if(_b_desc(i,j,k)==1)
			{
				_un(i,j,k) *= 0;
				_un(i+1,j,k) *= 0;

				_vn(i,j,k) *= 0;
				_vn(i,j+1,k) *= 0;

				_wn(i,j,k) *= 0 ;
				_wn(i,j,k+1) *= 0;
			}
		}
	});
}

void SmokeSolver3D::pcg_projection()
{

	int ni = _nx;
	int nj = _ny;
	int nk = _nz;

	int system_size = ni*nj*nk;
	if(rhs.size() != system_size) {
		rhs.resize(system_size);
		pressure.resize(system_size);
		matrix.resize(system_size);
	}

	matrix.zero();
	rhs.assign(rhs.size(), 0);
	pressure.assign(pressure.size(), 0);
	//write boundary velocity;
	int compute_num = ni*nj*nk;
	int slice = ni*nj;
	tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/ni;
		int i = thread_idx%ni;
		if ( _b_desc(i,j,k)==2 )//solid
		{
			_un(i,j,k) = 0;
			_un(i+1,j,k) = 0;
			_vn(i,j,k) = 0;
			_vn(i,j+1,k) = 0;
			_wn(i,j,k) = 0;
			_wn(i,j,k+1) = 0;
		}
	});


	//set up solver
	compute_num = ni*nj*nk;
	slice = ni*nj;
	tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/ni;
		int i = thread_idx%ni;
		if(i>=1 && i<ni-1 && j>=1 && j<nj-1 && k>=1 && k<nk-1)
		{
			int index = i + ni*j + ni*nj*k;

			rhs[index] = 0;
			pressure[index] = 0;

			if( _b_desc(i,j,k)==0 )//a fluid cell 
			{

				//right neighbour
				if( _b_desc(i+1,j,k)==0 ) {//a fluid cell
					matrix.add_to_element(index, index, 1.0/_hx/_hx);
					matrix.add_to_element(index, index + 1, -1.0/_hx/_hx);
				}
				else if( _b_desc(i+1,j,k)==1 )//an empty cell
				{
					matrix.add_to_element(index, index, 1.0/_hx/_hx);
				}
				rhs[index] -= _un(i+1,j,k) / _hx;

				//left neighbour
				if( _b_desc(i-1,j,k)==0 ) {
					matrix.add_to_element(index, index, 1.0/_hx/_hx);
					matrix.add_to_element(index, index - 1, -1.0/_hx/_hx);
				}
				else if( _b_desc(i-1,j,k)==1 ){

					matrix.add_to_element(index, index, 1.0/_hx/_hx);
				}
				rhs[index] += _un(i,j,k) / _hx;



				//top neighbour
				if( _b_desc(i,j+1,k)==0 ) {//a fluid cell
					matrix.add_to_element(index, index, 1.0/_hx/_hx);
					matrix.add_to_element(index, index + ni, -1.0/_hx/_hx);
				}
				else if( _b_desc(i,j+1,k)==1 )//an empty cell
				{
					matrix.add_to_element(index, index, 1.0/_hx/_hx);
				}
				rhs[index] -= _vn(i,j+1,k) / _hx;

				//bottom neighbour
				if( _b_desc(i,j-1,k)==0 ) {
					matrix.add_to_element(index, index, 1.0/_hx/_hx);
					matrix.add_to_element(index, index - ni, -1.0/_hx/_hx);
				}
				else if( _b_desc(i,j-1,k)==1 ){

					matrix.add_to_element(index, index, 1.0/_hx/_hx);
				}
				rhs[index] += _vn(i,j,k) / _hx;
				//rhs[index] += _burn_div(i,j,k);



				//back neighbour
				if( _b_desc(i,j,k+1)==0 ) {//a fluid cell
					matrix.add_to_element(index, index, 1.0/_hx/_hx);
					matrix.add_to_element(index, index + ni*nj, -1.0/_hx/_hx);
				}
				else if( _b_desc(i,j,k+1)==1 )//an empty cell
				{
					matrix.add_to_element(index, index, 1.0/_hx/_hx);
				}
				rhs[index] -= _wn(i,j,k+1) / _hx;

				//front neighbour
				if( _b_desc(i,j,k-1)==0 ) {
					matrix.add_to_element(index, index, 1.0/_hx/_hx);
					matrix.add_to_element(index, index - ni*nj, -1.0/_hx/_hx);
				}
				else if( _b_desc(i,j,k-1)==1 ){

					matrix.add_to_element(index, index, 1.0/_hx/_hx);
				}
				rhs[index] += _wn(i,j,k) / _hx;


				rhs[index] += _burn_div(i,j,k);
			}
		}
	});

	//Solve the system using a AMGPCG solver

	double tolerance;
	int iterations;
	//solver.set_solver_parameters(1e-6, 1000);
	//bool success = solver.solve(matrix, rhs, pressure, tolerance, iterations);
	bool success = AMGPCGSolve(matrix,rhs,pressure,1e-10,1000,tolerance,iterations,_nx,_ny,_nz);
	printf("Solver took %d iterations and had residual %e\n", iterations, tolerance);
	if(!success) {
		printf("WARNING: Pressure solve failed!************************************************\n");
	}

	//apply grad
	u_valid.assign(0);
	compute_num = _un._nx*_un._ny*_un._nz;
	slice = _un._nx*_un._ny;
	tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/_un._nx;
		int i = thread_idx%_un._nx;
		if(k<_un._nz && j<_un._ny && i<_un._nx-1 && i>0)
		{
			int index = i + j*ni + k*ni*nj;
			if(_b_desc(i,j,k) == 0 || _b_desc(i-1,j,k) == 0) {

				_un(i,j,k) -=  (float)(pressure[index] - pressure[index-1]) / _hx ; 
				u_valid(i,j,k) = 1;
			}

		}
	});

	v_valid.assign(0);
	compute_num = _vn._nx*_vn._ny*_vn._nz;
	slice = _vn._nx*_vn._ny;
	tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/_vn._nx;
		int i = thread_idx%_vn._nx;
		if(k<_vn._nz && j>0 && j<_vn._ny-1 && i<_vn._nx )
		{
			int index = i + j*ni + k*ni*nj;
			if(_b_desc(i,j,k) == 0 || _b_desc(i,j-1,k) == 0) {

				_vn(i,j,k) -=  (float)(pressure[index] - pressure[index-ni]) / _hx ; 
				v_valid(i,j,k) = 1;
			}

		}
	});

	w_valid.assign(0);
	compute_num = _wn._nx*_wn._ny*_wn._nz;
	slice = _wn._nx*_wn._ny;
	tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/_wn._nx;
		int i = thread_idx%_wn._nx;
		if(k>0 && k<_wn._nz-1 && j<_wn._ny && i<_wn._nx )
		{
			int index = i + j*ni + k*ni*nj;
			if(_b_desc(i,j,k) == 0 || _b_desc(i,j,k-1) == 0) {

				_wn(i,j,k) -=  (float)(pressure[index] - pressure[index-ni*nj]) / _hx ; 
				w_valid(i,j,k) = 1;
			}

		}
	});
	//write boundary velocity
	compute_num = ni*nj*nk;
	slice = ni*nj;
	tbb::parallel_for(0,compute_num,1,[&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/ni;
		int i = thread_idx%ni;
		if ( _b_desc(i,j,k)==2 )//solid
		{
			_un(i,j,k) = 0;
			u_valid(i,j,k) = 1;
			_un(i+1,j,k) = 0;
			u_valid(i+1,j,k) =1;
			_vn(i,j,k) = 0;
			v_valid(i,j,k) = 1;
			_vn(i,j+1,k) = 0;
			v_valid(i,j+1,k) = 1;
			_wn(i,j,k) = 0;
			w_valid(i,j,k) = 1;
			_wn(i,j,k+1) = 0;
			w_valid(i,j,k+1) =1;
		}
	});

	compute_num = _un._nx*_un._ny*_un._nz;
	slice = _un._nx*_un._ny;
	tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/_un._nx;
		int i = thread_idx%_un._nx;
		if(k<_un._nz&& j<_un._ny && i<_un._nx )
		{
			if(u_valid(i,j,k)==0)
			{
				_un(i,j,k) = 0;
			}
		}
	});

	compute_num = _vn._nx*_vn._ny*_vn._nz;
	slice = _vn._nx*_vn._ny;
	tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/_vn._nx;
		int i = thread_idx%_vn._nx;
		if(k<_vn._nz&& j<_vn._ny && i<_vn._nx )
		{
			if(v_valid(i,j,k)==0)
			{
				_vn(i,j,k) = 0;
			}
		}
	});

	compute_num = _wn._nx*_wn._ny*_wn._nz;
	slice = _wn._nx*_wn._ny;
	tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/_wn._nx;
		int i = thread_idx%_wn._nx;
		if(k<_wn._nz&& j<_wn._ny && i<_wn._nx )
		{
			if(w_valid(i,j,k)==0)
			{
				_wn(i,j,k) = 0;
			}
		}
	});
	//extrapolate(u_extrap,_un,u_valid);
	//extrapolate(v_extrap,_vn,v_valid);
	//extrapolate(w_extrap,_wn,w_valid);
}

void SmokeSolver3D::extrapolate(Array3f & grid, buffer3Df & u, Array3c & valid)
{
	grid.assign(0);
	int compute_num = u._nx*u._ny*u._nz;
	int slice = u._nx*u._ny;
	tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/u._nx;
		int i = thread_idx%u._nx;
		if(k<u._nz && j<u._ny && i<u._nx )
		{
			grid(i,j,k) = u(i,j,k);
		}
	});



	Array3f temp_grid = grid;
	Array3c old_valid(valid.ni,valid.nj,valid.nk);
	for(int layers = 0; layers < 10; ++layers) {
		old_valid = valid;
		int num = grid.ni*grid.nj*grid.nk;
		int s = grid.ni*grid.nj;
		tbb::parallel_for(0,num,1,[&](int thread_idx)
			//for(int k = 1; k < grid.nk-1; ++k) for(int j = 1; j < grid.nj-1; ++j) for(int i = 1; i < grid.ni-1; ++i) 
		{
			int k = thread_idx/s;
			int j = (thread_idx%s)/grid.ni;
			int i = thread_idx%grid.ni;
			if(k>=1 && k<grid.nk-1 && j>=1 && j<grid.nj-1 && i>=1 && i<grid.ni-1){

				float sum = 0;
				int count = 0;

				if(!old_valid(i,j,k)) 
				{

					if(old_valid(i+1,j,k)) {
						sum += grid(i+1,j,k);
						++count;
					}
					if(old_valid(i-1,j,k)) {
						sum += grid(i-1,j,k);
						++count;
					}
					if(old_valid(i,j+1,k)) {
						sum += grid(i,j+1,k);
						++count;
					}
					if(old_valid(i,j-1,k)) {
						sum += grid(i,j-1,k);
						++count;
					}
					if(old_valid(i,j,k+1)) {
						sum += grid(i,j,k+1);
						++count;
					}
					if(old_valid(i,j,k-1)) {
						sum += grid(i,j,k-1);
						++count;
					}

					//If any of neighbour cells were valid, 
					//assign the cell their average value and tag it as valid
					if(count > 0) {
						temp_grid(i,j,k) = sum /(float)count;
						valid(i,j,k) = 1;
					}

				}
			}
		});
		grid = temp_grid;

	}



	compute_num = u._nx*u._ny*u._nz;
	slice = u._nx*u._ny;
	tbb::parallel_for(0, compute_num, 1, [&](int thread_idx)
	{
		int k = thread_idx/slice;
		int j = (thread_idx%slice)/u._nx;
		int i = thread_idx%u._nx;
		if(k<u._nz && j<u._ny && i<u._nx )
		{
			u(i,j,k) = grid(i,j,k);
		}
	});


}