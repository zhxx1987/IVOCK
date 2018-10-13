#include <cmath>
#include "tbb/tbb.h"
//#include "Multigrid3D.h"
#include "array.h"
#include "fluid_buffer3D.h"
#include "Smoke_solver3D.h"
#include <iostream>

#include <ctime>
//using namespace gf;
double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks)/(CLOCKS_PER_SEC/1000);
	return diffms/1000.0;
}

int main(int argc, char** argv) {

	//int dimx=128, dimy=256, dimz=128;
	//int iter=1;
	//if(argc>3)
	//{
	//	sscanf(argv[1], "%d", &dimx);
	//	sscanf(argv[2], "%d", &dimy);
	//	sscanf(argv[3], "%d", &dimz);
	//	sscanf(argv[4], "%d", &iter);
	//}
	//buffer3Dd x1,b1, res1;
	//MGSolverd mg_solver;
	//mg_solver.m_InitialSystem(dimx, dimy, dimz);


	//x1.init(dimx, dimy, dimz);
	//b1.init(dimx, dimy, dimz);
	//res1.init(dimx, dimy, dimz);
	//x1.setZero();
	//b1.setZero();
	//res1.setZero();

	//for (int k=1; k<dimz-1; k++) for(int i=1;i<dimx-1;i++)
	//{
	//	b1(i,0,k) = -1.0f;
	//}
	//for (int k=1; k<dimz-1; k++) for(int i=1;i<dimx-1;i++)
	//{
	//	b1(i,dimy-1,k) = 1.0f;
	//}

	//double resd;
	//printf("problem size: %dx%dx%d\n", dimx, dimy, dimz);

	//clock_t start = clock();
	//mg_solver._b_Dirichlet = true;
	//for(int i=0; i<3; i++){
	//	mg_solver.m_applyOpenBoundaryCondition(&b1,&b1,1.0/128);
	//	mg_solver.m_Vcycle(&x1, &b1, 1e-5, resd,0);
	//}
	//

	//printf("full multi grid done!\n");
	//clock_t end = clock();
	//cout << diffclock(end,start)<<endl;


	//int modx = dimx-1;
	//int mody = (dimy-1)*dimx;
	//int modz = (dimz-1)*dimx*dimy;
	//for (int k=0; k<dimz; k++) for(int j=0; j<dimy; j++) for(int i=0;i<dimx;i++)
	//{
	//	
	//	int left_idx = (i==0)? i  : i -1;
	//	int right_idx = (i==(dimx-1))?i : i + 1;
	//	int up_idx = (j==(dimy-1))?j : j + 1;
	//	int down_idx = (j==0)?j : j - 1;
	//	int front_idx = (k==0)?k : k - 1;
	//	int back_idx = (k==(dimz-1))?k : k + 1;
	//	double lv = x1(left_idx,j,k);
	//	if(i==0) lv = 0;
	//	double rv = x1(right_idx,j,k);
	//	if(i==dimx-1) rv = 0;
	//	double tv = x1(i,up_idx,k);
	//	if(j==dimy-1) tv =0;
	//	double dv = x1(i,down_idx,k);
	//	if(j==0) dv = 0;
	//	double fv = x1(i,j,front_idx);
	//	if(k==0) fv = 0;
	//	double bv = x1(i,j,back_idx);
	//	if(k==dimz-1) bv = 0;

	//	res1(i,j,k) = b1(i,j,k)
	//	- (-6*x1(i,j,k)
	//	+ lv + rv + tv + dv + fv + bv);
	//}
	//double sum=0;
	//for (int k=0; k<dimz-0; k++) for(int j=0; j<dimy-0; j++) for(int i=0;i<dimx-0;i++)
	//{

	//	//if(b1.getHostPtr()[k*dimx*dimy+j*dimx+i]>0.0000001)
	//		if(fabs(res1(i,j,k))>sum)
	//			sum = fabs(res1(i,j,k));
	//	//printf("%f, ", res.getHostPtr()[k*dimx*dimy+j*dimx+i]);
	//}
	//printf("%e\n",sum);

	//mg_solver.m_FinalMemoryEachLevel();
	//x1.free();
	//b1.free();
	//res1.free();
	///////////////////////////////////////////////////////////////////////////////////
	if(argc<4)
	{
		printf("please specify output path, advection type(0:SL, 1:SL+iVOCK, 2:BFECC 3:BFECC+iVOCK, 4:MacCormack, 5:MacCormack+iVOCK, 6:FLIP, 7:FLIP+iVOCK, 8:Best Combination),vort_confine_str \n");
		return 0;
	}
	else{
		int advection_type;
		float vort_confine_str=0;

		char file_path[256];
		int n=sprintf(file_path,"%s",argv[1]);
		sscanf(argv[2],"%d",&advection_type);
		sscanf(argv[3], "%f", &vort_confine_str);
		if(advection_type<0 || advection_type>8)
		{
			printf("error!, advection type must take one of these values : 0:SL, 1:SL+iVOCK, 2:BFECC 3:BFECC+iVOCK, 4:MacCormack, 5:MacCormack+iVOCK, 6:FLIP, 7:FLIP+iVOCK, 8:Best Combination\n");
			return 0;
		}

		switch(advection_type)
		{
		case 0:
			cout<<"use Semi-Lagrangian advection\n";
			break;
		case 1:
			cout<<"use SL_iVOCK advection\n";
			break;
		case 2:
			cout<<"use BFECC advection\n";
			break;
		case 3:
			cout<<"use BFECC-IVOCK advection\n";
			break;
		case 4:
			cout<<"use MacCormack advection\n";
			break;
		case 5:
			cout<<"use MacCormack-IVOCK advection\n";
			break;
		case 6:
			cout<<"use FLIP advection\n";
			break;
		case 7:
			cout<<"use FLIP-IVOCK advection\n";
			break;
		case 8:
			cout<<"use Best Combination\n";
			break;
		default:
			break;
		}

		int nx=128,ny=256,nz=128;
		float L = 20.0;
		float g_h = L/(float)nx;
		SmokeSolver3D g_smokeSolver;
		g_smokeSolver.init(nx,ny,nz,L);
		g_smokeSolver.setSmoke(4.0,0.5,0.35,500,10.0);


		buffer3Dc bc_des;
		bc_des.init(nx,ny,nz);
		bc_des.setZero();
		for (int k=0;k<nz;k++)for(int j=0;j<ny;j++)for(int i=0;i<nx;i++)
		{
			//0:fluid;1:air;2:solid
			if(i<1) bc_des(i,j,k) = 1;
			if(j<1) bc_des(i,j,k) = 2;
			if(k<1) bc_des(i,j,k) = 1;

			if(i>=nx-1) bc_des(i,j,k) = 1;
			if(j>=ny-1) bc_des(i,j,k) = 1;
			if(k>=nz-1) bc_des(i,j,k) = 1;

			float x = g_h * i;
			float y = g_h * j;
			float z = g_h * k;
			float X = g_h * (float)nx;
			float Y = g_h * (float)ny;
			float Z = g_h * (float)nz;
			//if (sqrt((x-0.5*X)*(x-0.5*X)+(y-0.5*Y)*(y-0.5*Y)+(z-0.5*Z)*(z-0.5*Z))<=0.1*L)
//			if (sqrt((x-0.52*L)*(x-0.52*L)+(y-0.125*L)*(y-0.125*L)+(z-0.5*L)*(z-0.5*L))<=0.02*L)
//			{
//				bc_des(i,j,k) = 2;
//			}

		}
		g_smokeSolver.set_boundary(bc_des);

		bc_des.setZero();
		for (int k=0;k<nz;k++)for(int j=0;j<ny;j++)for(int i=0;i<nx;i++)
		{
			float x = g_h * i;
			float y = g_h * j;
			float z = g_h * k;
			//1:source; 2:clear
			if (sqrt((x-0.5*L)*(x-0.5*L)+(y-0.12*L)*(y-0.12*L)+(z-0.5*L)*(z-0.5*L))<=0.07*L)
			{
				bc_des(i,j,k) = 1;
			}
		}
		g_smokeSolver.set_heat(bc_des);
		g_smokeSolver.setEmitter(Vec3f(0.5*L, 0.12*L, 0.5*L),0.07*L,12288);
		



		clock_t start = clock();
		for (int frame = 0; frame<200;frame++)
		{
			//for(int subs=0; subs<2;subs++)
			float T = 0.1*(float)frame;


			//g_smokeSolver.setSmoke(0.2,0.25,0.175,500/*+10*sin(T/20.0*32*3.1415926)*/,10.0);
			g_smokeSolver.setSmoke(10.0,0.15,0.25,500,10.0);
			g_smokeSolver.set_vort_confine(vort_confine_str);



			//for(int subs=0; subs<2;subs++)
			g_smokeSolver.time_step(0.02, advection_type);

//			g_smokeSolver.write_bgeo(g_smokeSolver._nx,
//				g_smokeSolver._ny,
//				g_smokeSolver._nz,
//				frame,
//				file_path);
			g_smokeSolver.write_tracers_bgeo(file_path,frame);
			printf("frame %d done\n",frame);
		}
		clock_t end = clock();
		cout << diffclock(end,start)/200.0<<endl;




		return 0;
	}
	













	//buffer obj test
	//	buffer3Df a;
	//a.init(4,4,4);
	//tbb::parallel_for(0,64,1,[&a](int tidx)
	//{
	//	uint k = tidx/16;
	//	uint j = (tidx%16)/4;
	//	uint i = (tidx%16)%4;

	//	a(i,j,k)= 1.0f;
	//});

	//tbb::parallel_for(0,64,1,[&a](int tidx)
	//{
	//	uint k = tidx/16;
	//	uint j = (tidx%16)/4;
	//	uint i = (tidx%16)%4;

	//	a(i,j,k) = a(i,j,k) + (float)i;
	//});

	//for (uint k=0;k<4;k++)for(uint j=0;j<4;j++)for(uint i=0;i<4;i++)
	//{
	//	a(i,j,k) += 1.0;
	//	printf("%f,",a(i,j,k));
	//}


}