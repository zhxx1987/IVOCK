#include <cstdio>
#include "read_particles.h"
#include "render.h"
#include "shadow.h"
#include "write_sgi.h"
#include "obj_reader.h"
double pi = 3.14159265359;


double frand(double a, double b)
{
	return (((double)rand())/((double)RAND_MAX)*(b-a))+a;
}
int
main(int argc,
	 char **argv)
{
	if(argc<8){
		std::printf("usage: krak <particle_file_format> <radius> <density_scale> <opacity> <light_x> <light_y> <light_z> [quality] [start_frame] [frame_skip]\n");
		return 0;
	}


	const char *filename_format=argv[1];


	Array1<Vec3f> x, base_rgb, lit_rgb;
	Array1<float> d,temp;
	float radius, opacity;
	Vec3f light_position;
	float density_scale;
	std::sscanf(argv[2], "%f", &radius);
	std::sscanf(argv[3], "%f", &density_scale);
	std::sscanf(argv[4], "%f", &opacity);
	std::sscanf(argv[5], "%f", &light_position[0]);
	std::sscanf(argv[6], "%f", &light_position[1]);
	std::sscanf(argv[7], "%f", &light_position[2]);

	float quality_factor=1;
	int start_frame=1, frame_skip=1;
	if(argc>8) std::sscanf(argv[8], "%f", &quality_factor);
	if(argc>9) std::sscanf(argv[9], "%d", &start_frame);
	if(argc>10) std::sscanf(argv[10], "%d", &frame_skip);

	std::printf("radius %f, opacity %f, light (%f %f %f)\n", radius, opacity,
		light_position[0], light_position[1], light_position[2]);

	Array1<float> illum;
	Array2f shadow_map;
	Array2<Vec3f> image(1440,960);
	Array1<Vec3f> ball;

	/*char objfile[128];
	int n=sprintf(objfile, "C:/Users/xinxin/Desktop/bunny4.obj");

	vector<Vec3f> obj_vertices;
	vector<Vec3ui> obj_faces;

	load_obj(objfile, obj_vertices,obj_faces);
	normalize_object(obj_vertices,obj_faces,1.8);
	ball.resize(obj_vertices.size());

	for (int i=0; i<ball.size();i++)
	{
		ball[i] = obj_vertices[i];
	}*/
	//ball.resize(0);

	ball.resize(163840);
	int num_get = 0;
	while(num_get<163840)
	{
		float x = frand(-0.2,0.2);
		float y = frand(-0.2,0.2);
		float z = frand(-0.2,0.2);
		if(sqrt(x*x+y*y+z*z)<=0.105&&sqrt(x*x+y*y+z*z)>=0.095)
		{
			ball[num_get] = Vec3f(x,y+1.0,z);
			num_get++;
		}
	}
	for(int f=start_frame; f<frame_skip;f++){

		//Vec3f ball_position = Vec3f(0 + 0.005 * 4.0 - 0.01*4.0*(double)f, 0,0);
		Array1<Vec3f> y;


		
		while(!read_particles(x, d, "%s/Particle_data%04d.bin", filename_format, f)){
			
		}
		
		char filename[256];
		int n = sprintf(filename, "%s/temp_data%04d.bin", filename_format,f);

		FILE *temp_data = fopen(filename, "rb");
		if(temp_data!=NULL)
		{
			temp.resize(x.size());
			size_t result = fread(&(temp[0]),1,sizeof(float)*temp.size(), temp_data);
			fclose(temp_data);
		}
		y.resize(x.size());
		base_rgb.resize(x.size());
		lit_rgb.resize(x.size());
		for(int p=0;p<x.size();p++)
		{
			y[p]=x[p];
			//base_rgb[p]=Vec3f(.2, .4, .6);
			//lit_rgb[p]=Vec3f(1.0, 1, 1);


			base_rgb[p]=Vec3f(0.5,0.7,0.9);
			lit_rgb[p]=Vec3f(.9,.9,.99);

			/*base_rgb[p]=Vec3f(0.1,0.08,0.07);
			lit_rgb[p]=Vec3f(0.4,0.4,0.5);*/
			if(temp.size()>0)
			{
				float heat_index = (temp[p])/1000.0f;
				heat_index = max(min(heat_index, 1.0f),0.0f);

				float r = heat_index;
				float g = heat_index*0.4;
				float b = 0.01;

				base_rgb[p]=Vec3f(r, g, b)*3+Vec3f(0.08,0.1,0.15);
				lit_rgb[p]=Vec3f(r, g, b)*2+Vec3f(0.6,0.6,0.6);
			}
		}
		//for (int p = 0; p<163840;p++)
		//{
		//	y.push_back(ball[p]);

		//	base_rgb.push_back(Vec3f(0.5,0.2,0.2));
		//	lit_rgb.push_back(Vec3f(.8, .8, .8));
		//	d.push_back(2.0);

		//}
		for (unsigned int ii=0;ii<d.size();ii++)
		{
			d[ii] *= density_scale;
		}
		
		if(!compute_shadows(y, radius, opacity, d, light_position, quality_factor, illum, shadow_map)){
			break;
		}
		write_sgi(shadow_map, false, "%s/shadow%04d.sgi", filename_format, f);

		//for(int i=0;i<60;i++)
		{
			/*double theta=1.65*pi, phi=7*0.25+2.0/60.0*(double)i*pi, camr=30;
			double cam_x = camr*sin(theta)*cos(phi);
			double cam_y = camr*sin(theta)*sin(phi);
			double cam_z = camr*cos(theta);*/
			render_smoke(y, illum, base_rgb, lit_rgb,radius, opacity, d, light_position,  Vec3f(0.0,0.0,0.0), 
				Vec3f(0.0, 1.3,-2.5), Vec3f(0,1.1,0), 2, image);
			//render_smoke(y, illum, base_rgb, lit_rgb,radius, opacity, d, light_position,  Vec3f(0.0,0.0,0.0), 
			//	Vec3f(0.0, 8.0,-30.5), Vec3f(0,7,0), 2, image);

			//render_smoke(y, illum, base_rgb, lit_rgb,radius, opacity, d, light_position,  Vec3f(0.0,0.0,0.0), 
			//	Vec3f(0, 3.3,-2.6), Vec3f(0,0.9,0), 2, image);
			write_sgi(image, false, "%s/frame%04d.sgi", filename_format,f);

		}
		



		



	}

	return 0;
}

