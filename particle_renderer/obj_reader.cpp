#include "obj_reader.h"
#include<iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

void load_obj(const char* filename, vector<Vec3f> & vertices, vector<Vec3ui> & faces)
{
	if (vertices.size()>0)
	{
		vertices.resize(0);
	}
	if(faces.size()>0)
	{
		faces.resize(0);
	}
	std::ifstream ifs(filename, std::ifstream::in);
	if (!ifs) { cerr << "Cannot open " << filename << endl; exit(1); }


	string line;
	while (getline(ifs, line)) {
		if (line.substr(0,2) == "v ") {
			std::istringstream s(line.substr(2));
			Vec3f v; s>>v[0]; s>>v[1]; s>>v[2];
			vertices.push_back(v);
		}  else if (line.substr(0,2) == "f ") {
			std::istringstream s(line.substr(2));
			unsigned int a,b,c;
			s>>a; s>>b; s>>c;
			a--; b--; c--;
			Vec3ui f; f[0]=a; f[1]=b; f[2]=c;
			faces.push_back(f);
		}
		else if (line[0] == '#') { /* ignoring this line */ }
		else { /* ignoring this line */ }

	}


}

void normalize_object(vector<Vec3f> & vertices, vector<Vec3ui> & faces, float scale)
{
	if (vertices.size()>0 && faces.size()>0)
	{
		Vec3f v = Vec3f(0,0,0);
		for (int i=0;i<vertices.size();i++)
		{
			v += vertices[i];
		}
		v /= (float)(vertices.size());

		for (int i=0; i<vertices.size();i++)
		{
			vertices[i] -= v;
		}

		//find bounding box
		float xmin=vertices[0].v[0], xmax=vertices[0].v[0], ymin=vertices[0].v[1],ymax=vertices[0].v[1],zmin=vertices[0].v[2],zmax=vertices[0].v[2];
		for (int i=0;i<vertices.size(); i++)
		{
			if(vertices[i].v[0]<xmin)
			{
				xmin = vertices[i].v[0];
			}
			if(vertices[i].v[0]>xmax)
			{
				xmax = vertices[i].v[0];
			}
			if(vertices[i].v[1]<ymin)
			{
				ymin = vertices[i].v[1];
			}
			if(vertices[i].v[1]>ymax)
			{
				ymax = vertices[i].v[1];
			}
			if(vertices[i].v[2]<zmin)
			{
				zmin = vertices[i].v[2];
			}
			if(vertices[i].v[2]>zmax)
			{
				zmax = vertices[i].v[2];
			}
		}

		//compute scale of the object
		float zl = zmax - zmin;
		float yl = ymax - ymin;
		float xl = xmax - xmin;

		//do scaling so that the longest axis equals the scale factor
		float l = std::max(zl, std::max(yl, xl));

		float factor = scale/l;

		for (int i=0; i<vertices.size(); i++)
		{
			vertices[i] *= factor;
		}

	}
}