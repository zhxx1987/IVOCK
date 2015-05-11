#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include "read_particles.h"
#include "bfstream.h"

bool
read_particles(Array1<Vec3f> &x,
			   Array1<float> &d,
               const char *filename_format,
               ...)
{
    // open the file

    va_list ap;
    va_start(ap, filename_format);
#ifdef WIN32
    int len=_vscprintf(filename_format, ap) // _vscprintf doesn't count
                                        +1; // terminating '\0'
    char *filename=new char[len];
    vsprintf(filename, filename_format, ap);
	//sprintf(filename,"C:\\Users\\xinxin\\Desktop\\raising smoke surface data\\Particle_data%4d.bin", );
    FILE *in=std::fopen(filename, "rb");
	if(in==NULL) printf("file error\n");
    delete[] filename;
    va_end(ap);
#else
    char *filename;
    vasprintf(&filename, filename_format, ap);
    FILE *in=std::fopen(filename, "rb");
    std::free(filename);
    va_end(ap);
#endif

    if(!in) return false;

    // figure out length of file

    if(std::fseek(in, 0, SEEK_END)){
        std::fclose(in);
        return false;
    }

    long int filesize=std::ftell(in);
    if(filesize<0){
        std::fclose(in);
        return false;
    }

    if(std::fseek(in, 0, SEEK_SET)){
        std::fclose(in);
        return false;
    }

    // read the particles

    x.resize(filesize/(4*sizeof(float)));
	d.resize(filesize/(4*sizeof(float)));
    float zero=0;
    for(int i=0; i<x.n; ++i){
        // read x y z for particle i
        std::fread(&x[i], sizeof(float), 3, in);
        // skip 0.f
		x[i][0] = -x[i][0];
		//float z = x[i][2];
		x[i][2] = -x[i][2];
		//x[i][0] = -z;
        //std::fread(&zero, sizeof(float), 1, in);
		std::fread(&d[i],sizeof(float), 1, in);
    }

    std::fclose(in);
    return true;
}
