#ifndef READ_PARTICLES_H
#define READ_PARTICLES_H

#include "array1.h"
#include "vec.h"

bool
read_particles(Array1<Vec3f> &x,
			   Array1<float> &d,
               const char *filename_format,
               ...);

#endif
