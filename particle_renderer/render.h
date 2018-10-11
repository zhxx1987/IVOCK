#ifndef RENDER_H
#define RENDER_H

#include "array1.h"
#include "array2.h"
#include "vec.h"

void
render_smoke(Array1<Vec3f> const &x,
			 Array1<float> const &illum,
			 Array1<Vec3f> const &base_rgb, // particle colour with just ambient illumination
			 Array1<Vec3f> const &lit_rgb, // particle colour in full lighting
			 float radius,
			 float opacity,
			 Array1<float> const &den, 
			 Vec3f const &light_position,
			 Vec3f const &background_rgb,
			 Vec3f const &camera_position,
			 Vec3f const &camera_target, // position we are looking at (not a direction vector!)
			 float const camera_perspective, // 3 is reasonable, use higher numbers to zoom in
			 Array2<Vec3f> &image); // make sure to set the image size you want before calling render_smoke()

#endif
