#ifndef SHADOW_H
#define SHADOW_H

#include "array1.h"
#include "array2.h"
#include "vec.h"

// For now, light_position has to be well outside the bounding box of the points - will return false if it's too close or inside.
// The output "illum" is a scalar in [0,1] for each particle indicating how much light arrives at it (0=dark, 1=fully illuminated)
bool
compute_shadows(Array1<Vec3f> const &x,
                float radius,
                float opacity,
				Array1<float> const &dens,
                Vec3f const &light_position,
                float quality_factor,
                Array1<float> &illum,
                Array2f &shadow_map); // not necessary, but you might want to see it for debugging

#endif
