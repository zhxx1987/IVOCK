#ifndef WRITE_SGI_H
#define WRITE_SGI_H

#include "array2.h"
#include "vec.h"

// Monochrome
// Returns false if write fails.
bool write_sgi(const Array2f &img, bool high_precision, const char *filename_format, ...);

// RGB
// Returns false if write fails.
bool write_sgi(const Array2<Vec3f> &img, bool high_precision, const char *filename_format, ...);

#endif
