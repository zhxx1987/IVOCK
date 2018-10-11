#ifndef __obj_reader__
#define __obj_reader__
#include <stdio.h>
#include "vec.h"
#include <iostream>
#include <string>
using namespace std;





void load_obj(const char* filename, vector<Vec3f> & vertices, vector<Vec3ui> & faces);
void normalize_object(vector<Vec3f> & vertices, vector<Vec3ui> & faces, float scale);







#endif