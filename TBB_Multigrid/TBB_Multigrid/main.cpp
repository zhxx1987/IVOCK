#include <cmath>
#include "tbb/tbb.h"

double *output;
double *input;


int main() {
    const int size = 512*512*512;
    output = new double[size];
    input = new double[size];
   
    tbb::parallel_for(0, size, 1, [=](int i) {
        input[i] = i;
    });
      
    tbb::parallel_for(0, size, 1, [=](int i) {

        output[i] = sqrt(sin(input[i])*sin(input[i]) + cos(input[i])*cos(input[i]));
            
   });
    delete[] input;
    delete[] output;
    return 0;
}