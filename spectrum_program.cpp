#include <iostream>
using namespace std;
#include "spectrum.h"

int main()
{
    isignal<int> first(1);
    vector <int> myvec = {2,2,2,2,2,2,2,2}; //number of elements must be 2^n

    first.save_vector(myvec);
    first.print_vector();

    first.FFTAnalysis();
    first.print_spectrum();

    return 0;
}
