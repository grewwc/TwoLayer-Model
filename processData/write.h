//
// Created by grewwc on 27/12/2017.
//

#ifndef TWOLAYER_WRITE_H
#define TWOLAYER_WRITE_H

#include <fstream>

inline void write_data(double *x, double *y, int N, const char *fname = "/Users/grewwc/test/test.dat")
{
    std::ofstream fout{fname, std::ios_base::out};
    for (int i = 0; i < N; ++i)
    {
        fout << x[i] << " " << y[i] << std::endl;
    }

    fout.close();
}


#endif //TWOLAYER_WRITE_H
