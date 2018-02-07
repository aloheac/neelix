//
// Created by loheac on 2/6/18.
//

#ifndef NEELIX_FFTWINTERFACE_H
#define NEELIX_FFTWINTERFACE_H

#include <armadillo>
#include <fftw3.h>

class FFTWInterface {
public:
    FFTWInterface( int NX );

    arma::cx_mat fftw_forward( arma::cx_mat in );

    arma::cx_mat fftw_inverse( arma::cx_mat in );

private:
    int transformSize;

    fftw_plan p_forward;

    fftw_plan p_inverse;

    fftw_complex *in, *out;
};

#endif //NEELIX_FFTWINTERFACE_H
