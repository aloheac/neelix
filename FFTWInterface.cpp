//
// Created by loheac on 2/6/18.
//

#include "FFTWInterface.h"

using namespace std;
using namespace arma;

void cast_arma_mat_to_complex( fftw_complex* in, cx_mat* M, int NX ) {
    complex<double> v;
    for ( unsigned int i = 0; i < NX; i++ ) {
        v = as_scalar( (*M)( i, 0 ) );
        in[ i ][ 0 ] = v.real();
        in[ i ][ 1 ] = v.imag();
    }
}

void cast_complex_to_arma_mat( fftw_complex* out, cx_mat* M, int NX ) {
    for ( unsigned int i = 0; i < NX; i++ ) {
        (*M)( i, 0 ) = complex<double>( out[ i ][ 0 ], out[ i ][ 1 ] );
    }
}

FFTWInterface::FFTWInterface( int NX ) : transformSize( NX ) {
    in = (fftw_complex*) fftw_malloc( sizeof( fftw_complex ) * NX );
    out = (fftw_complex*) fftw_malloc( sizeof( fftw_complex ) * NX );

    p_forward = fftw_plan_dft_1d( NX, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
    p_inverse = fftw_plan_dft_1d( NX, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
}

cx_mat FFTWInterface::fftw_forward( cx_mat mat_in ) {
    cx_mat mat_out( transformSize, 1 );

    cast_arma_mat_to_complex( in, &mat_in, transformSize );
    fftw_execute( p_forward );
    cast_complex_to_arma_mat( out, &mat_out, transformSize );

    return mat_out;
}

cx_mat FFTWInterface::fftw_inverse( cx_mat mat_in ) {
    cx_mat mat_out( transformSize, 1 );

    cast_arma_mat_to_complex( in, &mat_in, transformSize );
    fftw_execute( p_forward );
    cast_complex_to_arma_mat( out, &mat_out, transformSize );

    return mat_out / transformSize;
}