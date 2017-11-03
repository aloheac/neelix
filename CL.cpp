//
// Created by loheac on 10/24/17.
//

#include "CL.h"
#include "FermionMatrix.h"
#include "Observables.h"

using namespace std;
using namespace arma;

CLEvolver::CLEvolver( MCParameters this_params, SigmaField* this_sigma ) : params( this_params ), sigma( this_sigma ) {
    rand_generator = default_random_engine( 23412341 );
    rand_normal = normal_distribution<double>( 0.0, 2.0 );
    rand_uniform = uniform_real_distribution<double>( -M_PI, M_PI );
}

cx_mat CLEvolver::calculateSigmaDot() {
    // Calculate the fermion matrix and its inverse for the current sigma field.
    FermionMatrix matM( params.NX, params.NTAU, params.g, params.dtau, params.mu, sigma );
    cx_mat M = matM.evaluate();
    cx_mat Minv = inv( M );

    cx_mat sigma_dot( params.NX, params.NTAU);    sigma_dot.zeros();
    complex<double> deltaS;
    for ( int i = 0; i < params.NX; i++ ) {
#pragma omp parallel for shared( M, Minv, sigma_dot ) private( deltaS )
        for ( int j = 0; j < params.NTAU; j++ ) {
            cx_mat deltaM = matM.evaluate_derivative(i, j);
            deltaS = 2.0 * as_scalar( trace( Minv * deltaM ) );
            sigma_dot( i, j ) = deltaS;
        }
    }

    complex<double> action = 2.0 * log( det( M ) );
    double magnitude = sqrt( action.real() * action.real() + action.imag() * action.imag() );
    cout << "         |  Action: " << magnitude << "    Cos: " << action.real() / magnitude << "    Sin: " << action.imag() / magnitude << endl;

    return sigma_dot;
}

void CLEvolver::integrateSigma() {
    cx_mat sigma_dot = calculateSigmaDot();
    complex<double> current;
    complex<double> action, regulator, noise, total;
    vector<double> vec_action_re;
    vector<double> vec_action_im;
    vector<double> vec_stn;

    for ( int i = 0; i < params.NX; i++ ) {
        for (int j = 0; j < params.NTAU; j++) {
            current = sigma->get( i, j );
            noise = complex<double>( rand_normal( rand_generator ), 0.0 ) * sqrt( params.dt );
            regulator = complex<double>( -2.0 * params.xi * sigma->get( i,j ).real() * params.dt, -2.0 * params.xi * sigma->get( i, j ).imag() * params.dt );
            action = sigma_dot( i, j ) * params.dt;
            total =  action + regulator + noise;

            vec_stn.push_back( ((action + regulator) / noise).real() );
            vec_action_re.push_back( total.real() );
            vec_action_im.push_back( total.imag() );

            sigma->set( i, j, current + total );

        }
    }

    cout << "         |  StN  -- Max: " << max( vec_stn ) << "    Min: " << min( vec_stn ) << "    Mean: " << mean( vec_stn ) << endl;
    cout << "         |  Real -- Max: " << max( vec_action_re ) << "    Min: " << min( vec_action_re ) << "    Mean: " << mean( vec_action_re ) << endl;
    cout << "         |  Imag -- Max: " << max( vec_action_im ) << "    Min: " << min( vec_action_im ) << "    Mean: " << mean( vec_action_im ) << endl;
    if ( max( vec_action_re ) > 0.5 or max( vec_action_im ) > 0.2 ) { params.dt /= 2; }
    else if ( min( vec_action_re ) < 0.002 ) params.dt *= 2;
}