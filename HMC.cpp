//
// Created by loheac on 9/3/17.
//

#include "HMC.h"
#include "FermionMatrix.h"

using namespace std;
using namespace arma;

HMCEvolver::HMCEvolver( MCParameters this_params, SigmaField* this_sigma, MomentumField* this_pi ) : params( this_params ), sigma( this_sigma ), pi( this_pi ) {
    sample_count = 0;
    rand_generator = default_random_engine( 23412341 );
    rand_metropolis = uniform_real_distribution<double>( 0.0, 1.0 );
}

cx_mat HMCEvolver::calculatePiDot() {
    // Calculate the fermion matrix and its inverse for the current sigma field.
    FermionMatrix matM( params, sigma );
    cx_mat M = matM.getMatrix();
    cx_mat Minv = inv( M );

    cx_mat pi_dot( params.NX, params.NTAU);    pi_dot.zeros();
    complex<double> deltaS;
#pragma omp parallel for shared( matM, M, Minv ) private( deltaS )
    for ( int i = 0; i < params.NX; i++ ) {
        for ( int j = 0; j < params.NTAU; j++ ) {
            cx_mat deltaM = matM.getDerivative(i, j);
            deltaS = 2.0 * as_scalar( trace( Minv * deltaM ) );
            pi_dot( i, j ) = deltaS;
        }
    }
    cout << "        | Action: " << 2.0 * log( det( M ) )  + 0.5 * pi->sum() << endl;
    return pi_dot;
}

void HMCEvolver::integrateSigma() {
    if ( sample_count == 20 ) {
        cout << "        | New momentum field initialized." << endl;
        sample_count = 0;
        pi->initialize();
    }

    sample_count++;

    complex<double> delta_sigma, delta_pi;
    FermionMatrix matM( params, sigma );
    cx_mat M = matM.getMatrix();
    double S_initial = 2.0 * log( det( M ) ).real() + 0.5 * pi->sum().real();

    // Update sigma field.
    for ( int i = 0; i < params.NX; i++ ) {
        for ( int j = 0; j < params.NTAU; j++ ) {
            delta_sigma = pi->get( i, j ) * params.dt;
            sigma->set( i, j, sigma->get( i, j ) + delta_sigma );
        }
    }

    cx_mat pi_dot = calculatePiDot();

    // Update momentum field.
    for ( int i = 0; i < params.NX; i++ ) {
        for ( int j = 0; j < params.NTAU; j++ ) {
            delta_pi = pi_dot( i, j ) * params.dt;
            pi->set( i, j, pi->get( i, j ) + delta_pi );
        }
    }

    matM = FermionMatrix( params, sigma );
    M = matM.getMatrix();



    double S_final = 2.0 * log( det( M ) ).real() + 0.5 * pi->sum().real();

    // Perform Metropolis accept-reject.
    double r = rand_metropolis( rand_generator );
    if ( exp( -S_final + S_initial ) > r ) {
        cout << "        | Accept. ( exp[ dH ] = " << exp( -S_final + S_initial ) << ", r = " << r <<  " )" << endl;
    } else {
        cout << "        | Reject. ( exp[ dH ] = " << exp( -S_final + S_initial ) << ", r = " << r <<  " )" << endl;
    }
}