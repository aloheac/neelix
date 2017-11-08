//
// Created by loheac on 9/3/17.
//

#include "HMC.h"
#include "FermionMatrix.h"

using namespace std;
using namespace arma;

HMCEvolver::HMCEvolver( MCParameters this_params, SigmaField* this_sigma, MomentumField* this_pi ) : params( this_params ), sigma( this_sigma ), pi( this_pi ) {
    sample_count = 0;
}

cx_mat HMCEvolver::calculatePiDot() {
    // Calculate the fermion matrix and its inverse for the current sigma field.
    FermionMatrix matM( params.NX, params.NTAU, params.g, params.dtau, params.mu, sigma );
    cx_mat M = matM.evaluate();
    cx_mat Minv = inv( M );

    cx_mat pi_dot( params.NX, params.NTAU);    pi_dot.zeros();
    complex<double> deltaS;
    for ( int i = 0; i < params.NX; i++ ) {
        for ( int j = 0; j < params.NTAU; j++ ) {
            cx_mat deltaM = matM.evaluate_derivative(i, j);
            deltaS = 2.0 * as_scalar( trace( Minv * deltaM ) );
            pi_dot( i, j ) = deltaS;
        }
    }
    cout << "        | Action: " << 2.0 * log( det( M ) )  + 0.5 * pi->sum() << endl;
    return pi_dot;
}

void HMCEvolver::integrateSigma() {
    if ( sample_count > 10 ) {
        cout << "        | New momentum field initialized." << endl;
        sample_count = 0;
        pi->initialize();
    }

    sample_count++;

    complex<double> current_sigma, current_pi, delta_sigma, delta_pi;

    cx_mat pi_dot_0 = calculatePiDot();

    // Update sigma field.
    for ( int i = 0; i < params.NX; i++ ) {
        for ( int j = 0; j < params.NTAU; j++ ) {
            delta_sigma = pi->get( i, j ) * params.dt + 0.5 * pi_dot_0( i, j ) * ( params.dt * params.dt );
            sigma->set( i, j, sigma->get( i, j ) + delta_sigma );
        }
    }

    cx_mat pi_dot_1 = calculatePiDot();

    // Update momentum field.
    for ( int i = 0; i < params.NX; i++ ) {
        for ( int j = 0; j < params.NTAU; j++ ) {
            delta_pi = 0.5 * ( pi_dot_0( i, j ) + pi_dot_1( i, j ) ) * params.dt;
            pi->set( i, j, pi->get( i, j ) + delta_pi );
        }
    }
}