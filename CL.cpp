//
// Created by loheac on 10/24/17.
//

#include "CL.h"
#include "FermionMatrix.h"

using namespace std;
using namespace arma;

CLEvolver::CLEvolver( MCParameters this_params, SigmaField* this_sigma ) : params( this_params ), sigma( this_sigma ) { }

cx_mat CLEvolver::calculateSigmaDot() {
    // Calculate the fermion matrix and its inverse for the current sigma field.
    FermionMatrix matM( params.NX, params.NTAU, params.g, params.dtau, params.mu, sigma );
    cx_mat M = matM.evaluate();
    cx_mat Minv = inv( M );

    cx_mat sigma_dot( params.NX, params.NTAU);    sigma_dot.zeros();
    complex<double> deltaS;
    for ( int i = 0; i < params.NX; i++ ) {
        for ( int j = 0; j < params.NTAU; j++ ) {
            cx_mat deltaM = matM.evaluate_derivative(i, j);
            deltaS = 2.0 * as_scalar( trace( Minv * deltaM ) );
            sigma_dot( i, j ) = deltaS;
        }
    }

    return sigma_dot;
}

void CLEvolver::integrateSigma() {
    cx_mat sigma_dot = calculateSigmaDot();
    complex<double> current;

    for ( int i = 0; i < params.NX; i++ ) {
        for (int j = 0; j < params.NTAU; j++) {
            current = sigma->get( i, j );
            sigma->set( i, j, current + sigma_dot( i, j ) * params.dt );

        }
    }
}