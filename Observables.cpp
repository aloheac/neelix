//
// Created by loheac on 10/24/17.
//

#include "CL.h"
#include <math.h>
#include "Observables.h"
#include "FermionMatrix.h"

using namespace std;
using namespace arma;

double freeGasDensity( int NX, double BETA, double MU ) {
    double freeGasDensity = 0.0;
    double p2, expFactor;
    int p_i;

    for (int i = 0; i < NX; i++) {
        p_i = i - ( NX - 1 ) / 2;
        p2 = pow(2.0 * p_i * M_PI / NX, 2);
        expFactor = exp(-BETA * (0.5 * p2 - MU));
        freeGasDensity += 2.0 * expFactor / (1.0 + expFactor);
    }

    return freeGasDensity;

}

complex<double> calculate_density( MCParameters params, SigmaField* sigma ) {
    FermionMatrix matM( params.NX, params.NTAU, params.g, params.dtau, params.mu, sigma );
    cx_mat M = matM.evaluate();
    cx_mat U = M - eye( params.NX, params.NX );
    cx_mat invM = inv( M );

    complex<double> n = 2.0 * as_scalar( trace( invM * U ) );
    return n / complex<double>( freeGasDensity( params.NX, params.NTAU * params.dtau, params.mu ), 0.0 );
}