//
// Created by loheac on 10/24/17.
//

#include "CL.h"
#include "FermionMatrix.h"

using namespace std;
using namespace arma;

CLEvolver::CLEvolver( int thisNX, int thisNTAU, double this_dt, SigmaField* this_sigma ) : NX( thisNX ),
                                                                                           NTAU( thisNTAU ),
                                                                                           dt( this_dt ),
                                                                                           sigma( this_sigma ) { }

cx_mat CLEvolver::calculateSigmaDot() {
    // Calculate the fermion matrix.
    double dtau = 0.05;
    double g = 0.35;
    double mu = 0.0;
    complex<double> deltaS;
    FermionMatrix matM( NX, NTAU, g, dtau, mu, sigma );
    cx_mat M = matM.evaluate();
    cx_mat Minv = inv( M );
    cx_mat sigma_dot(NX, NTAU);    sigma_dot.zeros();

    for ( int i = 0; i < NX; i++ ) {
        for ( int j = 0; j < NTAU; j++ ) {
            cx_mat deltaM = matM.evaluate_derivative(i, j);
            deltaS = 2.0 * as_scalar( trace( Minv * deltaM ) );
            sigma_dot( i, j ) = deltaS;
        }
    }


    return sigma_dot;
}