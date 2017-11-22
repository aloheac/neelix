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

    cx_mat pi_dot( params.NX, params.NTAU);

    for ( int i = 0; i < params.NX; i++ ) {
        for ( int j = 0; j < params.NTAU; j++ ) {
            cx_mat deltaM = matM.getDerivative(i, j);
            pi_dot( i, j ) =  2.0 * as_scalar( trace( Minv * deltaM ) );
        }
    }
    cout << "        | Action: " << 2.0 * log( det( M ) ) << endl;
    return pi_dot;
}

void HMCEvolver::integrateSigma() {
    if ( sample_count == 20 ) {
        cout << "        | New momentum field initialized." << endl;
        sample_count = 0;
        pi->initialize();
    }
    sample_count++;

    FermionMatrix M( params, sigma );
    complex<double> delta_sigma, delta_pi;
    double S_initial = -2.0 * log( det( M.getMatrix() ) ).real() + 0.5 * pi->sum().real();

    cx_mat pi_dot_0 = calculatePiDot();

   // Update sigma field.
   for ( int i = 0; i < params.NX; i++ ) {
       for ( int j = 0; j < params.NTAU; j++ ) {
           delta_sigma = pi->get( i, j ) * params.dt + 0.5 * pi_dot_0( i, j ) * pow( params.dt, 2 );
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

   M.reevaluate();
   double S_final = -2.0 * log( det( M.getMatrix() ) ).real() + 0.5 * pi->sum().real();

   // Perform Metropolis accept-reject.
   double r = rand_metropolis( rand_generator );
   if ( exp( -S_final + S_initial ) > r ) {
       cout << "        | Accept. ( exp[ dH ] = " << exp( -S_final + S_initial ) << ", r = " << r <<  " )" << endl;
   } else {
       cout << "        | Reject. ( exp[ dH ] = " << exp( -S_final + S_initial ) << ", r = " << r <<  " )" << endl;
   }
}