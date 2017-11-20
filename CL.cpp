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
    rand_normal = normal_distribution<double>( 0.0, sqrt( 2.0 ) );
}

cx_mat CLEvolver::calculateSigmaDot() {
    // Calculate the fermion matrix and its inverse for the current sigma field.
    FermionMatrix matM( params.NX, params.NTAU, params.g, params.dtau, params.mu, sigma );
    cx_mat M = matM.getMatrix();
    cx_mat Minv = inv( M );

    cx_mat sigma_dot( params.NX, params.NTAU);    sigma_dot.zeros();
    complex<double> deltaS;
    for ( int i = 0; i < params.NX; i++ ) {
#pragma omp parallel for shared( M, Minv, sigma_dot ) private( deltaS )
        for ( int j = 0; j < params.NTAU; j++ ) {
            cx_mat deltaM = matM.getDerivative(i, j);
            deltaS = 2.0 * as_scalar( trace( Minv * deltaM ) );
            sigma_dot( i, j ) = deltaS;
        }
    }

    complex<double> action = 2.0 * log( det( M ) );
    double magnitude = sqrt( action.real() * action.real() + action.imag() * action.imag() );
    cout << "         |  Action: " << magnitude << "    Cos: " << action.real() / magnitude << "    Sin: " << action.imag() / magnitude << endl;

    return sigma_dot;
}

cx_mat calculateAuxSigmaDot( MCParameters params, SigmaField* sigma ) {
    // Calculate the fermion matrix and its inverse for the current sigma field.
    FermionMatrix matM( params.NX, params.NTAU, params.g, params.dtau, params.mu, sigma );
    cx_mat M = matM.getMatrix();
    cx_mat Minv = inv( M );

    cx_mat sigma_dot( params.NX, params.NTAU);    sigma_dot.zeros();
    complex<double> deltaS;
    for ( int i = 0; i < params.NX; i++ ) {
#pragma omp parallel for shared( M, Minv, sigma_dot ) private( deltaS )
        for ( int j = 0; j < params.NTAU; j++ ) {
            cx_mat deltaM = matM.getDerivative(i, j);
            deltaS = 2.0 * as_scalar( trace( Minv * deltaM ) );
            sigma_dot( i, j ) = deltaS;
        }
    }

    return sigma_dot;
}

void CLEvolver::integrateSigma() {
    cx_mat sigma_dot = calculateSigmaDot();
    complex<double> current;
    complex<double> action, regulator, noise, total;
    vector<double> vec_action_re;
    vector<double> vec_action_im;
    vector<double> vec_stn;
    int METHOD = 1;

    if ( METHOD == 0 ) {  // Euler method.
        for ( int i = 0; i < params.NX; i++ ) {
            for (int j = 0; j < params.NTAU; j++) {
                current = sigma->get( i, j );
                noise = complex<double>( rand_normal( rand_generator ), 0.0 ) * sqrt( params.dt );
                regulator = complex<double>( 2.0 * params.xi * sigma->get( i,j ).real() * params.dt, 2.0 * params.xi * sigma->get( i, j ).imag() * params.dt );
                action = sigma_dot( i, j ) * params.dt;
                total =  action + regulator + noise;

                vec_stn.push_back( ( (action + regulator ) / noise).real() );
                vec_action_re.push_back( total.real() );
                vec_action_im.push_back( total.imag() );

                sigma->set( i, j, current + total );
            }
        }
    } else {
        SigmaField sigma_prime( 1, params.NX, params.NTAU );
        cx_mat k1, k2, k3, k4;
        complex<double> delta;

        noise = complex<double>( rand_normal( rand_generator ), 0.0 );  // Need to carefully consider sqrt(dt) integration step when computing RK4 coefficients.

        k1 = sigma_dot;

        // Compute auxiliary field sigma_prime for RK coefficient k2; determined as sigma_n + dt/2 * k1.
        for ( int i = 0; i < params.NX; i++ ) {
            for ( int j = 0; j < params.NTAU; j++ ) {
                current = sigma->get( i, j );
                regulator = complex<double>( -2.0 * params.xi * sigma->get( i, j ).real(), -2.0 * params.xi * sigma->get( i, j ).imag() );
                total = ( k1( i, j ) + regulator ) * params.dt / 2.0 + noise * sqrt( params.dt ) / 2.0;
                sigma_prime.set( i, j, current + total );
            }
        }

        k2 = calculateAuxSigmaDot( params, &sigma_prime );

        // Compute auxiliary field sigma_prime for RK coefficient k3; determined as sigma_n + dt/2 * k2.
        for ( int i = 0; i < params.NX; i++ ) {
            for ( int j = 0; j < params.NTAU; j++ ) {
                current = sigma->get( i, j );
                regulator = complex<double>( -2.0 * params.xi * sigma->get( i, j ).real(), -2.0 * params.xi * sigma->get( i, j ).imag() );
                total = ( k2( i, j ) + regulator ) * params.dt / 2.0  + noise * sqrt( params.dt ) / 2.0;;
                sigma_prime.set( i, j, current + total );
            }
        }

        k3 = calculateAuxSigmaDot( params, &sigma_prime );

        // Compute auxiliary field sigma_prime for RK coefficient k4; determined as sigma_n + dt * k3.
        for ( int i = 0; i < params.NX; i++ ) {
            for ( int j = 0; j < params.NTAU; j++ ) {
                current = sigma->get( i, j );
                regulator = complex<double>( -2.0 * params.xi * sigma->get( i, j ).real(), -2.0 * params.xi * sigma->get( i, j ).imag() );
                total = ( k3( i, j ) + regulator ) * params.dt + noise * sqrt( params.dt );
                sigma_prime.set( i, j, current + total );
            }
        }

        k4 = calculateAuxSigmaDot( params, &sigma_prime );

        // Update auxiliary field sigma according to RK4.
        for ( int i = 0; i < params.NX; i++ ) {
            for ( int j = 0; j < params.NTAU; j++ ) {
                current = sigma->get( i, j );
                delta = ( params.dt / 6.0 ) * ( k1( i, j ) + 2.0 * k2( i, j ) + 2.0 * k3( i, j ) + k4( i, j ) );
                sigma->set( i, j, current + delta );

                vec_stn.push_back( ( delta / noise ).real() );
                vec_action_re.push_back( delta.real() );
                vec_action_im.push_back( delta.imag() );
            }
        }
    }


    cout << "         |  StN  -- Max: " << max( vec_stn ) << "    Min: " << min( vec_stn ) << "    Mean: " << mean( vec_stn ) << endl;
    cout << "         |  Real -- Max: " << max( vec_action_re ) << "    Min: " << min( vec_action_re ) << "    Mean: " << mean( vec_action_re ) << endl;
    cout << "         |  Imag -- Max: " << max( vec_action_im ) << "    Min: " << min( vec_action_im ) << "    Mean: " << mean( vec_action_im ) << endl;
    if ( max( vec_action_re ) > 0.001 or max( vec_action_im ) > 0.001 ) { params.dt /= 2; }
    else if ( min( vec_action_re ) < 0.0001 ) params.dt *= 2;
}