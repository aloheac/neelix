//
// Created by loheac on 8/28/17.
//

#include "FermionMatrix.h"

// Neelix currently uses the Armadillo C++ Linear Algebra Library to implement matrix operations.
using namespace std;
using namespace arma;

UMatrix::UMatrix( unsigned int thisNX, unsigned int thisTau, double this_dtau, double this_mu, double this_g, SigmaField* thisSigma )
        : NX( thisNX ), tau( thisTau ), ptr_sigma( thisSigma ), dtau( this_dtau ), mu( this_mu ), g( this_g ) {
    U = cx_mat( NX, NX );  // cx_mat is an Armadillo typedef for Mat< std::complex<double> > (complex dense matrix).
    U.zeros();

    if ( NX != ptr_sigma->NX ) {
        throw AuxiliaryFieldException( "UMatrix: NX dimension of the passed SigmaField does not match the dimension of this UMatrix instance." );
    }

    // Generate kinetic energy matrix T in momentum basis. This is calculated once at construction and then stored
    // for all future calculations.
    T = cx_mat( NX, NX );
    T.zeros();
    int p_i;

    for ( unsigned int i = 0; i < NX; i++ ) {
        // Reorganize momentum basis so that it coincides with basis of FFT.
        if ( i <= ( NX - 1 ) / 2 ) { p_i = i; } else { p_i = -NX + i; }

        const double p2 = pow( 2.0 * p_i * M_PI / (double)NX, 2 );
        T( i, i ) = exp( -dtau * ( p2 / 2.0 - mu ) );
    }

}

void UMatrix::evaluateElements() {
    cx_mat S( NX, NX );    S.zeros();
    complex<double> A = sqrt( 2.0 * ( exp( dtau * g ) - 1.0 ) );

    // Generate potential energy (interaction) matrix S in coordinate basis.
    for ( unsigned int i = 0; i < NX; i++ ) {  // Since S is diagonal, only one loop required.
        S( i, i ) = 1.0 + A * sin( ptr_sigma->get( i, tau ) );
    }

    for ( unsigned int i = 0; i < NX; i++ ) {
        cx_mat basis_i( NX, 1 );
        basis_i.zeros();
        basis_i( i, 0 ) = 1.0;

        for ( unsigned int j = 0; j < NX; j++ ) {
            cx_mat basis_j( 1, NX );
            basis_j.zeros();
            basis_j( 0, j ) = 1.0;

            complex<double> u_ij;
            cx_mat partialProduct;

            partialProduct = ifft( S * basis_i );
            partialProduct =  fft( T * partialProduct );
            partialProduct = basis_j * partialProduct;
            u_ij = as_scalar( partialProduct );  // Retrieve scalar from single-element matrix.

            U( i, j ) = u_ij;
        }
    }
}

void UMatrix::evaluateElementsOfDerivative( int delta_x ) {
    cx_mat S( NX, NX );    S.zeros();
    complex<double> A = sqrt( 2.0 * ( exp( dtau * g ) - 1.0 ) );

    // Generate potential energy (interaction) matrix S in coordinate basis, where a derivative with respect to sigma
    // has been taken.
    for ( unsigned int i = 0; i < NX; i++ ) {  // Since S is diagonal, only one loop required.
        if ( i == delta_x ) {
            S( i, i ) = A * cos( ptr_sigma->get( i, tau ) );  // Derivative of V
        } else {
            S( i, i ) = 1.0 + A * sin( ptr_sigma->get( i, tau ) );
        }

    }

    for ( unsigned int i = 0; i < NX; i++ ) {
        cx_mat basis_i( NX, 1 );
        basis_i.zeros();
        basis_i( i, 0 ) = 1.0;

        for ( unsigned int j = 0; j < NX; j++ ) {
            cx_mat basis_j( 1, NX );
            basis_j.zeros();
            basis_j( 0, j ) = 1.0;

            complex<double> u_ij;
            cx_mat partialProduct;

            partialProduct = ifft( S * basis_i );
            partialProduct =  fft( T * partialProduct );
            partialProduct = basis_j * partialProduct;
            u_ij = as_scalar( partialProduct );  // Retrieve scalar from single-element matrix.

            U( i, j ) = u_ij;

        }
    }
}

cx_mat UMatrix::getMatrix() {
    return U;
}

string UMatrix::to_string() {
    stringstream ss;
    ss << "UMatrix[" << endl;
    ss << U;
    ss << "]";

    return ss.str();
}

FermionMatrix::FermionMatrix( unsigned int this_NX, unsigned int this_NTAU, double this_g, double this_dtau, double this_mu, SigmaField* sigma ) :
        NX( this_NX ), NTAU( this_NTAU ), g( this_g ), dtau( this_dtau ), ptr_sigma( sigma ), mu( this_mu ) {

    const int SPATIAL_DIMENSION = 1;
    UProduct = vector<UMatrix>();

    for ( unsigned int i = 0; i < NTAU; i++ ) {
        UProduct.push_back( UMatrix( NX, i, dtau, mu, g, ptr_sigma ) );
    }
}

cx_mat FermionMatrix::evaluateUProduct() {
    cx_mat product( NX, NX );
    product.eye();

    for ( int i = 0; i < NTAU; i++ ) {
        UProduct[ i ].evaluateElements();
        product *= UProduct[ i ].getMatrix();
    }

    return product;
}

cx_mat FermionMatrix::evaluateUProductDerivative( int delta_x, int delta_tau) {
    cx_mat product( NX, NX );
    product.eye();

    for ( int i = 0; i < NTAU; i++ ) {
        if ( i == delta_tau ) {
            UProduct[ i ].evaluateElementsOfDerivative( delta_x );
        } else {
            UProduct[ i ].evaluateElements();
        }
        product *= UProduct[ i ].getMatrix();
    }

    return product;
}

cx_mat FermionMatrix::evaluate() {
    cx_mat U = evaluateUProduct();
    U = eye( NX, NX ) + U;

    return U;
}

cx_mat FermionMatrix::evaluate_derivative( int delta_x, int delta_tau ) {
    cx_mat deltaU = evaluateUProductDerivative(delta_x, delta_tau);
    return deltaU;
}
