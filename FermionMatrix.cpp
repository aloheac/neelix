//
// Created by loheac on 8/28/17.
//

#include "FermionMatrix.h"

// Neelix currently uses the Armadillo C++ Linear Algebra Library to implement matrix operations.
using namespace std;
using namespace arma;

UMatrix::UMatrix( int this_tau, MCParameters this_params, SigmaField* thisSigma ) : tau( this_tau ), params( this_params), ptr_sigma( thisSigma ) {

    U = cx_mat( params.NX, params.NX );  // cx_mat is an Armadillo typedef for Mat< std::complex<double> > (complex dense matrix).
    dU = cx_mat( params.NX, params.NX );
    dU_delta_x = -1;  // -1 denotes dU is invalid.

    if (params.NX != ptr_sigma->NX ) {
        throw AuxiliaryFieldException( "UMatrix: NX dimension of the passed SigmaField does not match the dimension of this UMatrix instance." );
    }

    // Generate kinetic energy matrix T in momentum basis. This is calculated once at construction and then stored
    // for all future calculations. This matrix is defined such that it acts on both sides of the interaction matrix as
    // per the Trotter-Suzuki decomposition; note the additional 1 / 2 in the exponential.
    T = cx_mat( params.NX, params.NX );
    T.zeros();
    int p_i;

    for ( unsigned int i = 0; i < params.NX; i++ ) {
        // Reorganize momentum basis so that it coincides with basis of FFT.
        if ( i <= ( params.NX - 1 ) / 2 ) { p_i = i; } else { p_i = -params.NX + i; }

        const double p2 = pow( 2.0 * p_i * M_PI / (double)params.NX, 2 );
        T( i, i ) = exp( -params.dtau * ( p2 / 2.0 - params.mu ) / 2.0 );
    }

    // Evaluate elements of U. This evaluation will remain valid until the sigma field is modified.
    evaluateElements();
}

void UMatrix::evaluateElements() {
    cx_mat S( params.NX, params.NX );    S.zeros();
    complex<double> A = sqrt( 2.0 * ( complex<double>( exp( params.dtau * params.g ) ) - 1.0 ) );

    // Generate potential energy (interaction) matrix S in coordinate basis.
    for ( unsigned int i = 0; i < params.NX; i++ ) {  // Since S is diagonal, only one loop required.
        S( i, i ) = 1.0 + A * sin( ptr_sigma->get( i, tau ) );
    }

    cx_mat partialProduct;
    complex<double> u_ij;
    cx_mat basis_i( params.NX, 1 );
    cx_mat basis_j( 1, params.NX );

    for ( unsigned int i = 0; i < params.NX; i++ ) {
        basis_i.zeros();
        basis_i( i, 0 ) = 1.0;

        partialProduct = T * fft( basis_i );
        partialProduct = S * ifft( partialProduct );
        partialProduct = T * fft( partialProduct );
        partialProduct = ifft( partialProduct );

        for ( unsigned int j = 0; j < params.NX; j++ ) {
            basis_j.zeros();
            basis_j( 0, j ) = 1.0;

            u_ij = as_scalar( basis_j * partialProduct );  // Retrieve scalar from single-element matrix.
            U( i, j ) = u_ij;
        }
    }

    checksum = ptr_sigma->sum();
}

void UMatrix::evaluateDerivative( int delta_x ) {
    cx_mat S( params.NX, params.NX );    S.zeros();
    complex<double> A = sqrt( 2.0 * ( exp( complex<double>( params.dtau * params.g ) ) - 1.0 ) );

    // Generate potential energy (interaction) matrix S in coordinate basis, where a derivative with respect to sigma
    // has been taken.
    for ( unsigned int i = 0; i < params.NX; i++ ) {  // Since S is diagonal, only one loop required.
        if ( i == delta_x ) {
            S(i, i) = A * cos(ptr_sigma->get(i, tau));  // Derivative of V.
        }
    }

    cx_mat partialProduct;
    complex<double> u_ij;
    cx_mat basis_i( params.NX, 1 );
    cx_mat basis_j( 1, params.NX );

    for ( unsigned int i = 0; i < params.NX; i++ ) {
        basis_i.zeros();
        basis_i( i, 0 ) = 1.0;

        partialProduct = T * fft( basis_i );
        partialProduct = S * ifft( partialProduct );
        partialProduct = T * fft( partialProduct );
        partialProduct = ifft( partialProduct );

        for ( unsigned int j = 0; j < params.NX; j++ ) {
            basis_j.zeros();
            basis_j( 0, j ) = 1.0;

            u_ij = as_scalar( basis_j * partialProduct );  // Retrieve scalar from single-element matrix.
            dU( i, j ) = u_ij;
        }
    }

    dU_delta_x = delta_x;
}

cx_mat UMatrix::getMatrix() {
    if ( checksum != ptr_sigma->sum() ) {
        throw AuxiliaryFieldException( "UMatrix: Matrix requested for the current sigma field which does not match the stored checksum. Matrix must be reevaluated." );
    }

    return U;
}

cx_mat UMatrix::getDerivative( int delta_x ) {
    // If the currently stored derivative is invalid or is with repsect to a different coordinate, evaluate the
    // derivative.
    if ( dU_delta_x != delta_x ) {
        evaluateDerivative( delta_x );
    }

    return dU;
}

void UMatrix::reevaluate() {
    evaluateElements();
    dU_delta_x = -1;  // Any stored derivative is invalid and needs to be recalculated.
}

string UMatrix::to_string() {
    stringstream ss;
    ss << "UMatrix[" << endl;
    ss << U;
    ss << "]";

    return ss.str();
}

FermionMatrix::FermionMatrix( MCParameters this_params, SigmaField* sigma ) : params( this_params ),  ptr_sigma( sigma ) {

    const int SPATIAL_DIMENSION = 1;
    UProduct = vector<UMatrix>();

    for ( unsigned int i = 0; i < params.NTAU; i++ ) {
        UProduct.push_back( UMatrix( i, params, ptr_sigma ) );
    }
}

cx_mat FermionMatrix::getMatrix() {
    cx_mat product( params.NX, params.NX );
    product.eye();

    for ( int i = 0; i < params.NTAU; i++ ) {
        product *= UProduct[ i ].getMatrix();
    }

    return eye( params.NX, params.NX ) + product;
}

cx_mat FermionMatrix::getDerivative( int delta_x, int delta_tau ) {
    cx_mat product( params.NX, params.NX );
    product.eye();

    for ( int i = 0; i < params.NTAU; i++ ) {
        if ( i == delta_tau ) {
            product *= UProduct[ i ].getDerivative( delta_x );
        } else {
            product *= UProduct[ i ].getMatrix();
        }

    }

    return product;
}

void FermionMatrix::reevaluate() {
    for ( int i = 0; i < params.NTAU; i++ ) {
        UProduct[ i ].reevaluate();
    }
}
