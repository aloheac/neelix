//
// Created by loheac on 8/28/17.
//

#include "FermionMatrix.h"

// Neelix currently uses the Armadillo C++ Linear Algebra Library to implement matrix operations.
using namespace std;
using namespace arma;

// Define static matrices T and T_x. These defaults will be overwritten when UMatrix::initialize() is called.
cx_mat UMatrix::T, UMatrix::T_x;

UMatrix::UMatrix( int this_tau, MCParameters this_params, SigmaField* thisSigma ) : tau( this_tau ), params( this_params), ptr_sigma( thisSigma ) {
    U = cx_mat( params.NX, params.NX );  // cx_mat is an Armadillo typedef for Mat< std::complex<double> > (complex dense matrix).
    dU = new cx_mat[ params.NX ];  // Initialize array of matrix derivatives. Array deleted in destructor.

    if (params.NX != ptr_sigma->NX ) {
        throw AuxiliaryFieldException( "UMatrix: NX dimension of the passed SigmaField does not match the dimension of this UMatrix instance." );
    }

    // Evaluate elements of U. This evaluation will remain valid until the sigma field is modified.
    evaluateElements();

    // Evaluate derivatives of U for each spatial point in the lattice. This evaluation will remain valid until the
    // sigma field is modified.
    for ( int i = 0; i < params.NX; i++ ) {
        evaluateDerivative( i );
    }
}

UMatrix::UMatrix( const UMatrix &obj ) {
    // Copy primitive data fields.
    tau = obj.tau;
    params = obj.params;
    ptr_sigma = obj.ptr_sigma;
    checksum = obj.checksum;

    // Copy matrices using arma::mat copy constructor.
    U = cx_mat( obj.U );

    dU = new cx_mat[ params.NX ];  // Array of matrices; copy each element.
    for ( int i = 0; i < params.NX; i++ ) {
        dU[ i ] = cx_mat( obj.dU[ i ] );
    }
}

UMatrix::~UMatrix() {
    delete[] dU;
}

void UMatrix::evaluateElements() {
    cx_mat S( params.NX, params.NX );    S.zeros();
    complex<double> A = sqrt( 2.0 * ( complex<double>( exp( params.dtau * params.g ) ) - 1.0 ) );

    // Generate potential energy (interaction) matrix S in coordinate basis.
    for ( unsigned int i = 0; i < params.NX; i++ ) {  // Since S is diagonal, only one loop required.
        S( i, i ) = 1.0 + A * sin( ptr_sigma->get( i, tau ) );
    }

    U = T_x * S * T_x;

    if ( params.ENABLE_FIELD_CHECKSUM ) checksum = ptr_sigma->sum();
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

    dU[ delta_x ] = T_x * S * T_x;
}

cx_mat UMatrix::getMatrix() {
    if ( params.ENABLE_FIELD_CHECKSUM and checksum != ptr_sigma->sum() ) {
        throw AuxiliaryFieldException( "UMatrix: Matrix requested for the current sigma field which does not match the stored checksum. Matrix must be reevaluated." );
    }

    return U;
}

cx_mat UMatrix::getDerivative( int delta_x ) {
    return dU[ delta_x ];
}

void UMatrix::reevaluate() {
    // Evaluate matrix elements of U.
    evaluateElements();

    // Evaluate derivatives of U with respect to each spatial lattice point.
    for ( int i = 0; i < params.NX; i++ ) {
        evaluateDerivative( i );
    }
}

void UMatrix::initialize( MCParameters params ) {
    // Generate kinetic energy matrix T in momentum basis. This is calculated once at initialization and then stored
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

    // Generate transformation of the kinetic energy matrix T to the coordinate basis under a FFT, producing T_x.
    T_x = cx_mat( params.NX, params.NX );
    cx_mat basis_i( params.NX, 1 );
    cx_mat basis_j( 1, params.NX );

    for ( int i = 0; i < params.NX; i++ ) {
        basis_i.zeros();
        basis_i( i, 0 ) = 1.0;
        for ( int j = 0; j < params.NX; j++ ) {
            basis_j.zeros();
            basis_j( 0, j ) = 1.0;
            T_x( i, j ) = as_scalar( basis_j * ifft( T * fft( basis_i ) ) );
        }
    }
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

    // Note that the elements of UMatrix are evaluated upon construction.
    for ( unsigned int i = 0; i < params.NTAU; i++ ) {
        UProduct.push_back( UMatrix( i, params, ptr_sigma ) );
    }

    // If the use of partial products is enabled, these partial products must be evaluated upon construction of
    // FermionMatrix.
    PartialProducts = vector<cx_mat>();
    if ( params.ENABLE_PARTIAL_PRODUCTS ) reevaluate();
}

FermionMatrix::FermionMatrix( const FermionMatrix &obj ) {
    // Copy primitive data fields.
    params = obj.params;
    ptr_sigma = obj.ptr_sigma;

    // Copy vector of U matrices. Copy constructor of vector will call copy constructor of UMatrix.
    UProduct = vector<UMatrix>( obj.UProduct );

    // Copy vector of partial products.
    PartialProducts = vector<cx_mat>( obj.PartialProducts );
}

cx_mat FermionMatrix::getMatrix() {
    if ( params.ENABLE_PARTIAL_PRODUCTS ) {
        return eye( params.NX, params.NX ) + FullProduct;
    } else {
        cx_mat product( params.NX, params.NX );
        product.eye();

        for ( int i = 0; i < params.NTAU; i++ ) {
            product *= UProduct[ i ].getMatrix();
        }

        return eye( params.NX, params.NX ) + product;
    }
}

cx_mat FermionMatrix::getDerivative( int delta_x, int delta_tau ) {
    cx_mat product( params.NX, params.NX );
    product.eye();

    if ( params.ENABLE_PARTIAL_PRODUCTS ) {
        const int PARTIAL_PRODUCT_SIZE = params.NTAU / params.N_PARTIAL_PRODUCTS;
        const int DELTA_TAU_PP_INDEX = delta_tau / PARTIAL_PRODUCT_SIZE;

        // Multiply successive products of either PartialProducts[i] up to i = DELTA_TAU_PP_INDEX.
        for ( int i = 0; i < DELTA_TAU_PP_INDEX; i++ ) {
            product *= PartialProducts[ i ];
        }

        // Multiply successive products of UProduct[i] starting from where the last partial product ended, up to
        // the contribution UProduct[delta_tau - 1].
        for ( int i = DELTA_TAU_PP_INDEX * PARTIAL_PRODUCT_SIZE; i < delta_tau; i++ ) {
            product *= UProduct[ i ].getMatrix();
        }

        // Multiply by functional derivative of UProduct[delta_tau] at (delta_tau, delta_x).
        product *= UProduct[ delta_tau ].getDerivative( delta_x );

        // Multiply successive products of UProduct[i] starting from delta_tau + 1, up to where the next partial product
        // picks up.
        int upperBound = PARTIAL_PRODUCT_SIZE * ( DELTA_TAU_PP_INDEX + 1 );
        upperBound = ( upperBound < params.NTAU ) ? upperBound : params.NTAU;

        for ( int i = delta_tau + 1; i < upperBound; i++ ) {
            product *= UProduct[ i ].getMatrix();
        }

        // Multiply by remaining partial products.
        for ( int i = DELTA_TAU_PP_INDEX + 1; i < params.N_PARTIAL_PRODUCTS; i++ ) {
            product *= PartialProducts[ i ];
        }

        return product;
    } else {
        for ( int i = 0; i < params.NTAU; i++ ) {
            if ( i == delta_tau ) {
                product *= UProduct[ i ].getDerivative( delta_x );
            } else {
                product *= UProduct[ i ].getMatrix();
            }
        }

        return product;
    }
}

void FermionMatrix::reevaluate() {
    // Reevaluate U matrices for the updated sigma field.
    for ( int i = 0; i < params.NTAU; i++ ) {
        UProduct[ i ].reevaluate();
    }

    // Reevaluate partial products of the U matrices.
    if ( params.ENABLE_PARTIAL_PRODUCTS ) {
        const int partialProductsSize = params.NTAU / params.N_PARTIAL_PRODUCTS;
        PartialProducts.clear();

        FullProduct = cx_mat( params.NX, params.NX );
        FullProduct.eye();

        int tau = 0;
        while ( tau < params.NTAU ) {
            cx_mat partialProduct( params.NX, params.NX );
            partialProduct.eye();

            for ( int i = 0; i < partialProductsSize; i++ ) {
                partialProduct *= UProduct[ tau ].getMatrix();
                tau++;
                if ( tau == params.NTAU ) break;
            }

            FullProduct *= partialProduct;
            PartialProducts.push_back( partialProduct );
        }
    }
}
