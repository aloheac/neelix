//
// Created by loheac on 8/28/17.
//

#include "FermionMatrix.h"

// Neelix currently uses the Armadillo C++ Linear Algebra Library to implement matrix operations.
using namespace std;
using namespace arma;

UMatrix::UMatrix( unsigned int thisNX, unsigned int thisTau, SigmaField* thisSigma ) : NX( thisNX ), tau( thisTau ), ptr_sigma( thisSigma ) {
    U = cx_mat( NX, NX );  // cx_mat is an Armadillo typedef for Mat< std::complex<double> > (complex dense matrix).
    U.zeros();
    basis = MatrixBasis::COORDINATE;
}

void UMatrix::switchBasisFlag() {
    if ( basis ==  MatrixBasis::COORDINATE ) {
        basis = MatrixBasis::MOMENTUM;
    } else {
        basis = MatrixBasis::COORDINATE;
    }
}

void UMatrix::evaluateElements( double g, double dtau, double mu ) {
    cx_mat S( NX, NX );
    complex<double> A = sqrt( 2.0 * ( exp( dtau * g ) - 1.0 ) );

    for ( unsigned int i = 0; i < NX; i++ ) {
            S(i, i) = A * ( 1.0 + sin( ptr_sigma->get( i, tau ) ) );
    }

    U = S;
}

string UMatrix::to_string() {
    stringstream ss;
    ss << "UMatrix[" << endl;
    ss << U;
    ss << "]";

    return ss.str();
}
