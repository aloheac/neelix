//
// Created by loheac on 8/22/17.
//

#include <iostream>
#include <sstream>
#include <cstdlib>
#include "AuxiliaryField.h"

using namespace std;

std::default_random_engine AuxiliaryField::rand_generator = std::default_random_engine( 34563546 );

AuxiliaryFieldException::AuxiliaryFieldException( string thisMsg ) : exception(), msg( thisMsg ) { };

const char* AuxiliaryFieldException::what() const throw() {
    stringstream ss;
    ss << "AuxiliaryFieldException: " << msg;
    cout << "***ERROR: " << ss.str() << endl;
    return ss.str().c_str();
}

AuxiliaryField::AuxiliaryField( int thisDimension, int thisNx, int thisNtau ) : dimension( thisDimension ), NX( thisNx ), NTAU( thisNtau ) {
    if ( dimension == 1 ) {
        elements = new complex<double>*[ NX ];
        for ( int x = 0; x < NX; x++ ) {
            elements[ x ] = new complex<double>[ NTAU ];

            for ( int tau = 0; tau < NTAU; tau++ ) {  // Initialize values of array.
                elements[ x ][ tau ] = complex<double>( 0.0, 0.0 );
            }
        }
    } else {
        throw AuxiliaryFieldException( "AuxiliaryField: Requested field of spatial dimension other than 1; higher dimensions not yet supported." );
    }
}

AuxiliaryField::~AuxiliaryField() {
    for ( int i = 0; i < NX; i++ ) {
        delete[] elements[ i ];
    }

    delete[] elements;
}

void AuxiliaryField::set( unsigned int x, unsigned int tau, std::complex<double> value ) {
    if ( x >= NX ) throw AuxiliaryFieldException( "AuxiliaryField: Requested to set spatial element greater than NX dimension." );
    if ( tau >= NTAU ) throw AuxiliaryFieldException( "AuxiliaryField: Requested to set temporal element greater than NTAU dimension." );

    elements[ x ][ tau ] = value;
}

std::complex<double> AuxiliaryField::get( unsigned int x, unsigned int tau ) {
    if ( x >= NX ) throw AuxiliaryFieldException( "AuxiliaryField: Requested to get spatial element greater than NX dimension." );
    if ( tau >= NTAU ) throw AuxiliaryFieldException( "AuxiliaryField: Requested to get temporal element greater than NTAU dimension." );

    return elements[ x ][ tau ];
}

string AuxiliaryField::to_string() {
    stringstream ss;

    for ( int x = 0; x < NX; x++ ) {
        for ( int tau = 0; tau < NTAU; tau++ ) {
            ss << elements[ x ][ tau ] << "    ";
        }
        ss << endl;
    }
    ss << endl;

    return ss.str();
}

complex<double> AuxiliaryField::sum() {
    complex<double> s;

    for ( int x = 0; x < NX; x++ ) {
        for ( int tau = 0; tau < NTAU; tau++ ) {
            s += elements[ x ][ tau ];
        }
    }

    return s;
}

void AuxiliaryField::initialize() {  }

SigmaField::SigmaField( int thisDimension, int thisNx, int thisNtau ) : AuxiliaryField( thisDimension, thisNx, thisNtau ) { }

void SigmaField::initialize() {
    uniform_real_distribution<double> dist( -1.0, 1.0 );
    double re_rand, im_rand;

    for ( int x = 0; x < NX; x++ ) {
        for ( int tau = 0; tau < NTAU; tau++ ) {
            re_rand = dist( rand_generator );
            im_rand = 0.0;
            elements[ x ][ tau ]  = complex<double>( re_rand, im_rand );
        }
    }
}

MomentumField::MomentumField( int thisDimension, int thisNx, int thisNtau ) : AuxiliaryField( thisDimension, thisNx, thisNtau ) { }

void MomentumField::initialize() {
    normal_distribution<double> rand_normal( 0.0, 1.0 );
    double re_rand, im_rand;

    for ( int x = 0; x < NX; x++ ) {
        for ( int tau = 0; tau < NTAU; tau++ ) {
            re_rand = rand_normal( rand_generator );
            im_rand = 0.0;
            elements[ x ][ tau ]  = complex<double>( re_rand, im_rand );
        }
    }
}
