//
// Created by loheac on 8/22/17.
//

#include <iostream>
#include <sstream>
#include "AuxiliaryField.h"

using namespace std;

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
        throw AuxiliaryFieldException( "Field is of spatial dimension other than 1; higher dimensions not yet supported." );
    }
}

AuxiliaryField::~AuxiliaryField() {
    for ( int i = 0; i < NX; i++ ) {
        delete[] elements[ i ];
    }

    delete[] elements;
}

void AuxiliaryField::set( unsigned int x, unsigned int tau, std::complex<double> value ) {
    elements[ x ][ tau ] = value;
}

std::complex<double> AuxiliaryField::get( unsigned int x, unsigned int tau ) {
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