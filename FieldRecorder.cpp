//
// Created by loheac on 12/11/17.
//

#include <fstream>
#include <sstream>
#include "FieldRecorder.h"

using namespace std;

FieldRecorder::FieldRecorder( MCParameters this_params, SigmaField *this_sigma ) : params( this_params ), sigma( this_sigma ) {

}

void FieldRecorder::write( int n ) {
    stringstream ss;
    ss << "./" << "Field_" << n << ".dat";

    ofstream ofs( ss.str() );

    for ( int x = 0; x < params.NX; x++ ) {
        for ( int t = 0; t < params.NTAU; t++ ) {
            ofs << sigma->get( x, t ).real() << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}
