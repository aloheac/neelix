//
// Created by loheac on 12/11/17.
//

#ifndef NEELIX_FIELDRECORDER_H
#define NEELIX_FIELDRECORDER_H

#include "H5Cpp.h"
#include "AuxiliaryField.h"

class FieldRecorder {
public:
    FieldRecorder( MCParameters params, AuxiliaryField* sigma );

    void write( int n );

private:
    MCParameters params;

    AuxiliaryField* sigma;
};


class FieldSetRecorder {
public:
    FieldSetRecorder( const char *filename, const MCParameters &params );

    void write( AuxiliaryField* sigma, std::complex<double> action );

private:
    const char* filename;

    MCParameters params;
};
#endif //NEELIX_FIELDRECORDER_H
