//
// Created by loheac on 12/11/17.
//

#ifndef NEELIX_FIELDRECORDER_H
#define NEELIX_FIELDRECORDER_H

#include "AuxiliaryField.h"

class FieldRecorder {
public:
    FieldRecorder( MCParameters params, AuxiliaryField* sigma );

    void write( int n );

private:
    MCParameters params;

    AuxiliaryField* sigma;
};

#endif //NEELIX_FIELDRECORDER_H