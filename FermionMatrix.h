//
// Created by loheac on 8/28/17.
//

#ifndef NEELIX_FERMIONMATRIX_H
#define NEELIX_FERMIONMATRIX_H

#include <armadillo>
#include "AuxiliaryField.h"

enum MatrixBasis : int {
    COORDINATE = 0,
    MOMENTUM = 1
};

class UMatrix {
public:
    UMatrix( unsigned int NX, unsigned int tau, SigmaField* ptr_sigma );

    void switchBasisFlag();

    void evaluateElements( double g, double dtau, double mu );

    std::string to_string();

private:
    unsigned int NX;

    unsigned int tau;

    arma::cx_mat U;

    SigmaField* ptr_sigma;

    MatrixBasis basis;
};
#endif //NEELIX_FERMIONMATRIX_H
