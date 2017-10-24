//
// Created by loheac on 10/24/17.
//

#ifndef NEELIX_CL_H
#define NEELIX_CL_H

#include <armadillo>
#include "AuxiliaryField.h"

class CLEvolver {
public:
    CLEvolver( int NX, int NTAU, double dt, SigmaField* sigma );

    arma::cx_mat calculateSigmaDot();

    void integrateSigma();

private:

    double dt;

    int NX;

    int NTAU;

    SigmaField* sigma;

    SigmaField* sigma_prime;

};

#endif //NEELIX_CL_H
