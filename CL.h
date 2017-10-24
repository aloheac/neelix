//
// Created by loheac on 10/24/17.
//

#ifndef NEELIX_CL_H
#define NEELIX_CL_H

#include <armadillo>
#include "AuxiliaryField.h"

struct MCParameters {
    int NX;

    int NTAU;

    double g;

    double beta;

    double mu;

    double dtau;

    double dt;
};

class CLEvolver {
public:
    CLEvolver( MCParameters this_params, SigmaField* sigma );

    arma::cx_mat calculateSigmaDot();

    void integrateSigma();

private:

    MCParameters params;

    SigmaField* sigma;

    SigmaField* sigma_prime;

};

#endif //NEELIX_CL_H
