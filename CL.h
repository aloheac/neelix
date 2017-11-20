//
// Created by loheac on 10/24/17.
//

#ifndef NEELIX_CL_H
#define NEELIX_CL_H

#include <random>
#include <armadillo>
#include "AuxiliaryField.h"


class CLEvolver {
public:
    CLEvolver( MCParameters this_params, SigmaField* sigma );

    arma::cx_mat calculateSigmaDot();

    void integrateSigma();

    MCParameters params;

private:

    SigmaField* sigma;

    SigmaField* sigma_prime;

    std::normal_distribution<double> rand_normal;

    std::default_random_engine rand_generator;

};

#endif //NEELIX_CL_H
