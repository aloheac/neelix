//
// Created by loheac on 9/3/17.
//

#ifndef NEELIX_HMC_H
#define NEELIX_HMC_H

#include "AuxiliaryField.h"
#include "CL.h"

class HMCEvolver {
public:
    HMCEvolver( MCParameters params, SigmaField* sigma, MomentumField* pi );

    void integrateSigma();

    MCParameters params;

private:

    arma::cx_mat calculatePiDot();

    SigmaField* sigma;

    MomentumField* pi;
};

#endif //NEELIX_HMC_H
