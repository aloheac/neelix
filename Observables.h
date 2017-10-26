//
// Created by loheac on 10/24/17.
//

#ifndef NEELIX_OBSERVABLES_H
#define NEELIX_OBSERVABLES_H

#include "CL.h"
#include "Observables.h"
#include "FermionMatrix.h"

double mean( std::vector<double> vec );

std::complex<double> calculate_density( MCParameters params, SigmaField* sigma );

#endif //NEELIX_OBSERVABLES_H
