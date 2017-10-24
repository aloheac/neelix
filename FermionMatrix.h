//
// Created by loheac on 8/28/17.
//

#ifndef NEELIX_FERMIONMATRIX_H
#define NEELIX_FERMIONMATRIX_H

#include <armadillo>
#include <vector>
#include "AuxiliaryField.h"

class UMatrix {
public:
    UMatrix( unsigned int NX, unsigned int tau, SigmaField* ptr_sigma );

    void evaluateElements( double g, double dtau, double mu );

    arma::cx_mat getMatrix();

    std::string to_string();

private:
    unsigned int NX;

    unsigned int tau;

    arma::cx_mat U;

    SigmaField* ptr_sigma;
};

class FermionMatrix {
public:
    FermionMatrix( unsigned int NX, unsigned int NTAU, double g, double dtau, double mu );

    ~FermionMatrix();

    arma::cx_mat evaluateUProduct();

    std::complex<double> evaluateLogDet();

    SigmaField* ptr_sigma;

private:
    unsigned int NX;

    unsigned int NTAU;

    double g;

    double dtau;

    double mu;

    std::vector<UMatrix> UProduct;

};
#endif //NEELIX_FERMIONMATRIX_H
