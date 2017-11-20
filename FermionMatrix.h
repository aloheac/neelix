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
    UMatrix( unsigned int NX, unsigned int tau, double dtau, double mu, double g, SigmaField* ptr_sigma );

    arma::cx_mat getMatrix();

    arma::cx_mat getDerivative( int delta_x );

    void reevaluate();

    std::string to_string();

private:
    void evaluateElements();

    void evaluateDerivative( int delta_x );

    unsigned int NX;

    unsigned int tau;

    double dtau;

    double mu;

    double g;

    int dU_delta_x;

    arma::cx_mat U;

    arma::cx_mat dU;

    arma::cx_mat T;

    SigmaField* ptr_sigma;
};

class FermionMatrix {
public:
    FermionMatrix( unsigned int NX, unsigned int NTAU, double g, double dtau, double mu, SigmaField* sigma );

    arma::cx_mat getMatrix();

    arma::cx_mat getDerivative( int delta_x, int delta_tau );

    void reevaluate();

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
