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

    void evaluateElements();

    void evaluateElementsOfDerivative( int delta_x );

    arma::cx_mat getMatrix();

    std::string to_string();

private:
    unsigned int NX;

    unsigned int tau;

    double dtau;

    double mu;

    double g;

    arma::cx_mat U;

    arma::cx_mat T;

    SigmaField* ptr_sigma;
};

class FermionMatrix {
public:
    FermionMatrix( unsigned int NX, unsigned int NTAU, double g, double dtau, double mu, SigmaField* sigma );

    arma::cx_mat evaluate();

    arma::cx_mat evaluate_derivative( int delta_x, int delta_tau );

    SigmaField* ptr_sigma;



    // Methods to be made private after debugging:
    arma::cx_mat evaluateUProduct();

    arma::cx_mat evaluateUProductDerivative( int delta_x, int delta_tau );

private:
    unsigned int NX;

    unsigned int NTAU;

    double g;

    double dtau;

    double mu;

    std::vector<UMatrix> UProduct;

};
#endif //NEELIX_FERMIONMATRIX_H
