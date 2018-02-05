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
    UMatrix( int tau, MCParameters params, SigmaField* ptr_sigma );

    UMatrix( const UMatrix &obj );

    arma::cx_mat getMatrix();

    arma::cx_mat getDerivative( int delta_x );

    void reevaluate();

    std::string to_string();

private:
    void evaluateElements();

    void evaluateDerivative( int delta_x );

    MCParameters params;

    int tau;

    int dU_delta_x;

    std::complex<double> checksum;

    arma::cx_mat U;

    arma::cx_mat dU;

    arma::cx_mat T;

    arma::cx_mat T_x;

    SigmaField* ptr_sigma;
};

class FermionMatrix {
public:
    FermionMatrix( MCParameters params, SigmaField* sigma );

    FermionMatrix( const FermionMatrix &obj );

    arma::cx_mat getMatrix();

    arma::cx_mat getDerivative( int delta_x, int delta_tau );

    void reevaluate();

    SigmaField* ptr_sigma;

private:
    MCParameters params;

    std::vector<UMatrix> UProduct;

    std::vector<arma::cx_mat> PartialProducts;

    arma::cx_mat FullProduct;

};
#endif //NEELIX_FERMIONMATRIX_H
