//
// Created by loheac on 8/22/17.
//

#ifndef NEELIX_AUXILIARYFIELD_H
#define NEELIX_AUXILIARYFIELD_H

#include <complex>
#include <exception>
#include <string>
#include <random>

class AuxiliaryFieldException : public std::exception {
public:
    AuxiliaryFieldException( std::string msg );

    virtual const char* what() const throw();

private:
    std::string msg;
};

class AuxiliaryField {
public:

    AuxiliaryField( int thisDimension, int thisNx, int thisNtau );

    ~AuxiliaryField();

    void set( unsigned int x, unsigned int tau, std::complex<double> value );

    std::complex<double> get( unsigned int x, unsigned int tau );

    std::string to_string();

    std::complex<double> sum();

    virtual void initialize();

    const int dimension;

    const int NX;

    const int NTAU;

protected:
    std::default_random_engine rand_generator;

    std::complex<double>** elements;

};

class SigmaField : public AuxiliaryField {
public:
    SigmaField( int thisDimension, int thisNx, int thisNtau );

    void initialize();

};

class MomentumField : public AuxiliaryField {
public:
    MomentumField( int thisDimension, int thisNx, int thisNtau );

    void initialize();

};

#endif //NEELIX_AUXILIARYFIELD_H
