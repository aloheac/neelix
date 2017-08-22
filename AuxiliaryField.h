//
// Created by loheac on 8/22/17.
//

#ifndef NEELIX_AUXILIARYFIELD_H
#define NEELIX_AUXILIARYFIELD_H

#include <complex>
#include <exception>
#include <string>

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

private:
    const int dimension;

    const int NX;

    const int NTAU;

    std::complex<double>** elements;
};

class SigmaField : public AuxiliaryField {
public:
    SigmaField( int thisDimension, int thisNx, int thisNtau );

};

class MomentumField : public AuxiliaryField {
public:
    MomentumField( int thisDimension, int thisNx, int thisNtau );

};

#endif //NEELIX_AUXILIARYFIELD_H
