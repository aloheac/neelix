#include <iostream>
#include "AuxiliaryField.h"
#include "FermionMatrix.h"
#include "CL.h"

using namespace std;

int main() {
    int NX = 5;
    int NTAU = 5;
    double dt = 0.2;
    double dtau = 0.05;

    SigmaField sigma( 1, 5, 5 );
    sigma.initialize();
    //sigma.set( 2, 2, complex<double>( 1., 3.));
    //cout << sigma.get( 1, 1 ) << endl;

    //cout << sigma.to_string() << endl;

    UMatrix U( 5, 1, &sigma );
    U.evaluateElements( 0.1, 0.05, 0.0 );
    cout << U.to_string() << endl;

    FermionMatrix M( 5, 5, 1, 0.05, 0, &sigma );
    //cout << M.evaluate_derivative(1,1) << endl;

    CLEvolver cl( NX, NTAU, dt, &sigma );
    cout << cl.calculateSigmaDot() << endl;
    return 0;
}
