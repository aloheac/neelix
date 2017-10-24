#include <iostream>
#include "AuxiliaryField.h"
#include "FermionMatrix.h"
#include "CL.h"

using namespace std;

int main() {
    MCParameters params;
    params.NX = 5;
    params.NTAU = 5;
    params.dt = 0.2;
    params.g = 0.35;
    params.mu = 0.0;
    params.dtau = 0.05;

    SigmaField sigma( 1, params.NX, params.NTAU );
    sigma.initialize();

    CLEvolver cl( params, &sigma );
    FermionMatrix M( params.NX, params.NTAU, params.g, params.dtau, params.mu, &sigma );
    cout << log( pow( det( M.evaluate() ), 2 ) ) << endl;
    cout << sigma.to_string() << endl;
    cl.integrateSigma();
    cout << log( pow( det( M.evaluate() ), 2 ) ) << endl;
    cout << sigma.to_string() << endl;
    return 0;
}
