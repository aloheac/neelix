#include <iostream>
#include "AuxiliaryField.h"
#include "FermionMatrix.h"
#include "CL.h"

using namespace std;

int main() {
    MCParameters params;
    params.NX = 30;
    params.NTAU = 5;
    params.dt = 0.2;
    params.g = 0.5;
    params.mu = 0.5;
    params.dtau = 0.05;

    SigmaField sigma( 1, params.NX, params.NTAU );
    sigma.initialize();

    CLEvolver cl( params, &sigma );
    FermionMatrix M( params.NX, params.NTAU, params.g, params.dtau, params.mu, &sigma );

    for ( int i = 0; i < 100; i++ ) {
        cout << log( pow( det( M.evaluate() ), 2 ) ) << endl;
        cl.integrateSigma();
    }

    return 0;
}
