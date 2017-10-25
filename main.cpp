#include <iostream>
#include "AuxiliaryField.h"
#include "FermionMatrix.h"
#include "CL.h"
#include "Observables.h"

using namespace std;

int main() {
    MCParameters params;
    params.NX = 10;
    params.NTAU = 20;
    params.dt = 0.2;
    params.g = 0.2;
    params.mu = 1.0;
    params.dtau = 0.05;

    SigmaField sigma( 1, params.NX, params.NTAU );
    sigma.initialize();

    CLEvolver cl( params, &sigma );

    for ( int i = 0; i < 100; i++ ) {
        cl.integrateSigma();
        cout << calculate_density( params, &sigma ) << endl;
    }


    return 0;
}
