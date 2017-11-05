#include <iostream>
#include <vector>
#include "AuxiliaryField.h"
#include "FermionMatrix.h"
#include "CL.h"
#include "HMC.h"
#include "Observables.h"

using namespace std;

int main() {
    MCParameters params;
    params.NX = 21;
    params.NTAU = 10;
    params.dt = 0.001;
    params.g = 0.7071;
    params.mu = -5;
    params.dtau = 0.2;
    params.xi = 0.01;

    cout << "Input parameters:" << endl;
    cout << "    NX: " << params.NX << endl;
    cout << "    NTAU: " << params.NTAU << endl;
    cout << "    MU: " << params.mu << endl;
    cout << "    BARE_COUPLING: " << params.g << endl;
    cout << "    DT: " << params.dt << endl;
    cout << "    DTAU: " << params.dtau << endl;
    cout << "    XI: " << params.xi << endl;

    SigmaField sigma( 1, params.NX, params.NTAU );
    sigma.initialize();

    CLEvolver cl( params, &sigma );
    vector<double> density;

    for ( int i = 0; i < 100; i++ ) {
        cl.integrateSigma();
        complex<double> next_density = calculate_density( params, &sigma );
        density.push_back( next_density.real() );

        cout << i << "    " << next_density << "    " << mean( density ) << "    " << cl.params.dt << endl;
    }

    cout << sigma.to_string() << endl;
    return 0;
}
