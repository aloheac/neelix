#include <iostream>
#include <vector>
#include "AuxiliaryField.h"
#include "FermionMatrix.h"
#include "CL.h"
#include "HMC.h"
#include "Observables.h"
#include <complex.h>
using namespace std;

int main() {
    MCParameters params;
    params.NX = 21;
    params.NTAU = 60;
    params.dt = 0.01;
    params.g = 0.5;
    params.mu = 0.0;
    params.dtau = 0.1;
    params.xi = 10;

    cout << "Input parameters:" << endl;
    cout << "    NX: " << params.NX << endl;
    cout << "    NTAU: " << params.NTAU << endl;
    cout << "    MU: " << params.mu << endl;
    cout << "    BARE_COUPLING: " << params.g << endl;
    cout << "    DT: " << params.dt << endl;
    cout << "    DTAU: " << params.dtau << endl;
    cout << "    XI: " << params.xi << endl;
    cout << endl;

    cout << "Calculated parameters:" << endl;
    cout << "    BETA: " << params.NTAU * params.dtau << endl;
    cout << "    LAMBDA: " << sqrt( params.NTAU * params.dtau ) * params.g << endl;
    cout << "    BETA MU: " << params.NTAU * params.dtau * params.mu << endl;
    cout << endl;

    SigmaField sigma( 1, params.NX, params.NTAU );
    sigma.initialize();

    MomentumField pi( 1, params.NX, params.NTAU );
    pi.initialize();

    HMCEvolver hmc( params, &sigma, &pi );
    vector<double> density;

    int num_samples = 0;
    for ( int i = 0; i < 20000; i++ ) {
        complex<double> next_density = calculate_density( params, &sigma );
        if ( num_samples == 20 ) {
            num_samples = 0;
        density.push_back( next_density.real() );
        }
        num_samples++;

        cout << i << "    " << next_density << "    " << mean( density ) << "    " << hmc.params.dt << endl;
        hmc.integrateSigma();
    }

    cout << sigma.to_string() << endl;
    return 0;
}
