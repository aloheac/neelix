#include <iostream>
#include <vector>
#include "AuxiliaryField.h"
#include "FermionMatrix.h"
#include "CL.h"
#include "HMC.h"
#include "Observables.h"
#include <complex.h>
#include "ScimitarInput.h"

using namespace std;

int main() {
    ScimitarInputParser sip( "params.inp" );
    sip.parse();

    MCParameters params;
    params.NX = sip.getValueInt( 0 );
    params.NTAU = sip.getValueInt( 1 );
    params.dt = sip.getValueDouble( 2 );
    params.g = sip.getValueDouble( 3 );
    params.mu = sip.getValueDouble( 4 );
    params.dtau = sip.getValueDouble( 5 );
    params.xi = sip.getValueDouble( 6 );
    params.ENABLE_FIELD_CHECKSUM = false;
    params.N_PARTIAL_PRODUCTS = 10;
    params.ENABLE_PARTIAL_PRODUCTS = params.N_PARTIAL_PRODUCTS > 1;

    int N_STEPS = sip.getValueInt( 7 );
    int OBSERVE_FREQ = sip.getValueInt( 8 );

    cout << "Input parameters:" << endl;
    cout << "    NX: " << params.NX << endl;
    cout << "    NTAU: " << params.NTAU << endl;
    cout << "    MU: " << params.mu << endl;
    cout << "    BARE_COUPLING: " << params.g << endl;
    cout << "    DT: " << params.dt << endl;
    cout << "    DTAU: " << params.dtau << endl;
    cout << "    XI: " << params.xi << endl;
    cout << "    N_STEPS: " << N_STEPS << endl;
    cout << "    OBSERVE_FREQ: " << OBSERVE_FREQ << endl;
    cout << endl;

    cout << "Calculated parameters:" << endl;
    cout << "    BETA: " << params.NTAU * params.dtau << endl;
    cout << "    LAMBDA: " << sqrt( params.NTAU * params.dtau ) * params.g << endl;
    cout << "    A: " << sqrt( 2.0 * ( complex<double>( exp( params.dtau * params.g ) ) - 1.0 ) ) << endl;
    cout << "    BETA MU: " << params.NTAU * params.dtau * params.mu << endl;
    cout << endl;

    cout << "Performance settings: " << endl;
    cout << "    FIELD CHECKSUM ENABLED: " << params.ENABLE_FIELD_CHECKSUM << endl;
    cout << "    N_PARTIAL_PRODUCTS: " << params.N_PARTIAL_PRODUCTS << endl;
    cout << "    PARTIAL PRODUCTS ENABLED: " << params.ENABLE_PARTIAL_PRODUCTS << endl;
    cout << endl;

    SigmaField sigma( 1, params.NX, params.NTAU );
    sigma.initialize();

    MomentumField pi( 1, params.NX, params.NTAU );
    pi.initialize();

    HMCEvolver hmc( params, &sigma, &pi );
    vector<double> density;

    int num_samples = 0;
    for ( int i = 0; i < N_STEPS; i++ ) {
        complex<double> next_density = calculate_density( params, &sigma );
        if ( num_samples == OBSERVE_FREQ ) {
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
