#include <iostream>
#include "AuxiliaryField.h"
#include "FermionMatrix.h"

using namespace std;

int main() {
    SigmaField sigma( 1, 10, 10 );
    sigma.initialize();
    sigma.set( 5, 3, complex<double>( 1., 3.));
    cout << sigma.get( 5, 3 ) << endl;

    cout << sigma.to_string() << endl;

    UMatrix U( 10, 1, &sigma );
    U.evaluateElements( 0.1, 0.05, 0.0 );
    cout << U.to_string() << endl;
    return 0;
}
