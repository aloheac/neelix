#include <iostream>
#include "AuxiliaryField.h"

using namespace std;

int main() {
    AuxiliaryField sigma( 2, 10, 10 );
    sigma.set( 5, 3, complex<double>( 1., 3.));
    cout << sigma.get( 5, 3 ) << endl;

    cout << sigma.to_string() << endl;
    return 0;
}