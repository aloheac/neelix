//
// Created by loheac on 8/31/17.
//

#include <iostream>
#include <sstream>
#include <iomanip>
#include "AuxiliaryField.h"
#include "FermionMatrix.h"


using namespace std;

class UnitTest {
public:
    UnitTest( string name, string (*tst)(), string result ) {
        cout << "Running unit test '" << name << "'..." << endl;
        string ret = tst();

        if ( ret == result ) {
            passedTests++;
            cout << ">> Unit test PASSED." << endl;
        } else {
            failedTests++;
            cout << ">> Unit test FAILED. Result:" << endl;
            cout << ret << endl;
            cout << "   Correct result: " << endl;
            cout << result << endl;
        }
    }

    static int passedTests;

    static int failedTests;
};

int UnitTest::passedTests = 0;
int UnitTest::failedTests = 0;

string A01() {
    stringstream ss;
    AuxiliaryField A( 1, 5, 5 );
    ss << A.to_string();
    return ss.str();
}

string A02() {
    stringstream ss;
    AuxiliaryField A( 1, 5, 5 );
    A.initialize();
    ss << A.to_string();
    return ss.str();
}

int main( int argc, char** argv ) {
    cout << "**********************************************************************" << endl;
    cout << "  Neelix Complex Langevin Monte Carlo Unit Testing" << endl;
    cout << "**********************************************************************" << endl << endl;

    UnitTest( "A01: AuxiliaryField, Constructor, to_string()", &A01, "(0,0)    (0,0)    (0,0)    (0,0)    (0,0)    \n(0,0)    (0,0)    (0,0)    (0,0)    (0,0)    \n(0,0)    (0,0)    (0,0)    (0,0)    (0,0)    \n(0,0)    (0,0)    (0,0)    (0,0)    (0,0)    \n(0,0)    (0,0)    (0,0)    (0,0)    (0,0)    " );

    UnitTest( "A02: AuxiliaryField, initialize()", &A02, "" );

    cout << "----------------------------------------------------------------------" << endl;
    cout << UnitTest::passedTests << " tests PASSED, " << UnitTest::failedTests << " tests FAILED." << endl;
}
