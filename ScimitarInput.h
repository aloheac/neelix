/*
 * C++ Library for Parsing Scimitar-style Parameter Input Files.
 *
 * This library accepts the filename of an input file containing a series of
 * parameters in the Scimitar input format and provides a structure containing
 * field names and values.
 *
 * Andrew C. Loheac, Joaquin E. Drut
 * Department of Physics and Astronomy
 * University of North Carolina at Chapel Hill
 *
 * 18 January 2014 (original library)
 * rev. 26 January 2014, 1 October 2016
 */

#ifndef APOPHIS_SCIMITARINPUT_H
#define APOPHIS_SCIMITARINPUT_H

#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>

using namespace std;

class ScimitarInputParser {
private:
    // Character that delimits field value from the field name in a Scimitar
    // input file. Default value is '#'.
    char delimiter;

    // Input file stream for the Scimitar input file.
    ifstream ifs;

    // Boolean flag to indicate whether or not an input file was successfully
    // opened.
    bool fileIsOpen;

    // The number of valid fields read in from the input file.
    int numFields;

    // Container for holding the input values.
    vector<string> inputValues;

    // Container for holding the input field names.
    vector<string> inputFieldNames;

    /**
     * Trims tabs and spaces from the right hand side of an input string.
     * Method mutates the input and returns the modified string.
     *
     * @param String to be trimmed.
     * @return The trimmed string.
     */
    string rtrim( string str );

    /**
     * Trims tabs and spaces frmo the left hand side of an input string.
     * Method mutates the input and returns the modified string.
     *
     * @param String to be trimmed.
     * @return The trimmed string.
     */
    string ltrim( string str );

    /**
     * Trims tabs and spaces from the left hand and right hand side of an
     * input string. Method mutates the input and returns the modified string.
     *
     * @param String to be trimmed.
     * @return The trimmed string.
     */
    string trim( string str );

public:
    /**
     * Default constructor. Takes no arguments.
     */
    ScimitarInputParser() : delimiter( '#' ) , fileIsOpen( false ),
                            numFields( 0 ) { }

    /**
     * Constructor. Takes an input filename and opens the file.
     *
     * @param filename Scimitar input filename.
     */
    ScimitarInputParser( const char* filename ) : delimiter( '#' ),
                                                  fileIsOpen( false ),
                                                  numFields( 0 ) {
        openInputFile( filename );
    }

    /**
     * Opens the file stream for an input file.
     *
     * @param filename Input filename.
     * @return 1 if the operation is successful, -1 otherwise.
     */
    int openInputFile( const char* filename );

    /**
     * Parses the input file and enters the field values and names into their
     * respective containers.
     *
     * @return 1 if the operation is successful, -1 otherwise.
     */
    int parse();

    /**
     * Returns the value of a given field as a string.
     *
     * @param field Index of the field to be accessed.
     * @return The value of the field as the original string.
     */
    string getValue( int field );

    /**
     * Returns the name of a field.
     *
     * @param field Index of the field to be accessed.
     */
    string getName( int field );

    /**
     * Returns the value of a given field casted to an integer.
     *
     * @param field Index of the field to be accessed.
     * @return The value of the field as an integer.
     */
    int getValueInt( int field );

    /**
     * Returns the value of a given field casted to a double.
     *
     * @param field Index of the field to be accessed.
     * @return The value of the field as a double.
     */
    double getValueDouble( int field );
};

// ***** METHOD IMPLEMENTATIONS *****

string ScimitarInputParser::rtrim( string str ) {
    while ( str.length() > 0 && ( str.at( str.length() - 1 )  == ' ' || str.at( str.length() - 1 ) == '\t' ) ) {
        str.erase( str.length() - 1 );
    }

    return str;
}

string ScimitarInputParser::ltrim( string str ) {
    while ( str.length() > 0 && ( str.at( 0 ) == ' ' || str.at( 0 ) == '\t' ) ) {
        str.erase( 0 );
    }

    return str;
}

string ScimitarInputParser::trim( string str ) {
    return rtrim( ltrim( str ) );
}

int ScimitarInputParser::openInputFile( const char* filename ) {
    ifs.open( filename );

    if ( !ifs.good() ) {
        fileIsOpen = false;
        printf( "***ERROR (Scimitar Input Parser) : Could not open input file!\n" );
        return -1;
    } else {
        fileIsOpen = true;
        return 1;
    }
}

int ScimitarInputParser::parse() {
    string line;
    int i_field = 0;
    int delimiterPosition;

    if ( fileIsOpen && ifs.is_open() ) {
        while( getline( ifs, line ) ) {
            delimiterPosition = line.find( delimiter );
            inputValues.push_back( trim( line.substr( 0, delimiterPosition ) ) );
            inputFieldNames.push_back( trim( line.substr( delimiterPosition + 1, line.length() ) ) );
            numFields++;
        }
        ifs.close();
        return 1;
    } else {
        printf( "***ERROR (Scimitar Input Parser) : Could not open input file!\n" );
        return -1;
    }
}

string ScimitarInputParser::getValue( int field ) {
    if ( field > numFields ) {
        printf( "***ERROR (ScimitarInputParser) : Field index is out of range.\n" );
        return string();
    } else {
        return inputValues.at( field );
    }
}

string ScimitarInputParser::getName( int field ) {
    if ( field > numFields ) {
        printf( "***ERROR (ScimitarInputParser) : Field index is out of range.\n" );
        return string();
    } else {
        return inputFieldNames.at( field );
    }
}

int ScimitarInputParser::getValueInt( int field ) {
    if ( field > numFields ) {
        printf( "***ERROR (ScimitarInputParser) : Field index is out of range.\n" );
        return 0;
    } else {
        return atoi( getValue( field ).c_str() );
    }
}

double ScimitarInputParser::getValueDouble( int field ) {
    if ( field > numFields ) {
        printf( "***ERROR (ScimitarInputParser) : Field index is out of range.\n" );
        return 0;
    } else {
        return atof( getValue( field ).c_str() );
    }
}

#endif //APOPHIS_SCIMITARINPUT_H
