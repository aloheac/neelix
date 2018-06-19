//
// Created by loheac on 12/11/17.
//

#include <fstream>
#include <sstream>
#include "FieldRecorder.h"

using namespace std;

FieldRecorder::FieldRecorder( MCParameters this_params, AuxiliaryField *this_sigma ) : params( this_params ), sigma( this_sigma ) {

}

void FieldRecorder::write( int n ) {
    stringstream ss;
    ss << "./" << "Field_" << n << ".dat";

    ofstream ofs( ss.str() );

    for ( int x = 0; x < params.NX; x++ ) {
        for ( int t = 0; t < params.NTAU; t++ ) {
            ofs << x << "\t" << t << "\t" << sigma->get( x, t ).real() << "\n";
        }
    }

    ofs.close();
}

using namespace H5;

FieldSetRecorder::FieldSetRecorder( const char *filename, const MCParameters &params ) : filename( filename ), params( params ) {
    // Initialize HDF5 data file with field configuration count.
    H5File hdf( H5std_string( filename ), H5F_ACC_TRUNC );

    hsize_t integer[] = { 1 };
    DataSpace dataspace( 1, integer );
    DataSet dataset = hdf.createDataSet( "/Count", PredType::NATIVE_INT, dataspace );
    dataset.write( new int( 0 ), PredType::NATIVE_INT );

    // Create main groups.
    hdf.createGroup( "/FieldConfiguration" );
    hdf.createGroup( "/Action" );

    hdf.close();
}

void FieldSetRecorder::write( AuxiliaryField* sigma, std::complex<double> action ) {
    cout << "Writing field configurations to HDF5 data set..." << endl;

    /*
     *  INITIALIZE
     */

    // Open HDF5 file for writing.
    H5File hdf( H5std_string( filename ), H5F_ACC_RDWR );

    // Establish field dimensions.
    hsize_t field_dims[] = { params.NX, params.NTAU };
    hsize_t single_dims[] = { 1 };

    // Define data space for field configuration.
    DataSpace dataspace_field( 2, field_dims );  // 2-dimensional array.
    DataSpace dataspace_single( 1, single_dims );

    // Find last field configuration index.
    DataSet dataset_index = hdf.openDataSet( "/Count" );

    int *prev_index = new int();
    dataset_index.read( prev_index, PredType::NATIVE_INT, dataspace_single );
    int next_index = *prev_index + 1;

    // Update field configuration count.
    dataset_index.write( &next_index, PredType::NATIVE_INT );

    /*
     *  WRITE FIELD CONFIGURATION
     */

    // Generate real and imaginary components of fields.
    double sigma_real[ params.NX ][ params.NTAU ];
    double sigma_imag[ params.NX ][ params.NTAU ];

    for ( int x = 0; x < params.NX; x++ ) {
        for ( int y = 0; y < params.NTAU; y++ ) {
            sigma_real[ x ][ y ] = sigma->get( x, y ).real();
            sigma_imag[ x ][ y ] = sigma->get( x, y ).imag();
        }
    }

    // Generate required groups.
    stringstream ss_group, ss_re, ss_im;
    ss_group << "/FieldConfiguration/" << next_index;
    ss_re << "/FieldConfiguration/" << next_index << "/Re";  // Real component of field.
    ss_im << "/FieldConfiguration/" << next_index << "/Im"; // Imaginary component of field.

    hdf.createGroup( ss_group.str().c_str() );

    // Create dataset and write field configuration.
    DataSet dataset = hdf.createDataSet( ss_re.str().c_str(), PredType::NATIVE_DOUBLE, dataspace_field );
    dataset.write( sigma_real,  PredType::NATIVE_DOUBLE );

    dataset = hdf.createDataSet( ss_im.str().c_str(),  PredType::NATIVE_DOUBLE, dataspace_field );
    dataset.write( sigma_imag, PredType::NATIVE_DOUBLE );

    /*
     *  WRITE ACTION
     */

    double action_re = action.real();
    double action_im = action.imag();

    // Generate required groups.
    ss_group.str("");    ss_re.str("");    ss_im.str("");  // Reset stringstreams.
    ss_group << "/Action/" << next_index;
    ss_re << "/Action/" << next_index << "/Re";  // Real component of action.
    ss_im << "/Action/" << next_index << "/Im"; // Imaginary component of action.

    hdf.createGroup( ss_group.str() );

    // Create dataset and write value of action.
    dataset = hdf.createDataSet( ss_re.str(), PredType::NATIVE_DOUBLE, dataspace_single );
    dataset.write( &action_re, PredType::NATIVE_DOUBLE );

    dataset = hdf.createDataSet( ss_im.str(), PredType::NATIVE_DOUBLE, dataspace_single );
    dataset.write( &action_im, PredType::NATIVE_DOUBLE );

    hdf.close();
}
