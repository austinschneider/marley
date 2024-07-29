// Standard library includes
#include <fstream>
#include <iostream>
#include <vector>

// HepMC3 includes
#include "HepMC3/GenEvent.h"

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/EventFileReader.hh"
#include "marley/Generator.hh"
#include "marley/JSONConfig.hh"
#include "marley/Weighter.hh"

int main( int argc, char* argv[] ) {

  // If the user has not supplied enough command-line arguments, display the
  // standard help message and exit
  if (argc <= 3) {
    std::cout << "Usage: " << argv[0] << " OUTPUT_FILE CONFIG_FILE"
      << " INPUT_FILES\n";
    return 0;
  }

  // Check whether the output file exists and warn the user before
  // overwriting it if it does
  std::ifstream temp_stream( argv[1] );
  if ( temp_stream ) {
    bool overwrite = marley_utils::prompt_yes_no(
      "Really overwrite " + std::string(argv[1]) + '?' );
    if ( !overwrite ) {
      std::cout << "Action aborted.\n";
      return 0;
    }
  }

  // Get the configuration file name from the command line
  std::string config_file_name( argv[2] );

  // Set up the generator using the job configuration file. This will
  // also set up the owned Weighter object that will actually calculate
  // the weights.
  marley::JSONConfig jc( config_file_name );
  marley::Generator gen = jc.create_generator();

  // Don't include the (trivial) central-value weight in the output
  auto& weighter = gen.get_weighter();
  weighter.set_use_cv_weight( false );

  // Open the output file that will store the weight values
  std::ofstream outfile( argv[1] );
  if ( !outfile.good() ) {
    std::cerr << "Unable to open the output file \"" << argv[1]
      << "\" for writing\n";
    return 1;
  }

  // Write a header line in the output file describing the data format
  outfile << "event_num";
  for ( const auto& name : weighter.get_weight_names() ) {
    outfile << ' ' << name;
  }
  outfile << '\n';

  // Prepare to read the input file(s)
  std::vector< std::string > input_file_names;
  for ( int i = 3; i < argc; ++i ) input_file_names.push_back( argv[i] );

  // File loop
  for (const auto& file_name : input_file_names ) {

    // Open the current file for reading
    marley::EventFileReader efr( file_name );
    std::cout << "Opened file \"" << file_name << "\"\n";

    // Temporary object to use for reading in saved events
    HepMC3::GenEvent ev;

    int event_num = 0;

    // Event loop
    while ( efr >> ev ) {

      std::cout << "Event " << event_num << '\n';

      auto weight_vec = weighter.compute_weights( ev, gen );

      // Store the new row of weights in the output file
      outfile << event_num;
      for ( const auto& w : weight_vec ) {
        outfile << ' ' << w;
      }
      outfile << '\n';

      // Advance to the next event
      ++event_num;

    } // event loop

  } // input file loop

  return 0;
}
