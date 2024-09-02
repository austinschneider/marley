// Standard library includes
#include <fstream>
#include <iostream>
#include <vector>

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/EventFileReader.hh"
#include "marley/Generator.hh"
#include "marley/JSONConfig.hh"
#include "marley/WeightCalculator.hh"
#include "marley/Weighter.hh"

int main( int argc, char* argv[] ) {

  // If the user has not supplied enough command-line arguments, display the
  // standard help message and exit
  if ( argc != 3 ) {
    std::cout << "Usage: " << argv[0] << " CONFIG_FILE INPUT_FILE\n";
    return 0;
  }

  // Get the configuration file name from the command line
  std::string config_file_name( argv[1] );

  // Check for the presence of a "weights" key in the configuration file
  marley::JSON rw_config = marley::JSON::load_file( config_file_name );
  if ( !rw_config.has_key("weights") ) throw marley::Error( "Missing"
    " \"weights\" key in marreweight configuration file" );

  // Configure a Weighter object to do the actual weight calculation.
  const auto& json_weights = rw_config.at( "weights" );
  marley::Weighter weighter( json_weights );

  // We don't need to compute the CV weight since it should already be present
  // in the previously generated events
  weighter.set_use_cv_weight( false );

  // Open the input file and retrieve the first event so that we can get
  // access to the run information (and thus the previous weight names)
  std::string input_file_name( argv[2] );
  marley::EventFileReader efr( input_file_name );

  HepMC3::GenEvent ev;
  efr >> ev;

  auto run_info = ev.run_info();
  const std::vector< std::string > wgt_names = run_info->weight_names();

  // Prepend trivial weight calculators to the Weighter object for each of the
  // previous weights. This will preserve the old weights and allow the new
  // weights to be stored in the correct location in the output vector.
  auto& calc_vec = weighter.get_weight_calculators();
  size_t num_new_weights = calc_vec.size();

  // Loop over the previous weight names in reverse so that their original
  // order is preserved in the reweighting output
  for ( auto riter = wgt_names.crbegin();
    riter != wgt_names.crend(); ++riter )
  {
    const auto& w_name = *riter;
    marley::JSON temp_json;
    temp_json[ "name" ] = w_name;
    auto w_calc = std::make_shared< marley
      ::TrivialWeightCalculator >( temp_json );
    calc_vec.insert( calc_vec.begin(), w_calc );
  }

  // Make an updated list of weight names that includes the new ones
  auto full_name_vec = weighter.get_weight_names();

  // Extract the prior generator JSON configuration from the run information
  auto prior_config_str = run_info->attribute< HepMC3::StringAttribute >(
    "MARLEY.JSONconfig" );

  if ( !prior_config_str ) {
    throw marley::Error( "Failed to retrieve previous generator"
      " configuration from the input file \"" + input_file_name + "\"" );
  }

  // Construct a Generator object with the prior configuration. This will be
  // used by the weight calculators as needed
  auto prior_json_config = marley::JSON::load( prior_config_str->value() );
  marley::JSONConfig jc( prior_json_config );
  auto gen = std::make_unique< marley::Generator >( jc.create_generator() );

  // We're almost ready to go. Parse the configuration and set up the output
  // file(s) based on what the user has requested
  std::vector< std::shared_ptr<marley::OutputFile> > output_files;

  if ( rw_config.has_key("output") ) {
    marley::JSON output_set = rw_config.at( "output" );
    if ( !output_set.is_array() ) throw marley::Error( "The"
      " \"output\" key in the reweighting configuration must have a value that"
      " is a JSON array." );
    else for ( const auto& el : output_set.array_range() ) {
      if ( el.has_key("mode") ) {
        std::string mode_str = el.at( "mode" ).to_string();
        if ( mode_str != "overwrite" ) throw marley::Error( "Only the"
          " \"overwrite\" output file mode is allowed for a reweighting job." );
      }
      output_files.push_back( marley::OutputFile::make_OutputFile(el) );
    }
  }
  else {
    // If the user didn't specify anything for the output key, then
    // by default write to a single ASCII-format file.
    std::string out_config_str = "{ format: \"ascii\","
      " file: \"reweighted_events.hepmc3\", mode: \"overwrite\" }";
    auto out_config = marley::JSON::load( out_config_str );

    output_files.push_back( marley::OutputFile::make_OutputFile(out_config) );
  }

  // Now calculate the actual weights for each event. Use a "do" loop so that
  // we can start with the event that was already read in
  int event_count = 0;
  do {

    std::cout << "Event " << event_count << '\n';

    // Update the current list of weight names in the run information to
    // include the new ones
    run_info->set_weight_names( full_name_vec );

    // Append initial unit weights to the current event for each of the new
    // weights to be calculated in this reweighting job. These will be
    // multiplied by the results returned by the weight calculators when the
    // event is processed below.
    auto& ev_wgt_vec = ev.weights();
    for ( size_t w = 0u; w < num_new_weights; ++w ) {
      ev_wgt_vec.push_back( 1. );
    }

    // Process the event weights for the current event
    weighter.process_event( ev, *gen );

    // Write the updated events (including the new weights) to the output
    // file(s)
    for ( const auto& file : output_files ) {
      file->write_event( &ev );
    }

    // Set the weight names back to the old list to avoid problems with
    // reading in the new event
    run_info->set_weight_names( wgt_names );

    ++event_count;

  } while ( efr >> ev );

  return 0;
}
