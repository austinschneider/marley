/// @file
/// @copyright Copyright (C) 2016-2023 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

// Standard library includes
#include <array>
#include <iostream>

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"

// MARLEY includes
#include "marley/hepmc3_utils.hh"
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/JSONConfig.hh"
#include "marley/NucleusDecayer.hh"
#include "marley/OutputFile.hh"
#include "marley/Parity.hh"
#include "marley/Reaction.hh"

using ProcType = marley::Reaction::ProcessType;

// Assume that the target is an atom and thus has zero net charge
constexpr int TARGET_NET_CHARGE = 0;

// List of nuclear processes to consider. Other processes will not be
// considered by this program.
constexpr std::array< ProcType, 3 > nuclear_proc_types
  = { ProcType::NeutrinoCC, ProcType::AntiNeutrinoCC, ProcType::NC };

// Helper function that creates an event object containing the nucleus to be
// de-excited
std::shared_ptr< HepMC3::GenEvent > make_event_object(
  int pdg_a, int pdg_b, int pdg_c, ProcType process_type, double Ex,
  int twoJ, const marley::Parity& P,
  std::shared_ptr< HepMC3::GenParticle >& residue, int residue_net_charge )
{
  // NuHepMC E.R.4
  auto event = std::make_shared< HepMC3::GenEvent >( HepMC3::Units::MEV,
    HepMC3::Units::CM );

  // NuHepMC E.R.3
  int signal_process_id = marley_hepmc3::get_nuhepmc_proc_id( process_type );
  event->add_attribute( "signal_process_id",
    std::make_shared< HepMC3::IntAttribute >( signal_process_id )
  );

  // Create the primary vertex
  // NuHepMC E.R.6
  auto prim_vtx = std::make_shared< HepMC3::GenVertex >();
  prim_vtx->set_status( marley_hepmc3::NUHEPMC_PRIMARY_VERTEX );

  event->add_vertex( prim_vtx );

  // Create dummy particle objects for the projectile, target, and ejectile
  auto projectile = marley_hepmc3::make_particle( pdg_a, 0., 0., 0., 0.,
    marley_hepmc3::NUHEPMC_PROJECTILE_STATUS, 0. );

  auto target = marley_hepmc3::make_particle( pdg_b,
    marley_hepmc3::NUHEPMC_TARGET_STATUS, 0. );

  // Create particle objects representing the ejectile and residue in the CM
  // frame.
  auto ejectile = marley_hepmc3::make_particle( pdg_c, 0., 0.,
    0., 0., marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, 0. );

  // Attach the particles to the primary vertex
  prim_vtx->add_particle_in( projectile );
  prim_vtx->add_particle_in( target );

  prim_vtx->add_particle_out( ejectile );
  prim_vtx->add_particle_out( residue );

  // Add attributes needed to keep track of the nuclear de-excitation state
  residue->add_attribute( "Ex",
    std::make_shared< HepMC3::DoubleAttribute >(Ex) );
  residue->add_attribute( "twoJ",
    std::make_shared< HepMC3::IntAttribute >(twoJ) );
  residue->add_attribute( "parity",
    std::make_shared< HepMC3::IntAttribute >(static_cast<int>( P )) );

  // Assume that the target has net charge TARGET_NET_CHARGE
  marley_hepmc3::set_particle_charge( *target, TARGET_NET_CHARGE );

  // Assign the correct charge to the residue
  marley_hepmc3::set_particle_charge( *residue, residue_net_charge );

  return event;
}

int main( int argc, char* argv[] ) {

  // If the user has not supplied enough command-line arguments, display the
  // standard help message and exit
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " CONFIG_FILE\n";
    return 1;
  }

  // Get the configuration file name from the command line
  std::string config_file_name( argv[1] );

  // Set up the generator using the job configuration file
  marley::JSONConfig jc( config_file_name );

  marley::Generator gen = jc.create_generator();
  gen.set_up_run_info();

  // Parse the extra "decay sim" configuration parameters
  const marley::JSON& json = jc.get_json();

  // Flag to check that JSON parsing went all right
  bool ok;

  // Check that there is an object with the right label in the configuration
  // file
  const std::string decay_config_label( "decays" );

  marley::JSON decays;
  ok = get_from_json< marley::JSON >( decay_config_label, json, decays );
  if ( !ok ) throw marley::Error( "Missing key '" + decay_config_label
    + "' in job configuration file" );

  // Get the total number of events to be simulated (default to 1e3 events)
  long num_events = assign_from_json< long >( "events", decays, ok, 1000 );

  // Get the projectile PDG code. If it is missing, assume it's an electron
  // neutrino.
  int projectile_pdg = assign_from_json< int >( "projectile", decays,
    ok, marley_utils::ELECTRON_NEUTRINO );

  // Get the target proton and nucleon numbers
  int Zi = assign_from_json< int >( "target_Z", decays, ok );
  int Ai = assign_from_json< int >( "target_A", decays, ok );

  // Get the process type for the primary 2 --> 2 reaction to assume.
  // If it is absent, assume a NC process.
  auto proc_type = static_cast< ProcType >(
    assign_from_json< int >( "proc_type", decays, ok,
      static_cast<int>(ProcType::NC) )
  );

  auto iter = std::find( nuclear_proc_types.cbegin(),
    nuclear_proc_types.cend(), proc_type );

  bool process_is_a_nuclear_reaction = ( iter != nuclear_proc_types.cend() );

  if ( !process_is_a_nuclear_reaction ) {
    std::cerr << "This program handles nuclear reactions only\n";
    return 2;
  }

  // Decide whether Ex will be sampled for each event or fixed
  double Ex = 0.; // Excitation energy for the decay event
  double Ex_min = 0.;
  double Ex_max = 0.;
  bool sample_Ex = decays.has_key( "Ex_max" );
  if ( !sample_Ex ) {
    // Get the initial excitation energy
    Ex = assign_from_json< double >( "Ex", decays, ok, -1.0 );

  }
  else {
    Ex_min = assign_from_json< double >( "Ex_min", decays, ok, -1.0 );
    Ex_max = assign_from_json< double >( "Ex_max", decays, ok, -1.0 );
  }

  if ( Ex < 0. || Ex_min < 0. || Ex_max < 0. ) {
    throw marley::Error( "Negative excitation energy encountered" );
  }
  if ( Ex_max < Ex_min ) {
    throw marley::Error( "Upper Ex bound is less than lower Ex bound" );
  }

  // Get two times the initial nuclear spin
  std::vector< int > twoJ_vec;
  if ( !decays.has_key("twoJ") ) {
    throw marley::Error( "Missing \"twoJ\" key specifying the nuclear spin" );
  }
  else {
    const marley::JSON& twoJ_obj = decays.at( "twoJ" );
    if ( !twoJ_obj.is_array() ) {
     throw marley::Error( "The \"twoJ\" key must have a value that"
       " is a JSON array." );
    }
    else {
      convert_json< std::vector<int> >( twoJ_obj, twoJ_vec );
    }
  }

  for ( const auto& twoJ : twoJ_vec ) {
    if ( twoJ < 0 ) throw marley::Error( "Negative nuclear spin"
      " value encountered" );
  }

  // Sample twoJ values with equal probability
  std::vector< double > twoJ_sampling_weights( twoJ_vec.size(), 1. );
  const auto twoJ_begin = twoJ_sampling_weights.cbegin();
  const auto twoJ_end = twoJ_sampling_weights.cend();
  std::discrete_distribution< size_t > twoJ_dist( twoJ_begin, twoJ_end );

  // Get the initial nuclear parity
  auto parity_str = assign_from_json< std::string >( "parity", decays,
    ok, "+" );
  if ( parity_str != "+" && parity_str != "-" && parity_str != "random" ) {
    throw marley::Error( "Invalid parity setting \"" + parity_str + '\"' );
  }

  // Set up a discrete distribution for sampling parities (in case we need it)
  const std::vector< std::string > parity_strings = { "+", "-" };
  const std::vector< double > parity_sampling_weights = { 1., 1. };
  const auto par_begin = parity_sampling_weights.cbegin();
  const auto par_end = parity_sampling_weights.cend();
  std::discrete_distribution< size_t > par_dist( par_begin, par_end );

  // Get the name to use for the output file
  auto output_file_name = assign_from_json< std::string >( "output_file_name",
    decays, ok );

///////////////////////////////////

  // Determine the final-state particle PDG codes based on the process type
  int ejectile_pdg = marley::Reaction
    ::get_ejectile_pdg( projectile_pdg, proc_type );

  // Use conservation of electric charge to obtain the residue net charge
  int projectile_charge = marley_utils::get_particle_charge( projectile_pdg );
  int ejectile_charge = marley_utils::get_particle_charge( ejectile_pdg );

  // Determine the residue nucleon and proton numbers
  int Af = Ai;
  int Delta_Z = projectile_charge - ejectile_charge;
  int Zf = Zi + Delta_Z;

  int residue_net_charge = Delta_Z + TARGET_NET_CHARGE;

  // Prepare the output file(s)
  std::vector< std::shared_ptr<marley::OutputFile> > output_files;

  if ( decays.has_key("output") ) {
    marley::JSON output_set = decays.at( "output" );
    if ( !output_set.is_array() ) throw marley::Error( "The"
      " \"output\" key must have a value that is a JSON array." );
    else for ( const auto& el : output_set.array_range() ) {
      output_files.push_back( marley::OutputFile::make_OutputFile(el) );
    }
  }
  else {
    // If the user didn't specify anything for the output key, then
    // by default write to a single ASCII-format file.
    std::string out_config_str = "{ format: \"ascii\","
      " file: \"decay_events.hepmc3\", mode: \"overwrite\" }";
    auto out_config = marley::JSON::load( out_config_str );

    output_files.push_back( marley::OutputFile::make_OutputFile(out_config) );
  }

  for ( long evnum = 0; evnum < num_events; ++evnum ) {

    // Get access to the MassTable
    const auto& mt = marley::MassTable::Instance();

    // If excitation energy sampling is enabled, then choose a value uniformly
    // between the bounds
    if ( sample_Ex ) {
      Ex = gen.uniform_random_double( Ex_min, Ex_max, true );
    }

    // Sample an initial nuclear spin value for the current event
    size_t twoJ_index = gen.sample_from_distribution( twoJ_dist );
    int twoJ = twoJ_vec.at( twoJ_index );

    // Choose the initial nuclear parity, sampling a random value if needed
    std::string par_str = parity_str;
    if ( par_str == "random" ) {
      size_t par_index = gen.sample_from_distribution( par_dist );
      par_str = parity_strings.at( par_index );
    }
    marley::Parity parity;
    std::istringstream temp_iss( par_str );
    temp_iss >> parity;

    // If the requested excitation energy corresponds to a bound nuclear
    // state, match it to a known discrete level instead of taking
    // the input at face value
    // TODO: Revisit this. Maybe do something more sophisticated.
    if ( Ex <= mt.unbound_threshold(Zf, Af) ) {
      // Try to load a discrete level DecayScheme from the generator's owned
      // StructureDatabase. If we can find one, use it to assign a refined
      // excitation energy value. If not, give up and move on.
      auto* ds = gen.get_structure_db().get_decay_scheme( Zf, Af );

      if ( ds ) {
        // Replace our current excitation energy value with the closest
        // known discrete level
        auto* lev = ds->get_pointer_to_closest_level( Ex );
        // TODO: perhaps check that the spin-parity of the level matches
        // the expected one, issue a warning if it doesn't?
        Ex = lev->energy();
      }
    }

    int target_pdg = marley_utils::get_nucleus_pid( Zi, Ai );
    int residue_pdg = marley_utils::get_nucleus_pid( Zf, Af );

    // Approximate the ground-state mass of the (possibly ionized) residue by
    // subtracting the appropriate number of electron masses from its atomic
    // (i.e., neutral) ground state mass.
    double residue_gs_mass = mt.get_atomic_mass( residue_pdg )
      - ( residue_net_charge * mt.get_particle_mass(marley_utils::ELECTRON) );

    // Create a Particle object representing the residue. Since we're
    // simulating its decays without worrying about a primary 2 --> 2
    // interaction, let it be at rest in the lab frame.
    double m_residue = residue_gs_mass + Ex;

    auto residue = marley_hepmc3::make_particle( residue_pdg, 0., 0., 0.,
      m_residue, marley_hepmc3::NUHEPMC_UNDECAYED_RESIDUE_STATUS, m_residue );

    // Create the event object with a dummy projectile, target, and ejectile.
    auto event = make_event_object( projectile_pdg, target_pdg, ejectile_pdg,
      proc_type, Ex, twoJ, parity, residue, residue_net_charge );

    // Pass the event to the nuclear de-excitation simulation
    marley::NucleusDecayer nd;
    nd.process_event( *event, gen );

    // Add a little metadata to the event (attach run info, etc.)
    gen.finish_event_metadata( *event );

    // We're done, write the event to the output file(s)
    for ( const auto& file : output_files ) {
      file->write_event( event.get() );
    }

    std::cout << "Event " << evnum << "\n";
  }

}
