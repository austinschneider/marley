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

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/BackshiftedFermiGasModel.hh"
#include "marley/Error.hh"
#include "marley/FileManager.hh"
#include "marley/Fragment.hh"
#include "marley/JSON.hh"
#include "marley/KoningDelarocheOpticalModel.hh"
#include "marley/Logger.hh"
#include "marley/StandardLorentzianModel.hh"
#include "marley/StructureDatabase.hh"
#include "marley/TargetAtom.hh"

namespace {
  /// @brief JSON key for lookup of the default optical model parameters
  const std::string DEFAULT_OM_KEY( "Default" );
}

// Define static data members of the StructureDatabase class
std::map< int, std::pair<int, marley::Parity> >
  marley::StructureDatabase::jpi_table_;

std::map<int, marley::Fragment> marley::StructureDatabase::fragment_table_;

// Name of the data file which will be used to read in the ground-state nuclear
// spin-parity values
const std::string marley::StructureDatabase
  ::jpi_data_file_name_ = "gs_spin_parity_table.txt";

// Flag indicating whether the ground-state spin-parity data file has already
// been loaded
bool marley::StructureDatabase::initialized_gs_spin_parity_table_ = false;

marley::StructureDatabase::StructureDatabase() {}

void marley::StructureDatabase::add_decay_scheme( int pdg,
  std::unique_ptr<marley::DecayScheme>& ds )
{
  auto* temp_ptr = ds.release();
  decay_scheme_table_.emplace(pdg, std::unique_ptr<marley::DecayScheme>(temp_ptr));
}

void marley::StructureDatabase::emplace_decay_scheme(int pdg,
  const std::string& filename, DecayScheme::FileFormat format)
{
  int Z_ds = (pdg % 10000000)/10000;
  int A_ds = (pdg % 10000)/10;

  // Remove the previous entry (if one exists) for the given PDG code
  decay_scheme_table_.erase(pdg);

  // Add the new entry
  decay_scheme_table_.emplace(pdg, std::make_unique<marley::DecayScheme>(
    Z_ds, A_ds, filename, format));
}

std::set<int> marley::StructureDatabase::find_all_nuclides(
  const std::string& filename, DecayScheme::FileFormat format)
{
  if (format != DecayScheme::FileFormat::talys)
    throw marley::Error(std::string("StructureDatabase::")
      + "find_all_nuclides() is not implemented for"
      + " formats other than TALYS.");

  // First line in a TALYS level dataset has fortran
  // format (2i4, 2i5, 56x, i4, a2)
  // General regex for this line:
  static const std::regex nuclide_line("[0-9 ]{18} {56}[0-9 ]{4}.{2}");

  // Open the TALYs level data file for parsing
  std::ifstream file_in(filename);

  // If the file doesn't exist or some other error
  // occurred, complain and give up.
  if (!file_in.good()) throw marley::Error(std::string("Could not")
    + " read from the TALYS data file " + filename);

  std::string line; // String to store the current line
                    // of the TALYS file during parsing

  std::istringstream iss; // String stream used to parse the line

  double Z; // atomic number
  double A; // mass number

  // Particle Data Group codes for each nuclide in the file
  std::set<int> PDGs;

  // Loop through the data file, recording all nuclide PDGs found
  while (std::getline(file_in, line)) {
    if (std::regex_match(line, nuclide_line)) {
      // Load the new line into our istringstream object for parsing. Reset
      // the stream so that we start parsing from the beginning of the string.
      iss.str(line);
      iss.clear();

      // The first two entries on a TALYS nuclide line are Z and A.
      iss >> Z >> A;

      PDGs.insert(marley_utils::get_nucleus_pid(Z, A));
    }
  }

  file_in.close();

  return PDGs;
}

marley::DecayScheme* marley::StructureDatabase::get_decay_scheme(
  const int particle_id)
{
  // If we already have the DecayScheme object stored in the lookup table,
  // then just retrieve it
  auto iter = decay_scheme_table_.find( particle_id );
  if ( iter != decay_scheme_table_.end() ) return iter->second.get();
  // If not, see if a data file for it is listed in the index
  else {

    marley::TargetAtom ta_requested( particle_id );
    MARLEY_LOG_DEBUG() << "Looking up structure data for " << ta_requested;

    if ( !loaded_structure_index_ ) this->load_structure_index();
    auto ds_file_iter = decay_scheme_filenames_.find( particle_id );

    // If not, then just give up and return a null pointer
    if ( ds_file_iter == decay_scheme_filenames_.end() ) return nullptr;

    // If a file is available, load all of the decay schemes present in it
    // and add them to the lookup table. If we find the one we're looking
    // for, return a pointer to it. Otherwise, print a warning, give up,
    // and return a null pointer.
    std::string ds_file_name = ds_file_iter->second;
    auto& fm = marley::FileManager::Instance();
    std::string full_ds_file_name = fm.find_file( ds_file_name );
    std::ifstream ds_data_file( full_ds_file_name );
    bool found_it = false;
    auto temp_ds = std::make_unique< marley::DecayScheme >();
    int loaded_nuclide_count = 0;
    while ( ds_data_file >> *temp_ds ) {
      int ds_pdg = temp_ds->pdg();
      if ( particle_id == ds_pdg ) found_it = true;
      this->add_decay_scheme( ds_pdg, temp_ds );
      marley::TargetAtom ta( ds_pdg );
      MARLEY_LOG_DEBUG() << "Added decay scheme for " << ta << " from "
        << full_ds_file_name;
      ++loaded_nuclide_count;
      temp_ds = std::make_unique< marley::DecayScheme >();
    }
    if ( !found_it ) {
      MARLEY_LOG_WARNING() << "Failed to load nuclear structure"
        << " data for " << ta_requested << " from the"
        << " file " << ds_file_name;
      // Make a nullptr entry in the lookup table to avoid duplicate attempts
      // to load the missing data
      decay_scheme_table_[ particle_id ] = nullptr;
    }
    if ( loaded_nuclide_count > 0 ) {
      MARLEY_LOG_INFO() << "Loaded structure data for "
        << loaded_nuclide_count << " nuclides from the file "
        << full_ds_file_name;
    }

    // Return the pointer to the newly-loaded DecayScheme (or nullptr
    // if it could not be loaded) via recursion
    return this->get_decay_scheme( particle_id );
  }
}

marley::DecayScheme* marley::StructureDatabase::get_decay_scheme(const int Z,
  const int A)
{
  int particle_id = marley_utils::get_nucleus_pid(Z, A);
  return get_decay_scheme( particle_id );
}

marley::OpticalModel& marley::StructureDatabase::get_optical_model(
  int nucleus_pid )
{
  /// @todo add check for invalid nucleus particle ID value
  auto iter = optical_model_table_.find( nucleus_pid );

  if ( iter == optical_model_table_.end() ) {
    // The requested optical model wasn't found, so create it and add it
    // to the table, returning a reference to the stored optical model
    // afterwards.
    int Z = marley_utils::get_particle_Z( nucleus_pid );
    int A = marley_utils::get_particle_A( nucleus_pid );

    if ( om_config_map_.empty() ) this->load_optical_model_params();
    auto om_config = om_config_map_.at( DEFAULT_OM_KEY );

    return *( optical_model_table_.emplace( nucleus_pid,
      std::make_unique< marley::KoningDelarocheOpticalModel >(
      Z, A, om_config) ).first->second.get() );
  }
  else return *( iter->second.get() );
}

marley::OpticalModel& marley::StructureDatabase::get_optical_model(
  const int Z, const int A )
{
  int nucleus_pid = marley_utils::get_nucleus_pid( Z, A );
  auto iter = optical_model_table_.find( nucleus_pid );

  if ( iter == optical_model_table_.end() ) {

    if ( om_config_map_.empty() ) this->load_optical_model_params();
    auto om_config = om_config_map_.at( DEFAULT_OM_KEY );

    // The requested optical model wasn't found, so create it and add it
    // to the table, returning a reference to the stored optical model
    // afterwards.
    return *( optical_model_table_.emplace( nucleus_pid,
      std::make_unique< marley::KoningDelarocheOpticalModel >(
      Z, A, om_config) ).first->second.get() );
  }
  else return *( iter->second.get() );
}

marley::LevelDensityModel& marley::StructureDatabase::get_level_density_model(
  int nucleus_pid)
{
  auto iter = level_density_table_.find(nucleus_pid);

  if (iter == level_density_table_.end()) {
    // The requested level density model wasn't found, so create it and add it
    // to the table, returning a reference to the stored level density model
    // afterwards.
    int Z = marley_utils::get_particle_Z( nucleus_pid );
    int A = marley_utils::get_particle_A( nucleus_pid );
    return *(level_density_table_.emplace(nucleus_pid,
      std::make_unique<marley::BackshiftedFermiGasModel>(Z, A)).first
      ->second.get());
  }
  else return *(iter->second.get());
}

marley::LevelDensityModel& marley::StructureDatabase::get_level_density_model(
  const int Z, const int A)
{
  int pid = marley_utils::get_nucleus_pid(Z, A);
  return this->get_level_density_model( pid );
}

marley::GammaStrengthFunctionModel&
  marley::StructureDatabase::get_gamma_strength_function_model(
  const int nuc_pdg)
{
  int Z = marley_utils::get_particle_Z( nuc_pdg );
  int A = marley_utils::get_particle_A( nuc_pdg );
  return this->get_gamma_strength_function_model( Z, A );
}

marley::GammaStrengthFunctionModel&
  marley::StructureDatabase::get_gamma_strength_function_model(const int Z,
  const int A)
{
  int pid = marley_utils::get_nucleus_pid(Z, A);

  auto iter = gamma_strength_function_table_.find(pid);

  if (iter == gamma_strength_function_table_.end()) {
    // The requested gamma-ray strength function model wasn't found, so create
    // it and add it to the table, returning a reference to the stored strength
    // function model afterwards.
    return *(gamma_strength_function_table_.emplace(pid,
      std::make_unique<marley::StandardLorentzianModel>(Z, A)).first
      ->second.get());
  }
  else return *(iter->second.get());
}

void marley::StructureDatabase::remove_decay_scheme(int pdg)
{
  // Remove the decay scheme with this PDG code if it exists in the database.
  // If it doesn't, do nothing.
  decay_scheme_table_.erase( pdg );
}

void marley::StructureDatabase::clear() {
  decay_scheme_table_.clear();
}

const marley::Fragment* marley::StructureDatabase::get_fragment(
  const int fragment_pdg)
{
  // Before retrieving the fragment, make sure we've loaded the
  // necessary data tables
  if ( !initialized_gs_spin_parity_table_ ) initialize_jpi_table();
  auto iter = fragment_table_.find( fragment_pdg );
  if ( iter == fragment_table_.end() ) return nullptr;
  else return &iter->second;
}

const marley::Fragment* marley::StructureDatabase::get_fragment(
  const int Z, const int A)
{
  int fragment_pdg = marley_utils::get_nucleus_pid( Z, A );
  return get_fragment( fragment_pdg );
}

void marley::StructureDatabase::get_gs_spin_parity(
  const int Z, const int A, int& twoJ, marley::Parity& Pi)
{
  int nuc_pdg = marley_utils::get_nucleus_pid( Z, A );
  return get_gs_spin_parity( nuc_pdg, twoJ, Pi );
}

void marley::StructureDatabase::get_gs_spin_parity(int nuc_pdg,
  int& twoJ, marley::Parity& Pi)
{
  if ( !initialized_gs_spin_parity_table_ ) initialize_jpi_table();
  auto iter = jpi_table_.find( nuc_pdg );
  if ( iter == jpi_table_.end() ) throw marley::Error( "Unrecognized"
    " nuclear PDG code " + std::to_string(nuc_pdg) + " passed to"
    " marley::StructureDatabase::get_gs_spin_parity()" );
  auto pair = iter->second;
  twoJ = pair.first;
  Pi = pair.second;
}

void marley::StructureDatabase::initialize_jpi_table() {

  // Instantiate the file manager and use it to find
  // the data file containing the ground-state spin-parities
  // for many nuclei
  const auto& fm = marley::FileManager::Instance();
  std::string full_jpi_file_name
    = fm.find_file( jpi_data_file_name_ );

  if ( full_jpi_file_name.empty() ) {
    throw marley::Error( "Could not find the MARLEY nuclear ground-state"
      " spin-parity data file " + jpi_data_file_name_ + ". Please ensure that"
      " the folder containing it is on the MARLEY search path."
      " If needed, the folder can be appended to the MARLEY_SEARCH_PATH"
      " environment variable." );
  }

  MARLEY_LOG_INFO() << "Loading ground-state nuclear spin-parities from "
    << full_jpi_file_name;

  std::ifstream table_file( full_jpi_file_name );
  int nuc_pdg, twoJ;
  marley::Parity Pi;
  while ( table_file >> nuc_pdg >> twoJ >> Pi ) {
    MARLEY_LOG_DEBUG() << "Nucleus with PDG code " << nuc_pdg
      << " has spin-parity " << static_cast<double>( twoJ ) / 2.
      << Pi;
    jpi_table_[ nuc_pdg ] = std::pair<int, marley::Parity>( twoJ, Pi );
  }

  // Set the flag saying we've initialized the table of ground-state spin-parities.
  // This will avoid duplicate attempts at initialization.
  initialized_gs_spin_parity_table_ = true;

  // Also initialize the table of nuclear fragment properties now that we have
  // the needed information
  using namespace marley_utils;
  constexpr std::array< int, 6 > FRAGMENTS_TO_CONSIDER =
    { NEUTRON, PROTON,  DEUTERON, TRITON,  HELION,  ALPHA };

  for ( int f_pdg : FRAGMENTS_TO_CONSIDER ) {
    // Temporary storage
    int f_twoJ;
    marley::Parity f_Pi;

    // Look up the spin-parity of the nuclear fragment (assumed to be emitted in
    // its ground state) in the data table
    marley::StructureDatabase::get_gs_spin_parity( f_pdg, f_twoJ, f_Pi );

    // Create a new entry in the table of nuclear fragments that should be
    // considered during unbound nuclear de-excitations
    marley::StructureDatabase::fragment_table_.emplace(
      std::make_pair( f_pdg, marley::Fragment(f_pdg, f_twoJ, f_Pi) ));
  }

}

void marley::StructureDatabase::load_structure_index() {

  // Instantiate the file manager and use it to find
  // the index to the decay scheme data files
  const auto& fm = marley::FileManager::Instance();
  std::string full_index_file_name = fm.find_file( structure_index_filename_ );

  if ( full_index_file_name.empty() ) {
    throw marley::Error( "Could not find the MARLEY structure data index file "
      + structure_index_filename_ + ". Please ensure that"
      " the folder containing it is on the MARLEY search path."
      " If needed, the folder can be appended to the MARLEY_SEARCH_PATH"
      " environment variable." );
  }

  MARLEY_LOG_INFO() << "Loaded structure data index from "
    << full_index_file_name;

  std::ifstream index_file( full_index_file_name );
  int nuc_pdg;
  std::string data_file_name;
  while ( index_file >> nuc_pdg >> data_file_name ) {
    MARLEY_LOG_DEBUG() << "Nucleus with PDG code " << nuc_pdg
      << " has a tabulated decay scheme in the file " << data_file_name;
    decay_scheme_filenames_[ nuc_pdg ] = data_file_name;
  }

  // Avoid duplicate loading of the structure index by setting the
  // "already loaded" flag
  loaded_structure_index_ = true;
}

void marley::StructureDatabase::load_optical_model_params(
  const marley::JSON* om_config )
{
  // Boolean flag for checking that JSON parsing was successful
  bool ok;

  // If we've been passed a JSON configuration, then use it
  if ( om_config ) {
    om_config_map_ = assign_from_json< std::map<std::string, marley::JSON> >(
      *om_config, ok );
  }
  // Otherwise, fetch the default one from the standard data file
  else {
    const auto& fm = marley::FileManager::Instance();
    const std::string om_file_name( "optical_model.js" );
    std::string om_full_file_name = fm.find_file( om_file_name );

    if ( om_full_file_name.empty() ) {
      throw marley::Error( "Could not find the MARLEY nuclear optical model"
        " configuration file " + om_file_name + ". Please ensure that"
        " the folder containing it is on the MARLEY search path."
        " If needed, the folder can be appended to the MARLEY_SEARCH_PATH"
        " environment variable." );
    }

    MARLEY_LOG_INFO() << "Loading nuclear optical model parameters from "
      << om_full_file_name;

    auto temp_config = marley::JSON::load_file( om_full_file_name );
    om_config_map_ = assign_from_json< std::map<std::string, marley::JSON> >(
      temp_config, ok );
  }

  if ( !ok ) throw marley::Error( "Failed to parse optical model"
    " JSON configuration" );
}
