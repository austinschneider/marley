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
#include <limits>
#include <string>

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenParticle.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
// #include "TDirectory.h"
// #include "TVersionCheck.h"
// #include "TClass.h"
// #include "TObject.h"
// #include "TStorage.h"
// #include "TDataType.h"

// MARLEY includes
#include "marley/OutputFilePlainRoot.hh"
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/JSONConfig.hh"
#include "marley/hepmc3_utils.hh"
#include "marley/marley_utils.hh"
#include "marley/Logger.hh"

#ifdef USE_ROOT
  #include "marley/RootJSONConfig.hh"
#endif
#ifndef USE_ROOT
  #error TEMPORARY: Building this class (OutputFilePlainRoot) requires ROOT
#endif

namespace {

  constexpr char DUMMY_CHAR = 'a';

  // Moves to the end of an ASCII-format file containing HepMC3 events and
  // opened as the input std::fstream. Backs up until the start of the last
  // HepMC3 event stored in the file. The stream is left in a state that is
  // ready for reading in the last event using, e.g., HepMC3::ReaderAscii.
  void seek_to_last_genevent( std::fstream& stream ) {
    char c = DUMMY_CHAR;
    char old_c = DUMMY_CHAR;
    stream.seekg( 0, std::ios::end );
    std::streampos size = stream.tellg();
    for ( int i = 1; i <= size; ++i ) {
      stream.seekg( -i, std::ios::end );
      old_c = c;
      stream.get( c );
      if ( c == '\n' && old_c == 'E' ) {
        stream.seekg( -i + 1, std::ios::end );
        stream.clear();
        break;
      }
    }
  }

}

marley::OutputFilePlainRoot::OutputFilePlainRoot( const marley::JSON& output_config )
  : marley::OutputFile( output_config )
{
  format_ = Format::PLAIN_ROOT;

  // As established this isn't really opening anything
  this->open();

  // Create the TFile and TTree pointers
  out_tfile_ = new TFile( name_.c_str(), "recreate" );
  out_tree_ = new TTree( "mst", "MARLEY summary tree" );
  // Set the current ROOT directory to the TFile
  out_tfile_->cd();
 // Autosave after 100 kB of data is filled
 // (done to get updates on the TTree size before it is written to disk)
  //out_tree_->SetAutoSave(10000);

  // Define the branches
  // projectile branches
  out_tree_->Branch( "pdgv", &pdgv_, "pdgv/I" );
  out_tree_->Branch( "Ev", &Ev_, "Ev/D" );
  out_tree_->Branch( "KEv", &KEv_, "KEv/D" );
  out_tree_->Branch( "pxv", &pxv_, "pxv/D" );
  out_tree_->Branch( "pyv", &pyv_, "pyv/D" );
  out_tree_->Branch( "pzv", &pzv_, "pzv/D" );

  // target branches
  out_tree_->Branch( "pdgt", &pdgt_, "pdgt/I" );
  out_tree_->Branch( "Mt", &Mt_, "Mt/D" );

  // ejectile branches
  out_tree_->Branch( "pdgl", &pdgl_, "pdgl/I" );
  out_tree_->Branch( "El", &El_, "El/D" );
  out_tree_->Branch( "KEl", &KEl_, "KEl/D" );
  out_tree_->Branch( "pxl", &pxl_, "pxl/D" );
  out_tree_->Branch( "pyl", &pyl_, "pyl/D" );
  out_tree_->Branch( "pzl", &pzl_, "pzl/D" );

  // residue branches
  out_tree_->Branch( "pdgr", &pdgr_, "pdgr/I" );
  out_tree_->Branch( "Er", &Er_, "Er/D" );
  out_tree_->Branch( "KEr", &KEr_, "KEr/D" );
  out_tree_->Branch( "pxr", &pxr_, "pxr/D" );
  out_tree_->Branch( "pyr", &pyr_, "pyr/D" );
  out_tree_->Branch( "pzr", &pzr_, "pzr/D" );

  // Nuclear excitation energy branch
  out_tree_->Branch( "Ex", &Ex_, "Ex/D" );

  // Spin and parity branches
  out_tree_->Branch( "twoJ", &twoJ_, "twoJ/I" );
  out_tree_->Branch( "parity", &par_, "parity/I" );

  // De-excitation products (final-state particles other than the
  // ejectile and ground-state residue)
  out_tree_->Branch( "np", &np_, "np/I" );
  out_tree_->Branch( "pdgp", PDGs_.data(), "pdgp[np]/I" );
  out_tree_->Branch( "Ep",  Es_.data(), "Ep[np]/D" );
  out_tree_->Branch( "KEp", KEs_.data(), "KEp[np]/D" );
  out_tree_->Branch( "pxp", pXs_.data(), "pxp[np]/D" );
  out_tree_->Branch( "pyp", pYs_.data(), "pyp[np]/D" );
  out_tree_->Branch( "pzp", pZs_.data(), "pzp[np]/D" );

  // Flux-averaged total cross section
  out_tree_->Branch( "xsec", &flux_avg_tot_xsec_, "xsec/D" );

  // Clear the vectors (this is probably useless but here we go)
  PDGs_.clear();
  Es_.clear();
  KEs_.clear();
  pXs_.clear();
  pYs_.clear();
  pZs_.clear();

}

/// Implement the destructor
marley::OutputFilePlainRoot::~OutputFilePlainRoot() {
  // Write the TTree to the TFile
  out_tree_->Write();
  // Close the TFile
  out_tfile_->Close();
}

/// @todo This doesn't actually "open" anything, as in the ASCII case.
/// Maintaining the name for consistency with the ASCII case,
/// but this should be changed in the future.
void marley::OutputFilePlainRoot::open() {

  bool file_exists = check_if_file_exists( name_ );

  if ( mode_ == Mode::OVERWRITE && file_exists && !force_ ) {
    bool overwrite = marley_utils::prompt_yes_no( "Overwrite file " + name_ );
    if ( !overwrite ) {
      MARLEY_LOG_INFO() << "Cancelling overwrite of output file \""
        << name_ << '\"';
      mode_ = Mode::RESUME;
    }
  }

  if ( mode_ == Mode::RESUME ) {
    if ( !file_exists ) throw marley::Error( "Cannot resume run. Could"
      " not open the file \"" + name_ + '\"' );
  }
  else if ( mode_ != Mode::OVERWRITE ) {
    throw marley::Error( "Unrecognized file mode encountered in"
      " OutputFilePlainRoot::open()" );
  }

}

/// @todo Implement this function
bool marley::OutputFilePlainRoot::resume( std::unique_ptr<marley::Generator>& gen,
  long& num_previous_events )
{
  if ( mode_ != Mode::RESUME ) {
    throw marley::Error( "Cannot call OutputFilePlainRoot::resume() for an output"
      " mode other than \"resume\"" );
    return false;
  }

  return false;
}
 
int_fast64_t marley::OutputFilePlainRoot::bytes_written() {
  // If the ROOT file is open, then update the byte count.
  // We will use the byte count from the TTree as a (pretty accurate) estimate
  // for the file size as we don't write to the file until
  // the destructor is called.
  // Once the file is closed, the byte count will be fixed.
  // @todo: this last part is bullshit, fix it
  if ( out_tfile_->IsOpen() ) {
    byte_count_ = out_tree_->GetZipBytes();
  } else {
    std::filesystem::path p( name_ );
    byte_count_ = std::filesystem::file_size( p );
  }
  return byte_count_;
}

void marley::OutputFilePlainRoot::write_event( HepMC3::GenEvent* ev ) {

  if ( !ev ) throw marley::Error( "Null pointer passed to"
    " OutputFilePlainRoot::write_event()" );

  // Instead of writing the event into ASCII, we "convert" it and write
  // it in plain ROOT format.
  auto projectile = marley_hepmc3::get_projectile( *ev );
  pdgv_ = projectile->pid();

  double mv = projectile->generated_mass();
  const HepMC3::FourVector& p4v = projectile->momentum();

  Ev_ = p4v.e();
  KEv_ = std::max( 0., Ev_ - mv );
  pxv_ = p4v.px();
  pyv_ = p4v.py();
  pzv_ = p4v.pz();

  auto target = marley_hepmc3::get_target( *ev );
  pdgt_ = target->pid();
  Mt_ = target->generated_mass();

  auto ejectile = marley_hepmc3::get_ejectile( *ev );
  pdgl_ = ejectile->pid();

  // Object ID number for the ejectile in the event (*not* the same as the
  // PDG code)
  int ej_id = ejectile->id();

  double ml = ejectile->generated_mass();
  const HepMC3::FourVector& p4l = ejectile->momentum();

  El_ = p4l.e();
  KEl_ = std::max( 0., El_ - ml );
  pxl_ = p4l.px();
  pyl_ = p4l.py();
  pzl_ = p4l.pz();

  auto residue = marley_hepmc3::get_residue( *ev );
  pdgr_ = residue->pid();

  double mr = residue->generated_mass();
  const HepMC3::FourVector& p4r = residue->momentum();

  Er_ = p4r.e();
  KEr_ = std::max( 0., Er_ - mr );
  pxr_ = p4r.px();
  pyr_ = p4r.py();
  pzr_ = p4r.pz();

  // TODO: add error handling here for missing attributes
  auto Ex_attr = residue->attribute< HepMC3::DoubleAttribute >( "Ex" );
  Ex_ = Ex_attr->value();

  auto twoJ_attr = residue->attribute< HepMC3::IntAttribute >( "twoJ" );
  twoJ_ = twoJ_attr->value();

  auto parity_attr = residue->attribute< HepMC3::IntAttribute >( "parity" );
  par_ = parity_attr->value();

  /// @todo Move this to the destructor (as you only need it read once)
  auto run_info = ev->run_info();
  auto avg_xsec_attr = run_info->attribute< HepMC3::DoubleAttribute >(
      "NuHepMC.FluxAveragedTotalCrossSection" );
  flux_avg_tot_xsec_ = avg_xsec_attr->value()
      / marley_utils::fm2_to_picobarn
      * marley_utils::fm2_to_minus40_cm2 * 1e2;

  np_ = 0;
  const auto& particles = ev->particles();
  for ( const auto& p : particles ) {
    // Only save information about final-state particles in the array of
    // nuclear de-excitation products
    if ( p->status() != marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS ) {
      continue;
    }

    // Skip the ejectile since its information is already saved elsewhere in
    // the output branches. We use the unique object ID number in the event
    // (not the same as the PDG code) to check for this
    if ( p->id() == ej_id ) continue;

    ++np_;

    PDGs_.push_back( p->pid() );

    double mp = p->generated_mass();
    const HepMC3::FourVector& p4p = p->momentum();

    double Ep = p4p.e();
    Es_.push_back( Ep );

    double KEp = std::max( 0., Ep - mp );
    KEs_.push_back( KEp );

    pXs_.push_back( p4p.px() );
    pYs_.push_back( p4p.py() );
    pZs_.push_back( p4p.pz() );
  }

  // Update the branch addresses (manipulating the vectors may have
  // invalidated them)
  out_tree_->SetBranchAddress( "pdgp", PDGs_.data() );
  out_tree_->SetBranchAddress( "Ep",  Es_.data() );
  out_tree_->SetBranchAddress( "KEp", KEs_.data() );
  out_tree_->SetBranchAddress( "pxp", pXs_.data() );
  out_tree_->SetBranchAddress( "pyp", pYs_.data() );
  out_tree_->SetBranchAddress( "pzp", pZs_.data() );

  out_tree_->Fill();
  // MARLEY_LOG_INFO() << "Event written to the ROOT file with total bytes " << out_tree_->GetTotBytes();
  // MARLEY_LOG_INFO() << "Event written to the ROOT file with total Zipbytes " << out_tree_->GetZipBytes();
}