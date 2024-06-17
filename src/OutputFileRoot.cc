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

// HepMC3 includes
#include "HepMC3/Data/GenEventData.h"
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"

// MARLEY includes
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/OutputFileRoot.hh"
#include "marley/RootJSONConfig.hh"

marley::OutputFileRoot::OutputFileRoot( const marley::JSON& output_config )
  : marley::OutputFile( output_config )
{
  this->open();
}

marley::OutputFileRoot::~OutputFileRoot() {
  out_tfile_->cd();
  out_tree_->Write();
}

void marley::OutputFileRoot::open() {

  bool file_exists = check_if_file_exists( name_ );
  std::string open_mode_str( "update" );
  if ( mode_ == Mode::OVERWRITE ) open_mode_str = "recreate";

  // TODO: check if TTree already exists in the file
  out_tfile_ = std::make_unique< TFile >( name_.c_str(),
    open_mode_str.c_str() );
  out_tree_ = new TTree( "MARLEY_event_tree", "HepMC3-format MARLEY events" );

  temp_event_data_ = std::make_unique< HepMC3::GenEventData >();
  out_tree_->Branch( "event", temp_event_data_.get() );
}

bool marley::OutputFileRoot::resume(
  std::unique_ptr<marley::Generator>& /*gen*/, long& /*num_previous_events*/ )
{
  if ( mode_ != Mode::RESUME ) {
    throw marley::Error( "Cannot call OutputFileRoot::resume() for an output"
      " mode other than \"resume\"" );
    return false;
  }

  return false;
}

void marley::OutputFileRoot::write_event( const HepMC3::GenEvent* event ) {

  if ( !event ) throw marley::Error( "Null pointer passed to"
    " OutputFileRoot::write_event()" );

  // Clear out any pre-existing contents from the temporary event storage
  temp_event_data_->particles.clear();
  temp_event_data_->vertices.clear();
  temp_event_data_->links1.clear();
  temp_event_data_->links2.clear();
  temp_event_data_->attribute_id.clear();
  temp_event_data_->attribute_name.clear();
  temp_event_data_->attribute_string.clear();

  event->write_data( *temp_event_data_ );
  out_tree_->Fill();
}
