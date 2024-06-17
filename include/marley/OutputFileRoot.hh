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

#pragma once

//#include "HepMC3/Data/GenEventData.h"
//#include "HepMC3/Data/GenRunInfoData.h"
//// Needed for generating ROOT dictionaries for HepMC3
//namespace HepMC
//{
//typedef HepMC3::GenEventData GenEventData;
//typedef HepMC3::GenRunInfoData GenRunInfoData ;
//typedef HepMC3::GenParticleData GenParticleData ;
//typedef HepMC3::GenVertexData GenVertexData ;
//typedef HepMC3::FourVector FourVector ;
//typedef HepMC3::Units Units ;
//typedef HepMC3::Setup Setup ;
//}

// Standard library includes
#include <memory>

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// MARLEY includes
#include "marley/OutputFile.hh"

namespace HepMC3 {
  class GenEventData;
};

namespace marley {

  /// @brief Class that represents an output file in a ROOT-based binary format
  /// for HepMC3 events
  class OutputFileRoot : public OutputFile {

    public:

      OutputFileRoot( const JSON& output_config );

      virtual ~OutputFileRoot();

      virtual bool resume( std::unique_ptr<marley::Generator>& gen,
        long& num_previous_events ) override;

      inline virtual int_fast64_t bytes_written() override
        { return out_tfile_->GetBytesWritten(); }

      virtual void write_event( const HepMC3::GenEvent* event ) override;

    protected:

      /// Opens the owned ROOT file with the correct settings
      virtual void open();

      /// TFile object used to manage input/output
      std::unique_ptr< TFile > out_tfile_;

      /// TTree object used to store the events
      // This is a bare pointer, but ROOT will associate it with out_tfile_, so
      // we don't want to delete it ourselves or let a smart pointer do it.
      TTree* out_tree_ = nullptr;

      /// Temporary storage for the events
      std::unique_ptr< HepMC3::GenEventData > temp_event_data_;
  };

}
