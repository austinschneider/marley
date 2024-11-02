/// @file
/// @copyright Copyright (C) 2016-2021 Steven Gardiner
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

// Standard library includes
#include <string>

namespace marley {

  // Forward-declare the Event class
  class Event;

  class MacroEventFileReader {

    public:

      MacroEventFileReader( const std::string& file_name );

      ~MacroEventFileReader();

      bool next_event(marley::Event& ev);

      operator bool() const;

      /// @brief Stream operator for reading in the next event
      inline MacroEventFileReader& operator>>( marley::Event& ev ) {
        next_event( ev );
        return *this;
      }

      double flux_averaged_xsec( bool natural_units = false ) const;

    protected:

      void* event_file_reader_;
  };

}

