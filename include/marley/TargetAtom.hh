/// @file
/// @copyright Copyright (C) 2016-2024 Steven Gardiner
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
#include <ostream>
#include <string>

namespace marley {

  /// @brief An atomic target for a lepton scattering reaction
  class TargetAtom {

    public:

      /// @brief Create a TargetAtom object from an integer PDG code
      /// @param pdg Particle Data Group code for the atomic species
      explicit TargetAtom( int pdg );

      /// @brief Create a TargetAtom object from a proton number (Z)
      /// and mass number (A)
      /// @param Z The proton number of the target atom
      /// @param A The mass number of the target atom
      explicit TargetAtom( int Z, int A );

      /// @brief Convert the TargetAtom object to an integer PDG code
      /// @todo Make this "constexpr" eventually. There is a bug in GCC before
      /// version 7.2 that doesn't respect the C++14 standard. Clang 3.6
      /// correctly recognizes use of "constexpr" for this function as valid.
      /// See https://stackoverflow.com/a/36489493 for more details. For now,
      /// we'll drop "constexpr" to preserve compatibility with old versions of
      /// GCC.
      inline explicit operator int() const { return pdg_; }

      /// @brief Copy assignment operator
      inline marley::TargetAtom& operator=( const marley::TargetAtom& ta )
      {
        pdg_ = ta.pdg_;
        return *this;
      }

      /// @brief Two TargetAtom objects are considered equal if their
      /// nuclear PDG codes match
      inline bool operator==(const marley::TargetAtom& ta) const
      {
        return pdg_ == ta.pdg_;
      }

      inline bool operator!=(const marley::TargetAtom& ta) const
      {
        return pdg_ != ta.pdg_;
      }

      /// @brief Define the less-than operator so that TargetAtom objects
      /// can be used as the keys of a std::map
      inline bool operator<(const marley::TargetAtom& ta) const
      {
        return pdg_ < ta.pdg_;
      }

      /// @brief Converts the PDG code to a string representation (e.g.,
      /// "40Ar")
      std::string to_string() const;

      /// @brief Returns the proton number of the target atom
      int Z() const;

      /// @brief Returns the mass number of the target atom
      int A() const;

      /// @brief Returns the nuclear PDG code of the target atom
      inline int pdg() const { return pdg_; }

    protected:

      /// @brief Helper function for the constructors. Checks the atom's PDG
      /// code for validity.
      /// @details If the pdg_ member variable does not correspond to an
      /// atom (or a free neutron) when this function is called, then a
      /// marley::Error will be thrown.
      void check_pdg_validity() const;

      /// @brief PDG code for the nucleus of the target atom
      int pdg_;
  };

}

inline std::ostream& operator<<(std::ostream& out, const marley::TargetAtom& ta)
{
  out << ta.to_string();
  return out;
}
