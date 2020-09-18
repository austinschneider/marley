/// @file
/// @copyright Copyright (C) 2016-2020 Steven Gardiner
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
#include <fstream>

// MARLEY includes
#include "marley/MassTable.hh"
#include "marley/Reaction.hh"
#include "marley/TabulatedXSec.hh"
#include "marley/marley_utils.hh"

marley::TabulatedXSec::TabulatedXSec( int target_pdg,
  marley::Reaction::ProcessType p_type )
  : ta_( target_pdg ), proc_type_( p_type )
{
}

void marley::TabulatedXSec::add_table( const std::string& file_name )
{
  // TODO: add error handling
  std::ifstream in_file( file_name );

  // Multipole order
  unsigned J;
  in_file >> J;

  // Grid sizes
  unsigned num_w, num_q;

  // Temporary storage for grid points
  double w, q;

  // 1D grids
  auto wvec = std::make_shared< std::vector<double> >();
  auto qvec = std::make_shared< std::vector<double> >();

  in_file >> num_w;
  for ( unsigned iw = 0u; iw < num_w; ++iw ) {
    in_file >> w;
    wvec->push_back( w );
  }

  in_file >> num_q;
  for ( unsigned iq = 0u; iq < num_q; ++iq ) {
    in_file >> q;
    qvec->push_back( q );
  }

  // Natural parity is given by (-1)^J, while unnatural parity is the opposite
  marley::Parity natural;
  bool J_is_odd = (J % 2 == 1);
  if ( J_is_odd ) natural = -1;
  else natural = 1;

  marley::Parity unnatural = -natural;

  MultipoleLabel nat_ml( J, natural );
  MultipoleLabel unnat_ml( J, unnatural );

  auto nat_resp = std::make_shared<
    std::vector<marley::ResponseTable::NuclearResponses> >();

  auto unnat_resp = std::make_shared<
    std::vector<marley::ResponseTable::NuclearResponses> >();

  // Temporary storage for responses
  double rcc, rll, rcl, rt, rtprime;

  // Flag used to swap back and forth between natural and unnatural parity
  bool nat = true;
  while ( in_file >> rcc >> rll >> rcl >> rt >> rtprime ) {

    if ( nat ) nat_resp->emplace_back( rcc, rll, rcl, rt, rtprime );
    else unnat_resp->emplace_back( rcc, rll, rcl, rt, rtprime );

    // Flip the parity flag
    nat = !nat;
  }

  // Store the finished tables as new objects in the map
  responses_[ nat_ml ] = ResponseTable( wvec, qvec, nat_resp );
  responses_[ unnat_ml ] = ResponseTable( wvec, qvec, unnat_resp );

}


double marley::TabulatedXSec::diff_xsec( int pdg_a, double KEa, double omega,
  double cos_theta, const marley::TabulatedXSec::MultipoleLabel& ml )
{
  // TODO: check validity of pdg_a

  // Look up the masses of the projectile and ejectile
  int pdg_c = marley::Reaction::get_ejectile_pdg( pdg_a, proc_type_ );
  const auto& mt = marley::MassTable::Instance();
  double ma = mt.get_particle_mass( pdg_a );
  double mc = mt.get_particle_mass( pdg_c );

  // Determine their total energies and momenta
  double Ea = KEa + ma;
  double pa = marley_utils::real_sqrt( Ea*Ea - ma*ma );
  double Ec = Ea - omega;
  if ( Ec < mc ) return 0.;

  double pc = marley_utils::real_sqrt( Ec*Ec - mc*mc );

  // Get the magnitude of the 3-momentum transfer (q) from this information
  // and the scattering cosine. If the scattering cosine is unphysical, then
  // just return zero.
  if ( std::abs(cos_theta) > 1. ) return 0.;
  double q = marley_utils::real_sqrt( pa*pa + pc*pc - 2.*pa*pc*cos_theta );

  // Get the table of pre-computed nuclear responses for the requested
  // multipole
  const auto& rt = this->responses_.at( ml );

  // Interpolate a set of nuclear responses for the given omega and q values
  if ( omega < rt.w_min() || omega > rt.w_max()
    || q < rt.q_min() || q > rt.q_max() ) return 0.;
  auto nr = rt.interpolate( omega, q );

  // Compute the lepton factors
  double beta = pc / Ec;
  double sin_theta2 = marley_utils::real_sqrt( 1. - cos_theta*cos_theta );
  double q2 = q*q;
  double vcc = 1. + beta * cos_theta;
  double vll = vcc - 2.*Ea*Ec*sin_theta2/q2;
  double vcl = 2. * ( omega*vcc/q + mc*mc/Ec/q );
  double vT = 1. - beta*cos_theta + Ea*Ec*beta*beta*sin_theta2/q2;
  double vTprime = -2. * ( (Ea + Ec)*(1 - beta*cos_theta)/q - mc*mc/q/Ec );
  LeptonFactors lf( vcc, vll, vcl, vT, vTprime );

  // Compute the double-differential cross section with respect to
  // the energy transfer and the scattering cosine
  double xsec = lf * nr;
  xsec *= 2. * marley_utils::GF2 * marley_utils::Vud2 * Ec * pc;
  return xsec;
}

double marley::TabulatedXSec::integral( int pdg_a, double KEa,
  const marley::TabulatedXSec::MultipoleLabel& ml, double& diff_max )
{
  // Set the maximum differential cross section to minus infinity
  // to start
  diff_max = marley_utils::minus_infinity;

  // Get the table of nuclear responses for the requested multipole
  const auto& rt = this->responses_.at( ml );

  // Get the vectors of grid points
  const auto& wvec = rt.w_grid();
  const auto& qvec = rt.q_grid();

  // TODO: reduce code duplication with this and diff_xsec
  // Look up the masses of the projectile and ejectile
  int pdg_c = marley::Reaction::get_ejectile_pdg( pdg_a, proc_type_ );
  const auto& mt = marley::MassTable::Instance();
  double ma = mt.get_particle_mass( pdg_a );
  double mc = mt.get_particle_mass( pdg_c );

  // Projectile total energy
  double Ea = KEa + ma;

  // Get the grid step sizes
  // NOTE: this assumes that the grid is regularly spaced (I take advantage
  // of this to simplify the trapezoidal rule for integration)
  // TODO: revisit this assumption
  size_t num_w = wvec.size();
  size_t num_q = qvec.size();

  size_t num_w_minus_one = num_w - 1;
  size_t num_q_minus_one = num_q - 1;

  double dw = ( wvec.back() - wvec.front() ) / num_w_minus_one;
  double dq = ( qvec.back() - qvec.front() ) / num_q_minus_one;

  // Loop over every grid point. Do an angular integral at every omega
  // point, and use the results to integrate over omega.
  double integ = 0.;
  for ( size_t iw = 0u; iw < num_w; ++iw ) {

    // Get the energy transfer at the current grid point
    double w = wvec.at( iw );

    // Use it to compute the ejectile total energy, etc.
    double Ec = Ea - w;
    // We can skip unphysical terms for which the total energy is smaller than
    // the final lepton mass
    if ( Ec < mc ) continue;
    double pc = marley_utils::real_sqrt( Ec*Ec - mc*mc );
    double pa = marley_utils::real_sqrt( Ea*Ea - ma*ma );

    double w_integ = 0.;
    for ( size_t iq = 0u; iq < num_q; ++iq ) {

      // Get the magnitude of the 3-momentum transfer at the current grid point
      double q = qvec.at( iq );

      // Convert it into a value for the scattering cosine
      double cos_theta = ( pa*pa + pc*pc - q*q ) / ( 2.*pa*pc );

      // Compute the differential cross section at this 2D grid point
      double diff = this->diff_xsec( pdg_a, KEa, w, cos_theta, ml );

      // If it is larger than any value encountered so far, record it as
      // the new maximum
      if ( diff > diff_max ) diff_max = diff;

      // Apply a Jacobian to transform from dwdcostheta to dwdq (since we're
      // integrating on a grid of q values)
      diff *= q / pa / pc;

      // Add this term to the integral over omega according to the trapezoid
      // rule
      if ( iq == 0u || iq == num_q_minus_one ) diff /= 2.;
      w_integ += diff;
    }

    // We're done summing over q values. Scale the result to get the
    // integral over q for this w grid point.
    w_integ *= dq;

    // Now add this term to the integral over w according to the trapezoid rule
    if ( iw == 0u || iw == num_w_minus_one ) w_integ /= 2.;

    integ += w_integ;
  }

  // We're done. Scale the integral over w by the needed prefactor
  integ *= dw;

  return integ;
}