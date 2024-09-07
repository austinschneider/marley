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

// Standard library includes
#include <fstream>

// MARLEY includes
#include "marley/Error.hh"
#include "marley/FileManager.hh"
#include "marley/MassTable.hh"
#include "marley/Reaction.hh"
#include "marley/TabulatedXSec.hh"
#include "marley/marley_utils.hh"

marley::TabulatedXSec::TabulatedXSec( int target_pdg,
  marley::Reaction::ProcessType p_type, CoulombCorrector::CoulombMode mode,
  double Delta ) : ta_( target_pdg ), proc_type_( p_type ),
  coulomb_mode_( mode ), Delta_( Delta )
{
  // Set the residue PDG code (pdg_d_) based on the input information
  int dummy_charge;
  marley::Reaction::get_residue_pdg_and_charge( p_type, target_pdg,
    pdg_d_, dummy_charge );
}

void marley::TabulatedXSec::add_table( const std::string& file_name )
{
  const auto& fm = marley::FileManager::Instance();
  std::string full_file_name = fm.find_file( file_name );
  if ( full_file_name.empty() ) {
    throw marley::Error( "Could not open the nuclear response data file \""
      + file_name + '\"' );
  }

  std::ifstream in_file( full_file_name );

  // TODO: add error handling for file parsing

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
    std::vector<marley::NuclearResponses> >();

  auto unnat_resp = std::make_shared<
    std::vector<marley::NuclearResponses> >();

  // Temporary storage for responses
  double rcc, rll, rcl, rtVV, rtAA, rtprime;

  // Flag used to swap back and forth between natural and unnatural parity
  bool nat = true;
  while ( in_file >> rcc >> rll >> rcl >> rtVV >> rtAA >> rtprime ) {

    if ( nat ) nat_resp->emplace_back( rcc, rll, rcl, rtVV, rtAA, rtprime );
    else unnat_resp->emplace_back( rcc, rll, rcl, rtVV, rtAA, rtprime );

    // Flip the parity flag
    nat = !nat;
  }

  // Store the finished tables as new objects in the map
  responses_[ nat_ml ] = ResponseTable( wvec, qvec, nat_resp );
  responses_[ unnat_ml ] = ResponseTable( wvec, qvec, unnat_resp );

}

double marley::TabulatedXSec::diff_xsec( int pdg_a, double KEa, double omega,
  double cos_theta, const marley::TabulatedXSec::MultipoleLabel& ml ) const
{
  int helicity = marley_utils::get_particle_helicity( pdg_a );

  // Look up the masses of the projectile and ejectile
  int pdg_c = marley::Reaction::get_ejectile_pdg( pdg_a, proc_type_ );
  const auto& mt = marley::MassTable::Instance();
  double ma = mt.get_particle_mass( pdg_a );
  double mc = mt.get_particle_mass( pdg_c );

  // Determine the total energies and momenta of the initial and final leptons
  double Ea = KEa + ma;
  double pa = marley_utils::real_sqrt( Ea*Ea - ma*ma );
  double Ec = Ea - omega + this->Delta();
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
  double sin_theta2 = 1. - cos_theta*cos_theta;
  double q2 = q*q;
  double vcc = 1. + beta * cos_theta;
  double vll = vcc - 2.*Ea*Ec*sin_theta2*beta*beta/q2;
  double vcl = -2. * ( omega*vcc/q + mc*mc/Ec/q );
  double vT = 1. - beta*cos_theta + Ea*Ec*beta*beta*sin_theta2/q2;
  double vTprime = 2. * helicity * ( (Ea + Ec)*(1 - beta*cos_theta)/q
    - mc*mc/q/Ec );
  LeptonFactors lf( vcc, vll, vcl, vT, vTprime );

  // Compute the double-differential cross section with respect to
  // the energy transfer and the scattering cosine
  double xsec = lf * nr;
  xsec *= 2. * marley_utils::GF2 * marley_utils::Vud2 * Ec * pc;

  // For CC cross sections, apply the appropriate correction factor for
  // Coulomb corrections based on the active "Coulomb mode"
  bool is_charged_current = this->is_cc();

  if ( is_charged_current ) {

    // Get the mass of the atomic target
    double mb = mt.get_atomic_mass( ta_.pdg() );

    // Lab-frame total energy of the residue (via energy conservation)
    double Ed = omega + mb;

    // Mass of the residue (via 3-momentum conservation)
    double md = marley_utils::real_sqrt( Ed*Ed - q*q );

    // Dot product of the 4-momenta of the ejectile and residue
    double pc_dot_pd = Ec*Ed + pc*pc - pa*pc*cos_theta;

    // Manifestly Lorentz-invariant relative speed of particles c and d
    double beta_rel_cd = marley_utils::real_sqrt(
    std::pow(pc_dot_pd, 2) - mc*mc*md*md ) / pc_dot_pd;

    // Compute the Coulomb correction factor using the relative speed
    CoulombCorrector coul_corr( pdg_c, pdg_d_, coulomb_mode_ );
    double FC = coul_corr.coulomb_correction_factor( beta_rel_cd );

    // Apply it to the differential cross section
    xsec *= FC;
  }

  return xsec;
}

double marley::TabulatedXSec::compute_integral( int pdg_a, double KEa,
  const marley::TabulatedXSec::MultipoleLabel& ml, double& diff_max,
  std::vector< IntegralTerm >* integral_terms ) const
{
  // Set the maximum differential cross section to zero to start
  diff_max = 0.;

  // If the input vector is not null, then clear it in preparation for
  // accumulating the terms of the integral
  if ( integral_terms ) integral_terms->clear();

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
    double w_eff = w - Delta_;

    // Use it to compute the ejectile total energy, etc.
    double Ec = Ea - w_eff;
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

      // If the user supplied a vector to record the terms in the integral,
      // then store the contribution from the current iteration before
      // continuing to the next grid point
      if ( integral_terms ) {
        double value = diff * dq * dw;
        if ( iw == 0u || iw == num_w_minus_one ) value /= 2.;
        integral_terms->emplace_back( w, cos_theta, value );
      }

    } // q grid points

    // We're done summing over q values. Scale the result to get the
    // integral over q for this w grid point.
    w_integ *= dq;

    // Now add this term to the integral over w according to the trapezoid rule
    if ( iw == 0u || iw == num_w_minus_one ) w_integ /= 2.;

    integ += w_integ;

  } // w grid points

  // We're done. Scale the integral over w by the needed prefactor
  integ *= dw;

  return integ;
}

double marley::TabulatedXSec::integral( int pdg_a, double KEa,
  const marley::TabulatedXSec::MultipoleLabel& ml, double& diff_max ) const
{
  // First attempt to look up a pair of ChebyshevInterpolatingFunction
  // objects for the given projectile PDG code and multipole
  OptimizationMapKey key( pdg_a, ml );
  auto end = optimization_map_.end();
  auto iter = optimization_map_.find( key );
  // If we found one, attempt to compute the cross sections using them
  if ( iter != end ) {
    const auto& omv = iter->second;
    const auto& tot = omv.tot_xsec_;
    // If the requested kinetic energy is below threshold, then just
    // return zero
    if ( KEa < tot.x_min() ) return 0.;
    // If the kinetic energy is within the range covered by the
    // interpolating functions, then go ahead and use them
    else if ( KEa <= tot.x_max() ) {
      double tot_xsec = omv.tot_xsec_.evaluate( KEa );
      diff_max = omv.max_diff_xsec_.evaluate( KEa );

      // Ensure that interpolation problems will not give us a negative value
      // for the total or maximum differential cross sections
      tot_xsec = std::max( 0., tot_xsec );
      diff_max = std::max( 0., diff_max );

      return tot_xsec;
    }
  }

  // If we couldn't find a suitable interpolating function, then fall back to
  // brute-force integration
  return this->compute_integral( pdg_a, KEa, ml, diff_max );
}

double marley::TabulatedXSec::integral( int pdg_a, double KEa ) const
{
  double integ = 0.;
  for ( const auto& pair : responses_ ) {
    double dummy;
    const auto& ml = pair.first;
    double integ_ml = this->integral( pdg_a, KEa, ml, dummy );
    integ += integ_ml;
  }
  return integ;
}

void marley::TabulatedXSec::optimize( int pdg_a, double max_KEa ) {
  // Loop over each of the multipoles
  for ( const auto& pair : responses_ ) {
    const auto& ml = pair.first;

    double min_KEa = 0.;
    std::function<double(double)> tot_xsec_func = [](double)
      -> double { return 0.; };
    std::function<double(double)> max_diff_xsec_func = tot_xsec_func;

    // Verify that the total cross section for this multipole is non-vanishing
    // for the requested maximum projectile kinetic energy KEa. If it vanishes,
    // then just skip the current multipole.
    double dummy;
    double xsec_at_max = this->integral( pdg_a, max_KEa, ml, dummy );

    if ( xsec_at_max > 0. ) {
      // Do a binary search to find the threshold for this multipole. Continue
      // until we've found it within the given tolerance.
      constexpr double thresh_tol = 1e-6; // MeV
      // Set up the bounds of a bracketing interval that will contain the
      // kinetic energy threshold
      double low_KEa = 0.;
      double high_KEa = max_KEa;
      do {
        // Check the total cross section at the midpoint of the current
        // bracketing interval
        double cur_KEa = (low_KEa + high_KEa) / 2.;
        double xsec = this->integral( pdg_a, cur_KEa, ml, dummy );
        // If it vanishes, move the lower bound up
        if ( xsec <= 0. ) low_KEa = cur_KEa;
        // If it doesn't, move the upper bound down
        else high_KEa = cur_KEa;
        // Continue until the bracketing interval is no larger than the
        // tolerance defined above
      } while ( std::abs(high_KEa - low_KEa) > thresh_tol );

      // Adopt the lower bound of the bracketing interval as the threshold
      min_KEa = low_KEa;

      tot_xsec_func = [&, this](double KEa)
        -> double { return this->integral( pdg_a, KEa, ml, dummy ); };

      max_diff_xsec_func = [&, this](double KEa) -> double {
        double max_diff;
        this->integral( pdg_a, KEa, ml, max_diff );
        return max_diff;
      };
    }

    MARLEY_LOG_INFO() << "Optimizing total cross section for "
      << ml.J_ << ml.Pi_;

    // Now we're ready to build the Chebyshev interpolating functions for
    // this multipole. We need one for the total cross section and the
    // other one for the maximum of the differential cross section.
    ChebyshevInterpolatingFunction tot_xs_cif( tot_xsec_func, min_KEa,
      max_KEa, 64 );
    // TODO: do you want adaptive grid sizing here?

    MARLEY_LOG_DEBUG() << "Optimizing max diff for "
      << ml.J_ << ml.Pi_;

    ChebyshevInterpolatingFunction max_diff_xs_cif( max_diff_xsec_func,
      min_KEa, max_KEa, 64 );

    // Build the key and value we need for the map entry
    OptimizationMapKey key( pdg_a, ml );
    OptimizationMapValue value( tot_xs_cif, max_diff_xs_cif );

    optimization_map_.emplace( std::make_pair(key, value) );
  }
}

// For CC interactions, apply an energy shift to account for the difference in
// energy between the ground state of the initial nucleus and the isobaric
// analog state in the final nucleus. This leads to an effective energy
// transfer value omega_eff which is used when evaluating the final-state
// lepton energy.
double marley::TabulatedXSec::Delta() const {
  bool is_charged_current = this->is_cc();

  double result = 0.;
  if ( is_charged_current ) result = Delta_;
  return result;
}
