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

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>

#include "marley/marley_utils.hh"
#include "marley/AllowedNuclearReactionWithQ2.hh"
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/LeptonFactors.hh"
#include "marley/Level.hh"
#include "marley/Logger.hh"
#include "marley/MatrixElement.hh"
//#include "marley/NucleusDecayer.hh"
#include "marley/NuclearResponses.hh"
#include "marley/FormFactor.hh"
#include "marley/Integrator.hh"

using ME_Type = marley::MatrixElement::TransitionType;
using ProcType = marley::Reaction::ProcessType;

namespace {
  constexpr int BOGUS_TWO_J_VALUE = -99999;
}

marley::AllowedNuclearReactionWithQ2::AllowedNuclearReactionWithQ2( ProcType pt,
  int pdg_a, int pdg_b, int pdg_c, int pdg_d, int q_d,
  const std::shared_ptr< std::vector<marley::MatrixElement> >& mat_els,
  const std::pair< std::vector<int>, std::vector<double> > nucleon_radii,
  marley::CoulombCorrector::CoulombMode mode,
  marley::FormFactor::FFScalingMode ff_scaling_mode, bool superallowed )
  : marley::NuclearReaction( pt, pdg_a, pdg_b, pdg_c, pdg_d, q_d ),
  matrix_elements_( mat_els ), nucleon_radii_( nucleon_radii ),
  coulomb_corrector_( pdg_c, pdg_d, mode ),
  form_factor_( ff_scaling_mode ), superallowed_( superallowed )
{
}

// Creates an event object by sampling the appropriate quantities and
// performing kinematic calculations
std::shared_ptr< HepMC3::GenEvent > marley::AllowedNuclearReactionWithQ2
  ::create_event( int pdg_a, double KEa, marley::Generator& gen ) const
{
  // Check that the projectile supplied to this event is correct. If not, alert
  // the user that this event does not use the requested projectile.
  if ( pdg_a != pdg_a_ ) throw marley::Error( "Could not create this event."
    " The requested projectile particle ID, " + std::to_string( pdg_a )
    + ", does not match the projectile particle ID, " + std::to_string( pdg_a_ )
    + ", in the reaction dataset." );

  // Sample a final residue energy level. First, check to make sure the given
  // projectile energy is above threshold for this reaction.
  if ( KEa < KEa_threshold_ ) throw std::range_error( "Could not create"
    " this event. Projectile kinetic energy " + std::to_string( KEa )
    + " MeV is below the threshold value " + std::to_string( KEa_threshold_ )
    + " MeV." );

  /// @todo Add more error checks to AllowedNuclearReaction::create_event() as
  /// necessary

  // Create an empty vector of sampling weights (partial total cross
  // sections to each kinematically accessible final level)
  std::vector<double> level_weights;

  // Create a discrete distribution object for level sampling.
  // Its default constructor creates a single weight of 1.
  // We will always explicitly give it weights to use when sampling
  // levels, so we won't worry about its default behavior.
  static std::discrete_distribution< size_t > ldist;

  // Compute the total cross section for a transition to each individual
  // nuclear level, and save the results in the level_weights vector (which
  // will be cleared by summed_xs_helper() before being loaded with the cross
  // sections). The summed_xs_helper() method can also be used for differential
  // (d\sigma/d\cos\theta_c^{CM}) cross sections, so supply a dummy
  // cos_theta_c_cm value and request total cross sections by setting the last
  // argument to false.
  double dummy = 0.;
  double sum_of_xsecs = summed_xs_helper( pdg_a, KEa, dummy,
    &level_weights, false );

  // Note that the elements in matrix_elements_ are given in order of
  // increasing excitation energy (this is currently enforced by the reaction
  // data format and is checked during parsing). This ensures that we can
  // sample a matrix element index from level_weights (which is populated in
  // the same order by summed_xs_helper()) and have it refer to the correct
  // object.

  // Complain if none of the levels we have data for are kinematically allowed
  if ( level_weights.empty() ) {
    throw marley::Error( "Could not create this event. The DecayScheme object"
      " associated with this reaction does not contain data for any"
      " kinematically accessible levels for a projectile kinetic energy of "
      + std::to_string( KEa ) + " MeV (max E_level = "
      + std::to_string( max_level_energy(KEa) ) + " MeV)." );
  }

  // Complain if the total cross section (the sum of all partial level cross
  // sections) is zero or negative (the latter is just to cover all
  // possibilities).
  if ( sum_of_xsecs <= 0. ) {
    throw marley::Error( "Could not create this event. All kinematically"
      " accessible levels for a projectile kinetic energy of "
      + std::to_string( KEa ) + " MeV (max E_level = "
      + std::to_string( max_level_energy(KEa) )
      + " MeV) have vanishing matrix elements." );
  }

  // Create a list of parameters used to supply the weights to our discrete
  // level sampling distribution
  std::discrete_distribution<size_t>::param_type params( level_weights.begin(),
    level_weights.end() );

  // Sample a matrix_element using our discrete distribution and the
  // current set of weights
  size_t me_index = gen.sample_from_distribution( ldist, params );

  const auto& sampled_matrix_el = matrix_elements_->at( me_index );

  // Get the energy of the selected level.
  double E_level = sampled_matrix_el.level_energy();

  // Update the residue mass based on its excitation energy for the current
  // event
  md_ = md_gs_ + E_level;

  // Compute Mandelstam s, the ejectile's CM frame total energy, the magnitude
  // of the ejectile's CM frame 3-momentum, and the residue's CM frame total
  // energy.
  double s, Ec_cm, pc_cm, Ed_cm;
  two_two_scatter( KEa, s, Ec_cm, pc_cm, Ed_cm );

  // Determine the CM frame velocity of the ejectile
  double beta_c_cm = pc_cm / Ec_cm;

  // Sample a CM frame scattering cosine for the ejectile.
  double cos_theta_c_cm = sample_cos_theta_c_cm( sampled_matrix_el, KEa,
    beta_c_cm, gen );

  // Sample a CM frame azimuthal scattering angle (phi) uniformly on [0, 2*pi).
  // We can do this because the matrix elements are azimuthally invariant
  double phi_c_cm = gen.uniform_random_double( 0.,
    marley_utils::two_pi, false );

  // Load the initial residue twoJ and parity values into twoJ and P. These
  // variables are included in the event record and used by NucleusDecayer
  // to start the Hauser-Feshbach decay cascade.
  int twoJ = BOGUS_TWO_J_VALUE;
  marley::Parity P; // defaults to positive parity

  // Get access to the nuclear structure database owned by the Generator
  auto& sdb = gen.get_structure_db();

  // Retrieve the ground-state spin-parity of the initial nucleus
  int twoJ_gs;
  marley::Parity P_gs;

  sdb.get_gs_spin_parity( pdg_b_, twoJ_gs, P_gs );

  // For transitions to discrete nuclear levels, all we need to do is retrieve
  // these values directly from the Level object
  const marley::Level* final_lev = sampled_matrix_el.level();
  if ( final_lev ) {
    twoJ = final_lev->twoJ();
    P = final_lev->parity();
  }
  // For transitions to the continuum, we rely on the spin-parity selection
  // rules to determine suitable values of twoJ and P. In cases where more than
  // one value is allowed, assume equipartition, and sample a spin-parity based
  // on the relative nuclear level densities at the excitation energy of
  // interest.
  /// @todo Revisit the equipartition assumption made here in favor of
  /// something better motivated.
  else {

    // For a Fermi transition, the final spin-parity is always the same as the
    // initial one
    if ( sampled_matrix_el.type() == ME_Type::FERMI ) {
      twoJ = twoJ_gs;
      P = P_gs;
    }
    else if ( sampled_matrix_el.type() == ME_Type::GAMOW_TELLER ) {
      // For a Gamow-Teller transition, the final parity is the same as the
      // initial parity
      P = P_gs;

      // For a spin-zero initial state, take a shortcut: the final spin will
      // always be one.
      if ( twoJ_gs == 0 ) twoJ = 2;
      else {

        // For an initial state with a non-zero spin, make a vector storing
        // all of the spin values allowed by the GT selection rules.
        // Sample an allowed value assuming equipartition of spin. Use the
        // relative final nuclear level densities as sampling weights.
        std::vector<int> allowed_twoJs;
        std::vector<double> ld_weights;

        auto& ldm = sdb.get_level_density_model( pdg_d_ );

        for ( int myTwoJ = std::abs(twoJ_gs - 2); myTwoJ <= twoJ_gs + 2;
          myTwoJ += 2 )
        {
          allowed_twoJs.push_back( myTwoJ );
          ld_weights.push_back( ldm.level_density(E_level, myTwoJ, P) );
        }

        std::discrete_distribution<size_t> my_twoJ_dist( ld_weights.begin(),
          ld_weights.end() );

        size_t my_index = gen.sample_from_distribution( my_twoJ_dist );
        twoJ = allowed_twoJs.at( my_index );
      }
    }
    else throw marley::Error( "Unrecognized matrix element type encountered"
      " in marley::AllowedNuclearReaction::create_event()" );
  }

  MARLEY_LOG_DEBUG() << "Sampled a " << sampled_matrix_el.type_str()
    << " transition from the " << marley::TargetAtom( pdg_b_ )
    << " ground state (with spin-parity " << static_cast<double>( twoJ_gs ) / 2.
    << P_gs << ") to the " << marley::TargetAtom( pdg_d_ )
    << " level with Ex = " << E_level << " MeV and spin-parity "
    << static_cast<double>( twoJ ) / 2. << P;

  // Create the preliminary event object (after 2-->2 scattering, but before
  // de-excitation of the residual nucleus)
  auto event = this->make_event_object( KEa, pc_cm, cos_theta_c_cm, phi_c_cm,
    Ec_cm, Ed_cm, E_level, twoJ, P );

  // Return the preliminary event object (to be processed later by the
  // NucleusDecayer class)
  return event;
}

// Compute the total reaction cross section (summed over all final nuclear
// levels) in units of MeV^(-2) using the center of momentum frame.
double marley::AllowedNuclearReactionWithQ2::total_xs( int pdg_a, double KEa )
  const
{
  double dummy_cos_theta = 0.;
  return summed_xs_helper( pdg_a, KEa, dummy_cos_theta, nullptr, false );
}

// Compute the differential cross section d\sigma / d\cos\theta_c^{CM}
// summed over all final nuclear levels. This is done in units of MeV^(-2)
// using the center of momentum frame.
double marley::AllowedNuclearReactionWithQ2::diff_xs( int pdg_a, double KEa,
  double cos_theta_c_cm ) const
{
  return summed_xs_helper( pdg_a, KEa, cos_theta_c_cm, nullptr, true );
}

// Compute the differential cross section d\sigma / d\cos\theta_c^{CM} for a
// transition to a particular final nuclear level. This is done in units of
// MeV^(-2) using the center of momentum frame.
double marley::AllowedNuclearReactionWithQ2::diff_xs(
  const marley::MatrixElement& mat_el, double KEa, double cos_theta_c_cm,
  double& beta_c_cm, bool check_max_E_level ) const
{
  // Check that the scattering cosine is within the physically meaningful range
  if ( std::abs(cos_theta_c_cm) > 1. ) return 0.;

  // Don't bother to compute anything if the matrix element vanishes
  if ( mat_el.strength() == 0. ) return 0.;

  // Also don't proceed further if the reaction is below threshold
  // (equivalently, if the requested level excitation energy E_level exceeds
  // that maximum kinematically-allowed value). To avoid redundant checks of
  // the threshold, skip this check if check_max_E_level is set to false.
  if ( check_max_E_level ) {
    double max_E_level = max_level_energy( KEa );
    if ( mat_el.level_energy() > max_E_level ) return 0.;
  }

  // The final nuclear mass (before nuclear de-excitations) is the sum of the
  // ground state residue mass plus the excitation energy of the accessed level
  // This includes a correction to the ground-state atomic mass to account for
  // production of an ion by charged-current interactions (see the constructor
  // of the NuclearReaction class)
  double md2 = std::pow( md_gs_ + mat_el.level_energy(), 2 );

  // Compute Mandelstam s (the square of the total CM frame energy)
  double s = std::pow( ma_ + mb_, 2 ) + 2.*mb_*KEa;
  double sqrt_s = std::sqrt( s );

  // Compute some CM frame total energies and 3-momentum magnitudes
  double Eb_cm = ( s + mb_*mb_ - ma_*ma_ ) / ( 2. * sqrt_s );

  double Ea_cm = sqrt_s - Eb_cm;
  double pa_cm = marley_utils::real_sqrt( std::pow(Ea_cm, 2) - ma_*ma_ );

  double Ec_cm = ( s + mc_*mc_ - md2 ) / ( 2. * sqrt_s );
  double pc_cm = marley_utils::real_sqrt( std::pow(Ec_cm, 2) - mc_*mc_ );

  // Compute the CM frame value of the energy transfer (omega)
  double omega_cm = Ea_cm - Ec_cm;

  // Compute the CM frame value of the 3-momentum transfer magnitude (kappa)
  double kappa_cm = marley_utils::real_sqrt( std::pow( pa_cm, 2 )
    + std::pow( pc_cm, 2 ) - 2. * pa_cm * pc_cm * cos_theta_c_cm );

  // Negative square of the 4-momentum transfer
  double Q2 = kappa_cm*kappa_cm - omega_cm*omega_cm;

  // Compute the (dimensionless) speed of the ejectile in the CM frame
  beta_c_cm = pc_cm / Ec_cm;

  // CM frame total energy of the nuclear residue
  double Ed_cm = sqrt_s - Ec_cm;

  // Common factors for both CC and NC differential cross sections to a discrete
  // nuclear level in the CM frame
  double diff_xsec_prefactor = ( marley_utils::GF2 / ( 2 * marley_utils::pi ) )
    * ( Eb_cm * Ed_cm / s ) * Ec_cm * pc_cm;

  // Apply extra factors based on the current process type
  if ( process_type_ == ProcessType::NeutrinoCC_Discrete
    || process_type_ == ProcessType::AntiNeutrinoCC_Discrete )
  {
    // Dot product of the four-momenta of particles c and d
    double pc_dot_pd = Ed_cm*Ec_cm + std::pow( pc_cm, 2 );

    // Relative speed of particles c and d, computed with a manifestly
    // Lorentz-invariant expression
    double beta_rel_cd = marley_utils::real_sqrt(
      std::pow(pc_dot_pd, 2) - mc_*mc_*md2 ) / pc_dot_pd;

    // Calculate a Coulomb correction factor using either a Fermi function
    // or the effective momentum approximation
    double factor_C = coulomb_corrector_.coulomb_correction_factor(
      beta_rel_cd );
    diff_xsec_prefactor *= marley_utils::Vud2 * factor_C;
  }
  else if ( process_type_ == ProcessType::NC_Discrete )
  {
    // For NC, extra factors are only needed for Fermi transitions (which
    // correspond to CEvNS since they can only access the nuclear ground state)
    if ( mat_el.type() == ME_Type::FERMI ) {
      double Q_w = weak_nuclear_charge();
      diff_xsec_prefactor *= 0.25*std::pow( Q_w, 2 );
    }
  }
  else throw marley::Error( "Unrecognized or invalid process type encountered"
    " in marley::AllowedNuclearReaction::diff_xs()" );

  // We're done with the overall factors. Now compute the lepton part of the
  // tensor contraction
  // @todo Reduce code duplication with a similar implementation of these (in
  // the lab frame) within the TabulatedXSec class. Note that there are
  // different conventions used there.
  int helicity = marley_utils::get_particle_helicity( pdg_a_ );
  double sin2_theta_c_cm = 1. - cos_theta_c_cm*cos_theta_c_cm;
  double vcc = 1. + beta_c_cm * cos_theta_c_cm;
  double vll = vcc - 2.*Ea_cm*Ec_cm*sin2_theta_c_cm
    * beta_c_cm*beta_c_cm/kappa_cm/kappa_cm;
  double vcl = -1.* ( omega_cm*vcc/kappa_cm + mc_*mc_/Ec_cm/kappa_cm );
  double vT = 1. - beta_c_cm*cos_theta_c_cm + Ea_cm*Ec_cm
    * beta_c_cm*beta_c_cm*sin2_theta_c_cm/kappa_cm/kappa_cm;
  double vTprime = helicity * ( (Ea_cm + Ec_cm)
    * (1. - beta_c_cm*cos_theta_c_cm)/kappa_cm - mc_*mc_/kappa_cm/Ec_cm );
  marley::LeptonFactors lf( vcc, vll, vcl, vT, vTprime );

  // Now compute the nuclear responses. Copy the energy transfer, 3-momentum
  // transfer magnitude, and Q^2 to effective values. These will be set to zero
  // if we're working in the allowed approximation.
  double omega_cm_eff = superallowed_ ? 0. : omega_cm;
  double kappa_cm_eff = superallowed_ ? 0. : kappa_cm;
  double Q2_eff = superallowed_ ? 0. : Q2;

  // Calculate the adjusted transition strength appropriate for the chosen
  // nuclear level
  double strength_eff = mat_el.strength();
  strength_eff *= this->bessel_factor( kappa_cm_eff );

  // Scale the strength by the relevant nucleon form factor and divide by the
  // relevant coupling constant. Use the scaled value to compute the nuclear
  // responses relevant to the chosen transition type
  double rCC, rCL, rLL, rTvv, rTaa, rTprime;
  const double kM = kappa_cm_eff / marley_utils::m_nucleon;
  const double k2M = kappa_cm_eff * kappa_cm_eff / marley_utils::m_nucleon;

  if ( mat_el.type() == ME_Type::FERMI ) {
    double F1 = form_factor_.F1( Q2_eff, marley_utils::M_V );
    strength_eff *= F1*F1 / marley_utils::g_V2;

    rCC = strength_eff;
    rCL = strength_eff * kM;
    rLL = strength_eff * kM * kM / 4.;
    rTvv = 0.;
    rTaa = 0.;
    rTprime = 0.;
  }
  else if ( mat_el.type() == ME_Type::GAMOW_TELLER ) {
    double FA = form_factor_.FA( Q2_eff, marley_utils::M_A );
    strength_eff *= FA*FA / marley_utils::g_A2;

    double FP = form_factor_.FP( Q2_eff, marley_utils::M_A );
    double F1 = form_factor_.F1( Q2_eff, marley_utils::M_V );
    double F2 = form_factor_.F2( Q2_eff, marley_utils::M_V );

    double FPA = FP / FA;
    double F12A = ( F1 + 2. * marley_utils::m_nucleon * F2 ) / FA;

    rCC = strength_eff * kM * kM / 12. * ( 1. - 2.*omega_cm_eff*FPA
      + omega_cm_eff*omega_cm_eff*FPA*FPA );
    rCL = strength_eff * kM / 3. * ( 1. - (omega_cm_eff + 0.5*k2M)*FPA
      + 0.5*omega_cm_eff*k2M*FPA*FPA );
    rLL = strength_eff * ( 1./3. - k2M/3.*FPA + k2M*k2M/12.*FPA*FPA );
    rTvv = 0.;
    rTaa = strength_eff * ( 2./3. + kM*kM/6.*F12A*F12A );
    rTprime = -1. * strength_eff * 2./3. * kM * F12A;
  }
  else throw marley::Error( "Unrecognized matrix element type encountered in"
    " marley::AllowedNuclearReaction::diff_xs()" );

  marley::NuclearResponses nr( rCC, rLL, rCL, rTvv, rTaa, rTprime );

  double diff_xsec = diff_xsec_prefactor * ( lf * nr );
  return diff_xsec;
}

// Compute the total reaction cross section (in MeV^(-2)) for a transition to a
// particular nuclear level using the center of momentum frame
double marley::AllowedNuclearReactionWithQ2::total_xs(
  const marley::MatrixElement& mat_el, double KEa, double& beta_c_cm,
  bool check_max_E_level ) const
{
  if ( superallowed_ ) return 2. * this->diff_xs( mat_el, KEa, 0., beta_c_cm,
    check_max_E_level );

  // Integrator object to integrate over the scattering angle
  static marley::Integrator integrator;

  // Don't bother to compute anything if the matrix element vanishes
  if ( mat_el.strength() == 0. ) return 0.;

  // Integrate the differential cross section over cos_theta_c_cm
  double total_xsec = integrator.num_integrate(
    [ &mat_el, KEa, &beta_c_cm, check_max_E_level, this ](
      double cos_theta_cm ) -> double
    { return this->diff_xs( mat_el, KEa, cos_theta_cm, beta_c_cm,
      check_max_E_level ); }, -1., 1.
  );

  MARLEY_LOG_DEBUG() << "total xsec " << description_
    << " to level with energy " << mat_el.level_energy() << " MeV is "
    << total_xsec << " MeV^(-2).";

  return total_xsec;
}

// Helper function for total_xs and diff_xs()
double marley::AllowedNuclearReactionWithQ2::summed_xs_helper( int pdg_a,
  double KEa, double cos_theta_c_cm, std::vector<double>* level_xsecs,
  bool differential ) const
{
  // Check that the projectile supplied to this event is correct. If not,
  // return a total cross section of zero since this reaction is not available
  // for the given projectile.
  /// @todo Consider whether you should use an exception if pdg_a != pdg_a_
  if ( pdg_a != pdg_a_ ) return 0.;

  // If we're evaluating a differential cross section, check that the
  // scattering cosine is within the physically meaningful range. If it's
  // not, then just return 0.
  if ( differential && std::abs(cos_theta_c_cm) > 1. ) return 0.;

  // If the projectile kinetic energy is zero (or negative), then
  // just return zero.
  if ( KEa <= 0. ) return 0.;

  // If we've been passed a vector to load with the partial cross sections
  // to each nuclear level, then clear it before storing them
  if ( level_xsecs ) level_xsecs->clear();

  double max_E_level = max_level_energy( KEa );
  double xsec = 0.;
  for ( const auto& mat_el : *matrix_elements_ ) {

    // Get the excitation energy for the current level
    double level_energy = mat_el.level_energy();

    // Exit the loop early if you reach a level with an energy that's too high
    if ( level_energy > max_E_level ) break;

    // Check whether the matrix element (B(F) + B(GT)) is nonvanishing for the
    // current level. If it is, just set the weight equal to zero rather than
    // computing the total xs.
    if ( mat_el.strength() != 0. ) {

      // Set the check_max_E_level flag to false when calculating the total
      // cross section for this level (we've already verified that the
      // current level is kinematically accessible in the check against
      // max_E_level above)
      double beta_c_cm = 0.;
      double partial_xsec;

      // Compute either the total or differential (d\sigma / d\cos\theta_{CM})
      // cross section as requested
      if ( differential ) {
        partial_xsec = diff_xs( mat_el, KEa, cos_theta_c_cm, beta_c_cm, false );
      } else {
        partial_xsec = total_xs( mat_el, KEa, beta_c_cm, false );
      }

      if ( std::isnan(partial_xsec) ) {
        MARLEY_LOG_WARNING() << "Partial cross section for reaction "
          << description_ << " gave NaN result.";
        MARLEY_LOG_DEBUG() << "Parameters were level energy = "
          << mat_el.level_energy() << " MeV, projectile kinetic energy = "
          << KEa << " MeV, and reduced matrix element = " << mat_el.strength()
          << ". Differential was set to " << differential << ".";
        MARLEY_LOG_DEBUG() << "The partial cross section to this level"
          << " will be set to zero.";
        partial_xsec = 0.;
      }

      xsec += partial_xsec;

      // Store the partial cross section to the current individual nuclear
      // level if needed (i.e., if level_xsecs is not nullptr)
      if ( level_xsecs ) level_xsecs->push_back( partial_xsec );
    }
  }

  return xsec;
}

// Bessel factor calculation
double marley::AllowedNuclearReactionWithQ2::bessel_factor( double kappa_cm )
  const
{
  // If kappa is zero, return 1. to avoid division by zero.
  if ( kappa_cm == 0. ) return 1.;

  // If the scaling mode is set to flat, return 1.
  if ( form_factor_.ff_scaling_mode()
    == marley::FormFactor::FFScalingMode::FLAT ) return 1.;

  // If the nucleon radii pair is empty, give a warning and return 1.
  if ( nucleon_radii_.first.empty() || nucleon_radii_.second.empty() ) {
    MARLEY_LOG_WARNING() << "Nuclear radii vector is empty in"
      << " marley::AllowedNuclearReaction::bessel_factor().";
    return 1.;
  }

  double bessel_sum = 0.;
  int nucleon_count = 0;
  int nucleon_limit;

  // Set the nucleon type depending on the process type
  if ( process_type_ == ProcessType::NeutrinoCC_Discrete ) {
    nucleon_limit = marley_utils::get_particle_Z( pdg_b_ );
  }
  else if ( process_type_ == ProcessType::AntiNeutrinoCC_Discrete ) {
    nucleon_limit = marley_utils::get_particle_A( pdg_b_ )
      - marley_utils::get_particle_Z( pdg_b_ );
  }
  else {
    // Throw an error
    /// @todo Fix for NC reactions
    throw marley::Error( "Unrecognized or invalid process type encountered in"
      " marley::AllowedNuclearReaction::bessel_factor()" );
  }

  for ( int i = 0; i < nucleon_radii_.first.size(); i ++ ) {
    int degeneracy = nucleon_radii_.first.at(i);
    double radius = nucleon_radii_.second.at(i);

    int level_i_count = (nucleon_count + degeneracy) <= nucleon_limit
      ? degeneracy : nucleon_limit - nucleon_count;
    nucleon_count += level_i_count;
    bessel_sum += level_i_count * std::sin( kappa_cm * radius )
      / ( kappa_cm * radius );
  }

  return std::pow( bessel_sum / nucleon_limit, 2 );
}


// Sample an ejectile scattering cosine in the CM frame.
double marley::AllowedNuclearReactionWithQ2::sample_cos_theta_c_cm(
  const marley::MatrixElement& mat_el, double KEa, double beta_c_cm,
  marley::Generator& gen ) const
{
  // For now the max is unknown, so the rejection sample will calculate it
  // "on the fly"
  double max = marley_utils::UNKNOWN_MAX;

  // For the superallowed case, we know where it will be a priori
  if ( superallowed_ ) {
    if ( mat_el.type() == ME_Type::FERMI ) {
      max = this->diff_xs( mat_el, KEa, 1., beta_c_cm, false );
    }
    else if ( mat_el.type() == ME_Type::GAMOW_TELLER ) {
      max = this->diff_xs( mat_el, KEa, -1., beta_c_cm, false );
    }
    else throw marley::Error( "Unrecognized matrix element type encountered"
      " in marley::AllowedNuclearReactionWithQ2::sample_cos_theta_c_cm()" );
  }

  return gen.rejection_sample(
    [ &mat_el, KEa, &beta_c_cm, this ]( double cos_theta_cm ) -> double
    { return this->diff_xs( mat_el, KEa, cos_theta_cm, beta_c_cm, false ); },
    -1., 1., max );
}

std::shared_ptr< HepMC3::GenEvent > marley::AllowedNuclearReactionWithQ2
  ::make_event_object( double KEa, double pc_cm, double cos_theta_c_cm,
  double phi_c_cm, double Ec_cm, double Ed_cm, double E_level, int twoJ,
  const marley::Parity& P ) const
{
  auto event = marley::Reaction::make_event_object( KEa, pc_cm,
    cos_theta_c_cm, phi_c_cm, Ec_cm, Ed_cm, E_level, twoJ, P );

  this->set_charge_attributes( event );

  return event;
}

// Adds an indication of whether the reaction populates excited levels of the
// daughter nucleus or only accesses the ground state (e.g., CEvNS)
void marley::AllowedNuclearReactionWithQ2::set_description() {
  // Initialize the description_ member variable with the basic information
  // first
  marley::NuclearReaction::set_description();
  // Now decide whether this reaction can access excited levels in the daughter
  // nucleus
  bool has_excited_state = false;
  for ( const auto& me : *matrix_elements_ ) {
    if ( me.level_energy() > 0. ) {
      has_excited_state = true;
      break;
    }
  }
  if ( has_excited_state ) description_ += '*';
  else description_ += " (g.s.)";
}
