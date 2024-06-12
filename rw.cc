// Standard library includes
#include <algorithm>
#include <cmath>

// HepMC3 includes
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

// MARLEY includes
#include "marley/hepmc3_utils.hh"
#include "marley/marley_utils.hh"
#include "marley/Error.hh"
#include "marley/EventFileReader.hh"
#include "marley/HauserFeshbachDecay.hh"
#include "marley/JSON.hh"
#include "marley/Logger.hh"
#include "marley/StructureDatabase.hh"

constexpr double PRETTY_SMALL = 1e-6;

int main( int argc, char* argv[] ) {

  // If the user has not supplied enough command-line arguments, display the
  // standard help message and exit
  if ( argc <= 1 ) {
    std::cout << "Usage: " << argv[0] << " INPUT_FILE...\n";
    return 0;
  }

  //// Check whether the output file exists and warn the user before
  //// overwriting it if it does
  //std::ifstream temp_stream( argv[1] );
  //if ( temp_stream ) {
  //  bool overwrite = marley_utils::prompt_yes_no(
  //    "Really overwrite " + std::string(argv[1]) + '?');
  //  if ( !overwrite ) {
  //    std::cout << "Action aborted.\n";
  //    return 0;
  //  }
  //}

  // Create an alternative structure database with a custom set of optical
  // model parameters
  const std::string rw_om_config_file( "my_om_rw_config.js" );
  auto rw_config = marley::JSON::load_file( rw_om_config_file );
  marley::StructureDatabase sdb_alt;
  sdb_alt.load_optical_model_params( &rw_config );

  // Prepare to read the input file(s)
  std::vector< std::string > input_file_names;
  for ( int i = 1; i < argc; ++i ) input_file_names.push_back( argv[i] );

  // File loop
  for ( const auto& file_name : input_file_names ) {

    // Open the current file for reading
    marley::EventFileReader efr( file_name );
    std::cout << "Opened file \"" << file_name << "\"\n";

    // Temporary object to use for reading in saved events
    HepMC3::GenEvent ev;

    // Event loop
    int event_num = 0;
    while ( efr >> ev ) {

      //if ( event_num % 1000 == 0 )
      std::cout << "Event " << event_num << '\n';

      double weight = 1.;

      // Get a vector of pointers to all decay vertices in the event that
      // were handled using the Hauser-Feshbach model
      auto hf_vtx_vec = marley_hepmc3::get_vertices_with_status(
        marley_hepmc3::NUHEPMC_HF_DECAY_VERTEX, ev );

      // Include each vertex that was found in the reweighting calculation
      for ( const auto& vtx : hf_vtx_vec ) {

        // Start with the default assumption that the decay was to a discrete
        // nuclear level
        bool decayed_to_continuum = false;

        const auto& in_vec = vtx->particles_in();
        const auto& out_vec = vtx->particles_out();

        // Double-check that the expected numbers of particles are present
        if ( in_vec.size() != 1u || out_vec.size() != 2u ) {
          throw marley::Error( "Hauser-Feshbach vertex encountered that"
            " is not a binary decay" );
        }

        // Retrieve the decay widths computed in the original simulation
        auto width_tot_attr = vtx->attribute< HepMC3::DoubleAttribute >(
          "TotalWidth" );
        double width_tot = width_tot_attr->value();

        auto width_ec_attr = vtx->attribute< HepMC3::DoubleAttribute >(
          "ECWidth" );
        double width_ec = width_ec_attr->value();

        // This attribute is only present for decays to the continuum
        auto width_sp_attr = vtx->attribute< HepMC3::DoubleAttribute >(
          "SPWidth" );
        double width_sp = 0.;
        if ( width_sp_attr ) {
          decayed_to_continuum = true;
          width_sp = width_sp_attr->value();
        }

        // Get access to the mother nucleus for the current binary decay
        const auto& mother = in_vec.front();

        // Retrieve the starting nuclear excitation energy, spin, and parity
        // TODO: add error handling here for missing attributes
        auto Exi_attr = mother->attribute< HepMC3::DoubleAttribute >( "Ex" );
        double Exi = Exi_attr->value();

        auto twoJi_attr = mother->attribute< HepMC3::IntAttribute >( "twoJ" );
        int twoJi = twoJi_attr->value();

        auto Pi_attr = mother->attribute< HepMC3::IntAttribute >( "parity" );
        marley::Parity Pi( Pi_attr->value() );

        // NucleusDecayer populates the outgoing particles in the decay
        // vertex in a specific order. The first is the emitted nuclear
        // fragment or gamma-ray, while the second is the daughter nucleus.
        const auto& emitted_particle = out_vec.front();
        const auto& daughter = out_vec.back();

        int emitted_pdg = emitted_particle->pid();

        // Retrieve the final nuclear excitation energy, spin, and parity
        // TODO: add error handling here for missing attributes
        auto Exf_attr = daughter->attribute< HepMC3::DoubleAttribute >( "Ex" );
        double Exf = Exf_attr->value();

        auto twoJf_attr = daughter->attribute< HepMC3::IntAttribute >( "twoJ" );
        int twoJf = twoJf_attr->value();

        auto Pf_attr = daughter->attribute< HepMC3::IntAttribute >( "parity" );
        marley::Parity Pf( Pf_attr->value() );

        // Construct a set of decay widths calculated with alternative settings
        marley::HauserFeshbachDecay hf_alt( mother, Exi, twoJi, Pi, sdb_alt );

        // Get the total decay width under the alternative calculation
        double width_tot_alt = hf_alt.total_width();

        // Find the ExitChannel in the alternative calculation corresponding
        // to the original decay that was sampled
        const auto& ec_vec = hf_alt.exit_channels();
        auto ec_iter = std::find_if( ec_vec.cbegin(), ec_vec.cend(),
          [ emitted_pdg, decayed_to_continuum, Exf ](
            const std::unique_ptr< marley::ExitChannel >& test_ec ) -> bool
          {
            if ( emitted_pdg != test_ec->emitted_particle_pdg() ) return false;

            if ( decayed_to_continuum ) {
              if ( !test_ec->is_continuum() ) return false;
              else return true;
            }

            // If we get to here, then we're dealing with a discrete transition
            if ( test_ec->is_continuum() ) return false;

            const auto* dec = dynamic_cast<
              const marley::DiscreteExitChannel* >( test_ec.get() );
            if ( !dec ) return false;

            // Check the excitation energy of the final discrete level as
            // a last confirmation that we've found the correct ExitChannel
            const auto& lev = dec->get_final_level();
            double Ex_level = lev.energy();

            // To deal with possible numerical precision issues, this check of
            // the final excitation energy allows for small differences
            if ( std::abs(Exf - Ex_level) > PRETTY_SMALL ) return false;
            return true;
          }
        );

        if ( ec_iter == ec_vec.cend() ) {
          std::cout << "WARNING: could not find ExitChannel.\n";
          weight = 0.;
          continue;
        }

        // Partial width for the chosen ExitChannel under the alternative
        // calculation
        double width_ec_alt = ( *ec_iter )->width();

        // Partial differential width for the chosen spin-parity under the
        // alternative calculation (applies only to decays to the continuum)
        double width_sp_alt = 0.;
        if ( decayed_to_continuum ) {
          // Evaluate the differential width at the sampled excitation energy.
          // This will populate the table of SpinParityWidth objects with the
          // correct values.
          const auto& cec = dynamic_cast<
            const marley::ContinuumExitChannel& >( *(*ec_iter) );
          cec.differential_width( Exf, true );

          // Retrieve extra information based on the kind of emitted particle
          int mpol, two_j_frag, orb_l;
          bool emitted_gamma = false;
          auto mpol_attr = vtx->attribute< HepMC3::IntAttribute >(
            "multipolarity" );
          if ( mpol_attr ) {
            emitted_gamma = true;
            mpol = mpol_attr->value();
          }
          else {
            auto two_j_frag_attr = vtx->attribute< HepMC3::IntAttribute >(
              "two_j_frag" );
            auto orb_l_attr = vtx->attribute< HepMC3::IntAttribute >( "orb_l" );

            // TODO: add error handling for missing attributes here

            two_j_frag = two_j_frag_attr->value();
            orb_l = orb_l_attr->value();
          }

          double width_sp = 0.;
          if ( width_sp_attr ) {
            decayed_to_continuum = true;
            width_sp = width_sp_attr->value();
          }

          // Find the SpinParityWidth object corresponding to the spin-parity
          // value that was actually sampled
          const auto& spw_vec = cec.get_spw_table();
          auto spw_iter = std::find_if( spw_vec.cbegin(), spw_vec.cend(),
            [ twoJf, Pf, emitted_gamma, mpol, two_j_frag, orb_l ](
              const std::unique_ptr< marley::ContinuumExitChannel
                ::SpinParityWidth >& spw )
            -> bool
            {
              if ( Pf != spw->Pf ) return false;
              if ( twoJf != spw->twoJf ) return false;
              if ( emitted_gamma ) {
                const auto* g_spw = static_cast< const marley
                  ::GammaContinuumExitChannel::GammaSpinParityWidth* >(
                  spw.get() );
                if ( !g_spw ) return false;
                if ( mpol != g_spw->multipolarity ) return false;
              }
              else {
                // emitted fragment
                const auto* f_spw = static_cast< const marley
                  ::FragmentContinuumExitChannel::FragmentSpinParityWidth* >(
                  spw.get() );
                if ( !f_spw ) return false;
                if ( two_j_frag != f_spw->two_j_frag ) return false;
                if ( orb_l != f_spw->orb_l ) return false;
              }
              return true;
            }
          );

          if ( !emitted_gamma ) {

            std::cout << "VERTEX HAD twoJf = " << twoJf << ", Pf = "
                << Pf << ", two_j_frag = " << two_j_frag
                << ", orb_l = " << orb_l << '\n';

            std::cout << "*****Starting list of angular momenta*****\n";
            for ( const auto& spw : spw_vec ) {
              const auto* f_spw = static_cast< const marley
                ::FragmentContinuumExitChannel
                ::FragmentSpinParityWidth* >( spw.get() );

              std::cout << "DEBUG! twoJf = " << f_spw->twoJf << ", Pf = "
                << f_spw->Pf << ", two_j_frag = " << f_spw->two_j_frag
                << ", orb_l = " << f_spw->orb_l << ", diff_width = "
                << f_spw->diff_width << '\n';
            }
            std::cout << "*****Ending list of angular momenta*****\n";
          }

          if ( spw_iter == spw_vec.cend() ) {
            std::cout << "WARNING: could not find SpinParityWidth.\n";
            weight = 0.;
            continue;
          }

          // Store the partial differential width to this spin-parity state
          // under the alternative calculation
          width_sp_alt = spw_iter->get()->diff_width;
      }

        // We've accumulated all the width values we need to compute an
        // event weight. Do some sanity checks beforehand.
        std::string bad_width_name;
        if ( width_tot_alt <= 0. ) {
          bad_width_name = "Alternate total_width";
        }
        else if ( width_tot <= 0. ) {
          bad_width_name = "Original total width";
        }
        else if ( width_ec <= 0. ) {
          bad_width_name = "Original exit channel";
        }
        else if ( width_ec_alt <= 0. ) {
          bad_width_name = "Alternate exit channel";
        }
        else if ( decayed_to_continuum ) {
          if ( width_sp_alt <= 0. ) bad_width_name = "Alternate spin-parity";
          else if ( width_sp <= 0. ) bad_width_name = "Original spin-parity";
        }

        if ( !bad_width_name.empty() ) {
          std::cout << "WARNING: " << bad_width_name << " is non-positive.\n";
          weight = 0.;
          continue;
        }

        // Everything looks okay, so do the weight calculation using the
        // appropriate expression for a decay to the continuum or a discrete
        // nuclear level
        double w = width_tot / width_tot_alt;
        if ( decayed_to_continuum ) {
          w *= width_sp_alt / width_sp;
        }
        else {
          w *= width_ec_alt / width_ec;
        }

        // Multiply the weight for the current decay vertex into the overall
        // event weight
        weight *= w;

      } // decay vertex loop

      std::cout << "  weight = " << weight << '\n';

      ++event_num;

    } // event loop
  } // file loop

  return 0;
}
