#pragma once
#include <random>
#include "Gamma.hh"
#include "Parity.hh"

namespace marley {

  class Generator;
  
  /// A class that represents an ENSDF nuclear energy level
  
  class Level {
    public:
      /// @param jpi a string containing the spin
      /// and parity of the level (e.g., 0+)
      Level(double E, int twoJ, marley::Parity pi);
      marley::Gamma* add_gamma(const marley::Gamma& gamma);
      void clear_gammas();
      std::vector<marley::Gamma>* get_gammas();
      bool get_gamma_status() const;
      double get_energy() const;
      int get_two_J() const;
      marley::Parity get_parity() const;
      void set_energy(double E);
      void set_two_J(int twoJ);
      void set_parity(marley::Parity pi);
      const marley::Gamma* sample_gamma(marley::Generator& gen);
      std::string get_spin_parity_string() const;
  
    private:
      double energy; // MeV
      int two_J; // Two times the level spin (allows us to represent half-integer
                 // spins as integers)
      marley::Parity parity;
    
      bool gammas_known; // Determining whether or not the gammas are known
  
      std::vector<marley::Gamma> gammas;
      std::vector<double> gamma_intensities;
      std::discrete_distribution<int> gamma_dist;
  };

}