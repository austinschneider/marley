#pragma once
#include "LevelDensityModel.hh"

namespace marley {

  /// @brief Implementation of the back-shifted Fermi gas nuclear level
  /// density model
  /// @details %Level density parameters used in this class are based on the
  /// global fits originally performed in <a
  /// href="http://dx.doi.org/10.1016/j.nuclphysa.2008.06.005">A. J. Koning, et
  /// al., Nucl. Phys. A810 (2008) pp. 13-76</a>. In addition to being used in
  /// MARLEY, this model has been adopted as one of the available options in the
  /// <a href="http://www-nds.iaea.org/RIPL-3/">RIPL-3</a> data library
  /// and the <a href="http://talys.eu">TALYS</a> nuclear code.
  class BackshiftedFermiGasModel : public LevelDensityModel {

    public:

      /// Create a back-shifted Fermi gas model object for a specific nuclide.
      /// This constructor will use global fits to initialize all level density
      /// parameters.
      /// @param Z atomic number of the desired nuclide
      /// @param A mass number of the desired nuclide
      /// @todo Add other constructors to the BackshiftedFermiGasModel class to
      /// allow the user to set the level density parameters based on local
      /// fits
      BackshiftedFermiGasModel(int Z, int A);

      /// @copydoc marley::LevelDensityModel::level_density(double Ex)
      virtual double level_density(double Ex) override;

      /// @copydoc marley::LevelDensityModel::level_density(double Ex,
      ///   int two_J)
      virtual double level_density(double Ex, int two_J) override;

      /// @copydoc marley::LevelDensityModel::level_density(double Ex,
      ///   int two_j, marley::Parity Pi)
      /// @details The current implementation assumes parity equipartition.
      virtual double level_density(double Ex, int two_J, marley::Parity Pi)
        override;

    protected:

      int Z_; ///< atomic number for this nuclide
      int A_; ///< mass number for this nuclide

      double sigma_; ///< spin cut-off parameter

      double a_tilde_; ///< asymptotic level density parameter (MeV<sup> -1</sup>)
      double gamma_; ///< damping parameter (MeV<sup> -1</sup>)
      double delta_W_; ///< shell correction energy (MeV)
      double Delta_BFM_; ///< excitation energy shift (MeV)

      /// @brief global fit for discrete-region spin cut-off parameter
      double sigma_d_global_;
      double Sn_; ///< neutron separation energy (MeV)
  };
}