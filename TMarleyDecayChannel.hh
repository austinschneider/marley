#pragma once

#include "TMarleyFragment.hh"
#include "TMarleyGenerator.hh"
#include "TMarleyLevel.hh"
#include "TMarleyParity.hh"
#include "TMarleySphericalOpticalModel.hh"

class TMarleyDecayChannel {
  public:
    inline TMarleyDecayChannel(bool is_continuum): continuum(is_continuum)
    {
    }

    virtual int get_fragment_pid() = 0;
    virtual int get_fragment_two_s() = 0;
    virtual TMarleyParity get_fragment_parity() = 0;
    virtual int get_fragment_Z() = 0;
    virtual int get_fragment_A() = 0;
    virtual double get_fragment_mass() = 0;

    inline bool is_continuum() const {
      return continuum;
    }

    // Determine final values for several quantities of
    // interest for the residual nuclide after the decay occurs.
    // Input values are the initial state's excitation energy Ex,
    // 2*(total angular momentum) two_J, and parity Pi. These
    // variables are then loaded with the corresponding values
    // for the final state.
    virtual void get_post_decay_parameters(double& Ex, int& two_J,
      TMarleyParity& Pi) = 0;

    // Return the final excitation energy (at the center of the bin
    // if this is a continuum de-excitation channel)
    virtual double get_bin_center_Exf() = 0;

  protected:
    // Flag indicating whether this decay channel accesses the
    // unbound continuum (true) or a discrete level (false)
    bool continuum;
};

class TMarleyFragmentDecayChannel : public TMarleyDecayChannel {
  public:
    TMarleyFragmentDecayChannel(const TMarleyFragment& frag,
      bool is_continuum): TMarleyDecayChannel(is_continuum),
      fragment(frag)
    {
    }

    inline int get_fragment_pid() { return fragment.get_pid(); }
    inline int get_fragment_two_s() { return fragment.get_two_s(); }
    inline TMarleyParity get_fragment_parity() { return fragment.get_parity(); }
    inline int get_fragment_Z() { return fragment.get_Z(); }
    inline int get_fragment_A() { return fragment.get_A(); }
    inline double get_fragment_mass() { return fragment.get_mass(); }

    virtual void get_post_decay_parameters(double& Ex, int& two_J,
      TMarleyParity& Pi) = 0;

    // Return the final excitation energy (at the center of the bin
    // if this is a continuum de-excitation channel)
    virtual double get_bin_center_Exf() = 0;

  protected:
    const TMarleyFragment& fragment;
};

class TMarleyGammaDecayChannel : public TMarleyDecayChannel {
  public:
    TMarleyGammaDecayChannel(bool is_continuum):
      TMarleyDecayChannel(is_continuum)
    {
    }

    inline int get_fragment_pid() { return marley_utils::PHOTON; }
    inline int get_fragment_two_s() { return 2; }
    inline TMarleyParity get_fragment_parity() { return TMarleyParity(-1); }
    inline int get_fragment_Z() { return 0; }
    inline int get_fragment_A() { return 0; }
    inline double get_fragment_mass() { return 0.; }

    virtual void get_post_decay_parameters(double& Ex, int& two_J,
      TMarleyParity& Pi) = 0;

    // Return the final excitation energy (at the center of the bin
    // if this is a continuum de-excitation channel)
    virtual double get_bin_center_Exf() = 0;
};

class TMarleyDiscreteFragmentDecayChannel : public TMarleyFragmentDecayChannel {
  public:

    inline TMarleyDiscreteFragmentDecayChannel(const TMarleyFragment& frag,
      TMarleyLevel& f_lev): TMarleyFragmentDecayChannel(frag, false),
      final_level(f_lev)
    {
    }

    inline TMarleyLevel& get_final_level() {
      return final_level;
    }

    inline void get_post_decay_parameters(double& Ex, int& two_J,
      TMarleyParity& Pi)
    {
      Ex = final_level.get_energy();
      two_J = final_level.get_two_J();
      Pi = final_level.get_parity();
    }

    // Return the final excitation energy
    inline double get_bin_center_Exf() {
      return final_level.get_energy();
    }

  private:
    TMarleyLevel& final_level; // Pointer to final level in a nuclear decay scheme,
                               // nullptr if in the unbound continuum
    //int l;      // Orbital angular momentum of the outgoing fragment
};

class TMarleyDiscreteGammaDecayChannel : public TMarleyGammaDecayChannel {
  public:

    inline TMarleyDiscreteGammaDecayChannel(TMarleyLevel& f_lev):
      TMarleyGammaDecayChannel(false), final_level(f_lev)
    {
    }

    inline TMarleyLevel& get_final_level() {
      return final_level;
    }

    inline void get_post_decay_parameters(double& Ex, int& two_J,
      TMarleyParity& Pi)
    {
      Ex = final_level.get_energy();
      two_J = final_level.get_two_J();
      Pi = final_level.get_parity();
    }

    // Return the final excitation energy
    inline double get_bin_center_Exf() {
      return final_level.get_energy();
    }

  private:
    TMarleyLevel& final_level; // Pointer to final level in a nuclear decay scheme,
                               // nullptr if in the unbound continuum
    //int l;      // Orbital angular momentum of the outgoing fragment
};


// Continuum bin with boundaries Emin and Emax
class TMarleyContinuumBin {
  public:
    inline TMarleyContinuumBin(double E_min, double E_max,
      TMarleyGenerator& gener, bool use_upper_edge = false):
      gen(gener)
    {
      Emin = E_min;
      Emax = E_max;
      include_upper_edge = use_upper_edge;
    }

  protected:
    double Emin; // Minimum energy for this continuum bin
    double Emax; // Maximum energy for this continuum bin

    // Flag that instructs this object whether or not to include Emax
    // in the distribution of possible energies to sample. This is intended
    // to be used so that the last energy bin in the continuum sets it to
    // true and all others use false.
    bool include_upper_edge;

    TMarleyGenerator& gen;
    //int l;      // Orbital angular momentum of the outgoing fragment
};

// Fragment emission decay channel to a continuum bin with boundaries Emin and Emax
class TMarleyContinuumFragmentDecayChannel : public TMarleyFragmentDecayChannel,
  public TMarleyContinuumBin
{
  public:
    inline TMarleyContinuumFragmentDecayChannel(const TMarleyFragment& frag,
      TMarleyGenerator& gener, const TMarleySphericalOpticalModel& optmod,
      double E_min, double E_max, double mconst, double mfgs, double migs,
      bool use_upper_edge = false): TMarleyFragmentDecayChannel(frag, true),
      TMarleyContinuumBin(E_min, E_max, gener, use_upper_edge), om(optmod)
    {
      Mconst = mconst;
      Mfgs = mfgs;
      Migs = migs;
    }

    void get_post_decay_parameters(double& Ex, int& two_J, TMarleyParity& Pi);

    // Return the final excitation energy at the center of this bin
    inline double get_bin_center_Exf() {
      return Emin + ((Emax - Emin) / 2.0);
    }

  private:
    // Parameters used to determine outgoing projectile kinetic energy
    double Mconst;
    double Mfgs;
    double Migs;

    // Optical model to use for sampling final spin-parity
    const TMarleySphericalOpticalModel& om;
};

class TMarleyContinuumGammaDecayChannel : public TMarleyGammaDecayChannel,
  public TMarleyContinuumBin
{
  public:
    inline TMarleyContinuumGammaDecayChannel(TMarleyGenerator& gener,
      int Zinitial, int Ainitial, double E_min, double E_max,
      bool use_upper_edge = false): TMarleyGammaDecayChannel(true),
      TMarleyContinuumBin(E_min, E_max, gener, use_upper_edge)
    {
      Zi = Zinitial;
      Ai = Ainitial;
    }

    void get_post_decay_parameters(double& Ex, int& twoJ, TMarleyParity& Pi);

    // Return the final excitation energy at the center of this bin
    inline double get_bin_center_Exf() {
      return Emin + ((Emax - Emin) / 2.0);
    }

  private:
    int Zi, Ai; // Initial state atomic and mass numbers
};