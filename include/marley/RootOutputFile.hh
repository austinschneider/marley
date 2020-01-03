#pragma once

// standard library includes
#include "marley/OutputFile.hh"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

namespace marley {

  class RootOutputFile : public OutputFile {
    public:

      RootOutputFile(const std::string& name, const std::string& format,
        const std::string& mode, bool force = false);

      virtual ~RootOutputFile() = default;

      inline int_fast64_t bytes_written() override {
        return file_->GetBytesWritten();
      }

      // This function is a no-op for the ROOT format (we will take
      // care of saving this information when we write the generator
      // state variables)
      inline virtual void write_flux_avg_tot_xsec(double /*avg_tot_xsec*/)
        override {}

    private:

      virtual void open() override;

      void load_marley_headers();

      // Write the internal state string of the random number generator to
      // disk, as well as the generator's JSON configuration. This will allow
      // MARLEY to resume event generation from where it left off with no
      // loss of consistency. This trick is based on
      // http://tinyurl.com/hb7rqsj
      void write_generator_state(const marley::JSON& json_config,
        const marley::Generator& gen, const long /*num_events*/) override;

      // Clean up once we're done with the ROOT file. Save some information
      // about the current generator configuration as we clean things up.
      virtual void close(const marley::JSON& json_config,
        const marley::Generator& gen, const long dummy) override;

      virtual bool resume(std::unique_ptr<marley::Generator>& gen,
        long& num_previous_events) override;

      virtual void write_event(const marley::Event* event) override;

      // TFile object that will be used to access the ROOT file
      std::unique_ptr<TFile> file_ = nullptr;

      // This is a bare pointer, but ROOT will associate it with file_, so we
      // don't want to delete it ourselves or let a smart pointer do it.
      TTree* tree_ = nullptr;
  };

}
