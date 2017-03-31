//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Rory Pond and Richard F L Evans 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "config.hpp"

// config module headers
#include "internal.hpp"

namespace config{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{
      //------------------------------------------------------------------------
      // Shared variables inside config module
      //------------------------------------------------------------------------
      bool output_atoms_config=false;
      int output_atoms_config_rate=1000;

      int output_rate_counter_coords = 0;
      int total_output_atoms = 0;

      bool output_cells_config=false;
      int output_cells_config_rate=1000;

      double field_output_min_1=-10000.0;
      double field_output_max_1=-0.0;
      double field_output_min_2=0.0;
      double field_output_max_2=10000.0;

      double atoms_output_min[3]={0.0,0.0,0.0};
      double atoms_output_max[3]={1.0,1.0,1.0};

      data_format output_data_format = text;
      bool output_new = false;
      bool mpi_io = false;

      std::vector<int> local_output_atom_list(0);

   } // end of internal namespace

} // end of config namespace