//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "environment.hpp"
#include "random.hpp"
#include "../cells/internal.hpp"
#include "../micromagnetic/internal.hpp"
#include "micromagnetic.hpp"

// environment module headers
#include "internal.hpp"
#include "sim.hpp"
#include "cells.hpp"
#include "vmpi.hpp"

#include <iostream>
#include <math.h>



using namespace std;

namespace env = environment::internal;

void neighbours(int env, int mm, double overlap);

namespace environment{

   //----------------------------------------------------------------------------
   // Function to initialize environment module
   //----------------------------------------------------------------------------
   void initialize(   double system_dimensions_x, double system_dimensions_y, double system_dimensions_z){

         std::cout << "Initialising environment module" << std::endl;

         env::env_field_uv[0] = env::env_field_uv[0]*env::env_field;
         env::env_field_uv[1] = env::env_field_uv[1]*env::env_field;
         env::env_field_uv[2] = env::env_field_uv[2]*env::env_field;

         double size_x = system_dimensions_x + env::shift[0];
         double size_y = system_dimensions_y + env::shift[1];
         double size_z = system_dimensions_z + env::shift[2];

         double dx,dy,dz;

         if (env::dim[0] + 0.01 < size_x) dx = size_x + 0.01;
         else dx = env::dim[0] + 0.01;
         if (env::dim[1] + 0.01 < size_y) dy = size_y + 0.01;
         else dy = env::dim[1] + 0.01;
         if (env::dim[2] + 0.01 < size_z) dz = size_z + 0.01;
         else dz = env::dim[2] + 0.01;
         //convert to the total system size

         // determine number of stacks in x and y (global)
         env::num_cells_x =  dx/env::cell_size[0];
         env::num_cells_y =  dy/env::cell_size[1];
         env::num_cells_z =  dz/env::cell_size[2];

         //total number of cells and cell volume
         env::num_cells = env::num_cells_x*env::num_cells_y*env::num_cells_z;
         env::cell_volume = env::cell_size[0]*env::cell_size[1]*env::cell_size[2];

         std::cout << "Number of environment cells: " << env::num_cells << std::endl;

         //convert Ms from input to Ms = ms.V and Ku = ku.V

         env::Ms = env::Ms*env::cell_volume;
         env::ku = -env::ku*env::cell_volume;

         //resize arrays
         env::x_mag_array.resize(env::num_cells,0.0);
         env::y_mag_array.resize(env::num_cells,0.0);
         env::z_mag_array.resize(env::num_cells,0.0);

         env::cell_coords_array_x.resize(env::num_cells,0.0);
         env::cell_coords_array_y.resize(env::num_cells,0.0);
         env::cell_coords_array_z.resize(env::num_cells,0.0);

         env::dipole_field_x.resize(env::num_cells,0.0);
         env::dipole_field_y.resize(env::num_cells,0.0);
         env::dipole_field_z.resize(env::num_cells,0.0);

         env::neighbour_list_start_index.resize(env::num_cells,0.0);
         env::neighbour_list_end_index.resize(env::num_cells,0.0);

         env::list_env_cell_atomistic_cell.resize(cells::num_cells,0);
         env::env_cell_is_in_atomistic_region.resize(env::num_cells,0);

         environment_field_x.resize(cells::num_cells,0.0);
         environment_field_y.resize(cells::num_cells,0.0);
         environment_field_z.resize(cells::num_cells,0.0);

         //stores the max and min coordinates of all env cells
         std::vector <double > x_max(env::num_cells,0.0);
         std::vector <double > y_max(env::num_cells,0.0);
         std::vector <double > z_max(env::num_cells,0.0);
         std::vector <double > x_min(env::num_cells,0.0);
         std::vector <double > y_min(env::num_cells,0.0);
         std::vector <double > z_min(env::num_cells,0.0);


         //calculates the minimum and maximum positions of each cell for calcualtions of which cell the atomistic atoms are in later
         int cell = 0;
         for (int x = 0; x < env::num_cells_x; x++){
            double pos_x = env::cell_size[0]*x + env::cell_size[0]/2.0;
            for (int y = 0; y < env::num_cells_y; y++){
               double pos_y = env::cell_size[1]*y + env::cell_size[1]/2.0;
               for (int z = 0; z < env::num_cells_z; z++){
                  double pos_z = env::cell_size[2]*z + env::cell_size[2]/2.0;
                  env::cell_coords_array_x[cell] = pos_x;
                  env::cell_coords_array_y[cell] = pos_y;
                  env::cell_coords_array_z[cell] = pos_z;
                  x_max[cell] = pos_x + env::cell_size[0]/2.0;
                  x_min[cell] = pos_x - env::cell_size[0]/2.0;
                  y_max[cell] = pos_y + env::cell_size[1]/2.0;
                  y_min[cell] = pos_y - env::cell_size[1]/2.0;
                  z_max[cell] = pos_z + env::cell_size[2]/2.0;
                  z_min[cell] = pos_z - env::cell_size[2]/2.0;
            //      std::cout << x_min[cell] << '\t' << x_max[cell] << '\t' << y_min[cell] << "\t" << y_max[cell] << '\t' << z_min[cell] << '\t' << z_max[cell] << std::endl;
                  cell++;
               }
            }
         }

         std::cout << "Identifying environment cells which are atomistic" << std::endl;

         //loops over all atomistic cells to determine if the atomsitic simulation lies within an environment cell

         for(int lc=0; lc<cells::num_local_cells; lc++){
            int cell = cells::cell_id_array[lc];

            //calcaultes the x,y,z position for each atomsitic cell
            double x = cells::pos_and_mom_array[4*cell + 0] + env::shift[0];
            double y = cells::pos_and_mom_array[4*cell + 1] + env::shift[1];
            double z = cells::pos_and_mom_array[4*cell + 2] + env::shift[2];


            //loops over all environment cells to determine which cell this atomistic cell lies within
            for (int env_cell = 0; env_cell < env::num_cells; env_cell++){

               //if atom is within environment

               if (x < x_max[env_cell] && x > x_min[env_cell] && y < y_max[env_cell] && y > y_min[env_cell] && z < z_max[env_cell] && z > z_min[env_cell]){
                  //then the magnetisation of the enciroment cell is a sum of all atomistic atoms in that cell
                  //adds the cell to the list of atomistic cells
                  env::list_env_cell_atomistic_cell[cell] = env_cell;
                  std::cout << cell << '\t' << env_cell << std::endl;
                  //this cell is an atomistic cell.
                  env::env_cell_is_in_atomistic_region[env_cell] = 1;
                  env::x_mag_array[env_cell] += cells::mag_array_x[cell];
                  env::y_mag_array[env_cell] += cells::mag_array_y[cell];
                  env::z_mag_array[env_cell] += cells::mag_array_z[cell];
               }
            }
         }




         //calcualtes the neighbour lists for each cell.
         //if cells are neighbours add them to the neighbour list array_index
         //each cell has a start and end index index. the neighbours for each cell are within these index

         std::cout << "Calculating neighbour list for environment cells" << std::endl;
         // This is massively inefficient for large numbers of cells
         // better to store x,y,z cell associations and calculate neighbours directly - RE
         int array_index = 0;
         for (int celli = 0; celli < env::num_cells; celli ++){
            double xi = env::cell_coords_array_x[celli];
            double yi = env::cell_coords_array_y[celli];
            double zi = env::cell_coords_array_z[celli];
            env::neighbour_list_start_index[celli] = array_index;
            int b = 0;
            for (int cellj = 0; cellj < env::num_cells; cellj ++){
               int a = 0;
               double xj = env::cell_coords_array_x[cellj];
               double yj = env::cell_coords_array_y[cellj];
               double zj = env::cell_coords_array_z[cellj];
               double dx = sqrt((xi - xj)*(xi - xj));
               double dy = sqrt((yi - yj)*(yi - yj));
               double dz = sqrt((zi - zj)*(zi - zj));

               //calcualtes how many sides of the cube are touching neaarest neighbours have 1.
               if (dx == env::cell_size[0] && dy ==0 && dz == 0) a++;
               if (dy == env::cell_size[1] && dz ==0 && dx == 0) a++;
               if (dz == env::cell_size[2] && dy ==0 && dx == 0) a++;
               if (a == 1){
                  b++;
                  env::neighbour_list_array.push_back(cellj);                                        //if the interaction is non zero add the cell to the neighbourlist
                  env::neighbour_list_end_index[celli] = array_index;                                //the end index is updated for each cell so is given the value for the last cell.
                  array_index ++;
               }
            }
         }


         //calcualtes me
         double m_e = pow((env::Tc-sim::temperature)/(env::Tc),0.365);

         //initialise the direction of the cell magnetisation
         //if the initial spin configuration is set to random - give each cell a random magnetisation
         if (env::random_spins){

            for (int cell = 0; cell < env::num_cells; cell++){
               //if within the system
               if (env::cell_coords_array_x[cell] < env::dim[0] && env::cell_coords_array_y[cell] < env::dim[1] && env::cell_coords_array_z[cell] < env::dim[2] ) {
                  env::x_mag_array[cell] = m_e*mtrandom::gaussian()*env::Ms;
                  env::y_mag_array[cell] = m_e*mtrandom::gaussian()*env::Ms;
                  env::z_mag_array[cell] = m_e*mtrandom::gaussian()*env::Ms;
               }

               else{
                  env::x_mag_array[cell] = 0.0;
                  env::y_mag_array[cell] = 0.0;
                  env::z_mag_array[cell] = 0.0;
               }

            }
         }
         //if initial cell mangetisation is set as a direction set to that direction (normalised)
         else{
            const double normal = sqrt(env::initial_spin[0]*env::initial_spin[0] + env::initial_spin[1]*env::initial_spin[1]  + env::initial_spin[2]*env::initial_spin[2]);

            for (int cell = 0; cell < env::num_cells; cell++){
               //if within the system
               if (env::cell_coords_array_x[cell] < env::dim[0] && env::cell_coords_array_y[cell] < env::dim[1] && env::cell_coords_array_z[cell] < env::dim[2] ) {
                  env::x_mag_array[cell] = env::initial_spin[0]*env::Ms;
                  env::y_mag_array[cell] = env::initial_spin[1]*env::Ms;
                  env::z_mag_array[cell] = env::initial_spin[2]*env::Ms;

               }
               else{
                  env::x_mag_array[cell] = 0.0;
                  env::y_mag_array[cell] = 0.0;
                  env::z_mag_array[cell] = 0.0;
               }
            }
         }


         //adds all cells which are not within the atomistic section to the the none atomsitic cells list or the atomistic cells list
         for (int cell = 0; cell < env::num_cells; cell++){
            std::cout << cell << '\t' << env::x_mag_array[cell] << '\t' << env::y_mag_array[cell] << '\t' << env::z_mag_array[cell] <<std::endl;
            // calculate if cell is part of shield structure
            bool included = true;//internal::in_shield(env::cell_coords_array_x[cell], env::cell_coords_array_y[cell], env::cell_coords_array_z[cell]);
            if (env::x_mag_array[cell] + env::y_mag_array[cell] +env::z_mag_array[cell]  < 1e-27) included =false;
               // check environment cell is in shields and not atomistic

               if (env::env_cell_is_in_atomistic_region[cell] == 0 && included){
                  env::none_atomistic_cells.push_back(cell);
                  env::num_env_cells ++;
               }
               // if it is atomistic then add it to list of atomistic cells
               else if(env::env_cell_is_in_atomistic_region[cell] == 1){
                  env::atomistic_cells.push_back(cell);
               }
               // otherwise don't add it to anything
         }


         //initalise the demag fields
         int a = env::initialise_demag_fields();



         return;

      }

   } // end of environment namespace

void neighbours(int env, int mm, double overlap){
   environment::list_of_mm_cells_with_neighbours.push_back(mm);
   environment::list_of_env_cells_with_neighbours.push_back(env);
   environment::list_of_overlap_area.push_back(overlap);
   environment::num_interactions++;
   }
