//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: sw766@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <fstream>
#include <cmath>

// Vampire headers
#include "anisotropy.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //---------------------------------------------------------------------------------
      // Function to add fourth order uniaxial anisotropy along vector e
      //
      //  Higher order anisotropies generally need to be described using spherical harmonics. The usual form (a
      //  series in S leads to cross pollution of terms, giving strange temperature dependencies.
      //
      //  The harmonics are described with Legendre polynomials with even order, which for 2nd, 4th and 6th are:
      //  ( http://en.wikipedia.org/wiki/Legendre_polynomials )
      //
      //  k2(sz) = - (1/2)  * (3sz^2 - 1)
      //  k4(sz) = - (1/8)  * (35sz^4 - 30sz^2 + 3)
      //  k6(sz) = - (1/16) * (231sz^6 - 315*sz^4 + 105sz^2 - 5)
      //
      //  The harmonics feature an arbritrary 2/3 factor compared with the usual form, and so in VAMPIRE these are
      //  renormalised to maintain consistency for the 2nd order terms.
      //
      //  The field induced by the harmonics is given by the first derivative w.r.t. sz. This can be projected onto
      //  any arbritrary direction ex,ey,ez allowing higher order anisotropy terms along any direction. This
      //  direction is shared with the other uniaxial anisotropy coefficients since they should not be used
      //  simultaneously.
      //
      //--------------------------------------------------------------------------------------------------------------
      void rot_fourth_order_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index){

         //// if not enabled then do nothing
         ////if(!internal::enable_uniaxial_fourth_order) return;

         //// constant factors
         //const double oneo8 = 1.0/8.0;

         //// rescaling prefactor
         //const double scale = oneo8*2.0/3.0; // Factor to rescale anisotropies to usual scale

         //// Loop over all atoms between start and end index
         //for(int atom = start_index; atom < end_index; atom++){

         //   // get atom material
         //   const int mat = atom_material_array[atom];

         //   const double sx = spin_array_x[atom]; // store spin direction in temporary variables
         //   const double sy = spin_array_y[atom];
         //   const double sz = spin_array_z[atom];

         //   const double ex = internal::ku_vector[mat].x; // phi is relative to this
         //   const double ey = internal::ku_vector[mat].y;
         //   const double ez = internal::ku_vector[mat].z;


         //   // get reduced anisotropy constant ku/mu_s
         //   const double ku4 = internal::ku4[mat];
         //   // hard coding the constant to test (TEMP)
         //   const double k4prime = 2.1E-21;//0.1*ku4;

         //   // cheat and assume that always aligned correctly (TEMP)
         //   //const double phi = atan2(sy,sx);
         //   const double xx = 1;
         //   const double xy = 0;
         //   const double xz = 0;
         //   const double cos_phi = sx*xx + sy*xy + sz*xz;
         //   const double cos_4phi = 1 - 8*cos_phi*cos_phi + 8*cos_phi*cos_phi*cos_phi*cos_phi;
         //   const double mod_const = k4prime * cos_4phi;


         //   const double sdote  = (sx*ex + sy*ey + sz*ez);
         //   const double sdote3 = sdote*sdote*sdote;


         //   // calculate field (double negative from scale factor and negative derivative)
         //   const double k4 = scale*mod_const*(140.0*sdote3 - 60.0*sdote);

         //   //field_array_x[atom] += ex*k4;
         //   //field_array_y[atom] += ey*k4;
         //   //field_array_z[atom] += ez*k4;


         //}

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add fourth order uniaxial anisotropy
      // E = 2/3 * - (1/8)  * (35sz^4 - 30sz^2 + 3)
      //---------------------------------------------------------------------------------
      double rot_fourth_order_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz){

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double ku4 = internal::ku4[mat];

         const double ex = internal::ku_vector[mat].x;
         const double ey = internal::ku_vector[mat].y;
         const double ez = internal::ku_vector[mat].z;

         // hard coding the constant to test (TEMP)
         //const double k4prime = 0.2*ku4;
         const double k4prime = internal::k4prime[mat];

         // cheat and assume that always aligned correctly (TEMP)
         const double xx = 1.0;
         const double xy = 0.0;
         const double xz = 0.0;
         const double cos_phi = sx*xx + sy*xy + sz*xz;
         //std::cout<<"COS_PHI:"<<cos_phi<<std::endl;
         const double cos_4phi = 1.0 - 8.0*cos_phi*cos_phi + 8.0*cos_phi*cos_phi*cos_phi*cos_phi;
         const double mod_const = k4prime * cos_4phi;


         const double sdote  = (sx*ex + sy*ey + sz*ez);//cos theta
         const double sdote2 = sdote*sdote;

         const double scrosse_x =  (sy*ez - sz*ey);
         const double scrosse_y = -(sx*ez - sz*ex);
         const double scrosse_z =  (sx*ey - sy*ex);
         const double sin2_theta = scrosse_x * scrosse_x
                                + scrosse_y * scrosse_y
                                + scrosse_z * scrosse_z;
         const double sin4_theta = sin2_theta * sin2_theta;
         // factor = 2/3 * -1/8 = -1/12 = -0.08333333333
         //return -0.08333333333*mod_const*(35.0*cos_phi*cos_phi*cos_phi*cos_phi - 30.0*cos_phi*cos_phi); //(35.0*sdote2*sdote2 - 30.0*sdote2);
         // cos_4phi from e^(4i phi) + e^(-4i phi)
         return -0.08333333333*4.1833*2*mod_const*sin4_theta*cos_4phi;

      }

   } // end of internal namespace

} // end of anisotropy namespace
