#ifndef VOPENCL_LLG_HEUN_HPP_
#define VOPENCL_LLG_HEUN_HPP_

#include "internal.hpp"
#include "opencl_include.hpp"
#include "typedefs.hpp"

namespace vopencl
{
   namespace internal
   {
      struct heun_parameter_t
      {
         vcl_real_t prefactor;
         vcl_real_t lambda_times_prefactor;
      };

      namespace llg
      {
         extern bool initialized;

         extern cl::Buffer x_spin_array;
         extern cl::Buffer y_spin_array;
         extern cl::Buffer z_spin_array;

         extern cl::Buffer dS_x_array;
         extern cl::Buffer dS_y_array;
         extern cl::Buffer dS_z_array;

         extern cl::Buffer heun_parameters_device;

         void init(void);
         void step(void);
      }
   }
}

#endif // VOPENCL_LLG_HEUN_HPP_