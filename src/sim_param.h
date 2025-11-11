#ifndef SIM_PARAM_H_
#define SIM_PARAM_H_

#include "biodynamo.h"


namespace bdm {

// My simulation parameters
struct SimParam : public ParamGroup {
  BDM_PARAM_GROUP_HEADER(SimParam, 1);

  // Biological Parameters
  //  Growth rate 
  real_t mu = 3.3e-4;
  //  EPS production fraction
  real_t f = 0.30;
  //  Threshold for division
  real_t th_g = 1.3;
  //  Threshold for EPS excretion 
  real_t d_eps = 1.2;
  //  Enable quorum-sensing for cells
  bool QS = true;
  //  Quorum-sensing search radius
  real_t r_sense = 3;
  //  Quorum-sensing threshold
  int th_sense = 52;


  // Physical Parameters
  //  Spring stiffness
  real_t kn = 1e+17;
  //  Coefficient of restitution
  real_t e = 1e-30;
  //  Agent-agent Hamaker constant
  real_t cell_H = 1.6e+9;
  //  Agent-surface Hamaker constant
  real_t surface_H = 2.3e+9;
  //  Minimum separation distance
  real_t hmin = 5e-3;

  // Fluid Parameters
  //  Bulk fluid velocity
  real_t Ub = 1.86e3;
  //  Fluid dynamic viscosity
  real_t visc = 1e+6;
  //  Fluid density
  real_t rhof = 1000;

  // Resolution Parameters
  //  Biological timestep
  real_t biology_time_step = 550;

};

}  // namespace bdm

#endif  // SIM_PARAM_H_

