#ifndef BOUNDARY_H_
#define BOUNDARY_H_


#include "biodynamo.h"
#include "cell_class/mycell.h"

namespace bdm {

// Calculate boundary force from simulation walls. Note that maximum x boundary is open
// and minimum z boundary is the surface, where the force is calculated in the custom
// agent class.
inline Real3 BoundaryForce(Agent* agent) {

  // Retrieve agent information
  auto* cell = dynamic_cast<MyCell*>(agent);
  
  // Simulation boundaries
  real_t lb = 0;
  real_t rb = 50;

  // Retrieve agent position and radius
  auto pos = cell->GetPosition();
  real_t radius = cell->GetDiameter() / 2.;

  // Initialise force vector
  Real3 fn = {0,0,0};

  // Overlap with minimum x boundary
  if ( pos[0] - radius < lb) {
      // Retrieve agent velocity and mass
      auto vel = cell->GetVelocity();
      real_t mass = cell->GetMass();

      // Unit distance between agent and boundary
      Real3 r_unit = {1, 0, 0};

      // Normal velocity between agent and boundary
      Real3 vn = DotProduct(vel, r_unit)*r_unit;

      // Overlap distance with boundary 
      real_t delta = radius - pos[0] + lb;

      // Retrieve spring stiffness constant and coefficient of restitution
      real_t kn = 1e+17;
      real_t e = 1e-30;

      // Calculate spring damping constant
      real_t gamman  = -2*std::log(e)/std::sqrt( ((Math::kPi*Math::kPi) + (std::log(e)*std::log(e))) * (mass/kn) );

      // Spring-dashpot force equation
      //      Contact force      damping force
      fn += (kn*delta*r_unit) - (gamman*mass*vn);
  }

  // Minimum y boundary
  if ( pos[1] - radius < lb) {
      // Retrieve agent velocity and mass
      auto vel = cell->GetVelocity();
      real_t mass = cell->GetMass();

      // Unit distance between agent and boundary
      Real3 r_unit = {0, 1, 0};

      // Normal velocity between agent and boundary
      Real3 vn = DotProduct(vel, r_unit)*r_unit;

      // Overlap distance with boundary
      real_t delta = radius - pos[1] + lb;

      // Retrieve spring stiffness constant and coefficient of restitution
      real_t kn = 1e+17;
      real_t e = 1e-30;

      // Calculate spring damping constant
      real_t gamman  = -2*std::log(e)/std::sqrt( ((Math::kPi*Math::kPi) + (std::log(e)*std::log(e))) * (mass/kn) );

      // Spring-dashpot force equation
      //      Contact force      damping force
      fn += (kn*delta*r_unit) - (gamman*mass*vn);
  }

  // Maximum y boundary
  if ( pos[1] + radius > rb) {
      // Retrieve agent velocity and mass
      auto vel = cell->GetVelocity();
      real_t mass = cell->GetMass();

      // Unit distance between agent and boundary
      Real3 r_unit = {0, -1, 0};

      // Normal velocity between agent and boundary
      Real3 vn = DotProduct(vel, r_unit)*r_unit;

      // Overlap distance with boundary
      real_t delta = radius + pos[1] - rb;

      // Retrieve spring stiffness constant and coefficient of restitution
      real_t kn = 1e+17;
      real_t e = 1e-30;

      // Calculate spring damping constant
      real_t gamman  = -2*std::log(e)/std::sqrt( ((Math::kPi*Math::kPi) + (std::log(e)*std::log(e))) * (mass/kn) );

      // Spring-dashpot force equation
      //      Contact force      damping force
      fn += (kn*delta*r_unit) - (gamman*mass*vn);
  }

  // Maximum z boundary
  if ( pos[2] + radius > rb) {
      // Retrieve agent velocity and mass
      auto vel = cell->GetVelocity();
      real_t mass = cell->GetMass();

      // Unit distance between agent and boundary
      Real3 r_unit = {0, 0, -1};

      // Normal velocity between agent and boundary
      Real3 vn = DotProduct(vel, r_unit)*r_unit;

      // Overlap distance with boundary
      real_t delta = radius + pos[2] - rb;

      // Retrieve spring stiffness constant and coefficient of restitution
      real_t kn = 1e+17;
      real_t e = 1e-30;

      // Calculate spring damping constant
      real_t gamman  = -2*std::log(e)/std::sqrt( ((Math::kPi*Math::kPi) + (std::log(e)*std::log(e))) * (mass/kn) );

      // Spring-dashpot force equation
      //      Contact force      damping force
      fn += (kn*delta*r_unit) - (gamman*mass*vn);
  }

  return fn;

}

// Determine if agent is outside simulation space
inline int ExitSim(Agent* agent) {

  // Retrieve agent information
  auto* cell = dynamic_cast<MyCell*>(agent);

  // Simulation dimensions
  Real3 lb = {0, 0, 0};
  Real3 rb = {100, 50, 50};

  // Retrieve agent position
  auto pos = cell->GetPosition();

  // Agent is outside if position is beyond the simulation dimensions
  if (pos[0] < lb[0] or pos[0] > rb[0]
     or pos[1] < lb[1] or pos[1] > rb[1]
     or pos[2] > rb[2]) {
      return 1;
  }
  else {
      return 0;
  }
}

}; // namespace bdm

#endif  // MOVE_H_

