#ifndef MYCELL_H_
#define MYCELL_H_

#include "biodynamo.h"
#include "core/agent/cell.h"
#include "sim_param.h"
#include "CFDgrid.h"
#include "MyMath.h"
#include <cmath>
#include <string>

namespace bdm {

// Base agent class used for the particle-laden flow. The bacteria and
// EPS cell classes are extensions of this class.
class MyCell : public Cell {
  BDM_AGENT_HEADER(MyCell, Cell, 1);

 public:
  MyCell() {}
  explicit MyCell(const Real3& position) : Base(position) {}
  virtual ~MyCell() {}

  // Setter and getter for agent type.
  void SetCellType(int cell_type) { cell_type_ = cell_type; }
  int GetCellType() const { return cell_type_; }

  // Setter and getter for agent velocity
  void SetVelocity(Real3 velocity) {
    velocity_ = velocity;
  }
  Real3 GetVelocity() { return velocity_; }

  // Update velocity and position based on net force, physical time step
  // and maximum velocity cap
  void UpdateVelocity(Real3 F_net, real_t dt, real_t vel_cap) {
    
    // Calculate change in velocity	  
    real_t mass = GetMass();
    Real3 deltaV = F_net/mass * dt;
    velocity_ += deltaV;

    // If velocity is above velocity cap, set velocity to velocity cap
    // while preserving direction
    if (mod(velocity_) > vel_cap) {
        velocity_ *= vel_cap/mod(velocity_);
    }

    // Update position based on new velocity
    UpdatePosition(velocity_ * dt);

    // Check if agent is considered static. Used for finding physical 
    // steady state.
    // Counter for number of physical time steps of particle-laden flow
    static_iter += 1;
    // Once number of physical time steps reach 1000.
    if (static_iter == 1000) {

	// Calculate distance between current position previous position
	// from 1000 time steps before.    
        Real3 pos = GetPosition();
        real_t distance = mod(pos - pos_old);

	// If distance is less than 0.1 micrometers, agent is considered static
        if (distance < 0.1) {
            static_ = 1;
        }
        else {
            static_ = 0;
        }

	// Set old position to current position and reset counter. 
        pos_old = pos;
        static_iter = 0;
    }
  }

  // Getter and resetter for static check
  int GetStatic() { return static_; }
  void ResetStatic() {
      static_ = 0;
      static_iter = 0;
      pos_old = {0,0,0};
  }

  // Setter and getter for check if agent is outside simulation space
  void SetOutside(int outside) {
      outside_ = outside;
  }
  int GetOutside() { return outside_; }

  // Setter and getter for agent velocity
  void SetFnet(Real3 Fnet) {
      fnet_ = Fnet;
  }
  Real3 GetFnet() { return fnet_; }


  // Drag force calculation
  Real3 CalculateDragForce() {
    
    // Gather simulation pointers
    auto* sim = Simulation::GetActive();                   // Current simulation
    const auto* sparam = sim->GetParam()->Get<SimParam>(); // My simulation parameters

    // Retrieve fluid dynamics grid
    auto* cfdgrid = CFDgrid::GetInstance();
    // Retrieve fluid velocity and solid volume fraction at agent position
    Real3 fluid_vel = cfdgrid->GetFluidVelocity(GetPosition());
    real_t alpha = cfdgrid->GetAlpha(GetPosition());

    // Retrieve fluid viscosity, density and bulk velocity
    real_t visc = sparam->visc;
    real_t rhof = sparam->rhof;
    real_t Ub = sparam->Ub;

    // Retrieve agent diameter
    real_t d = GetDiameter();

    // The fluid velocity field is normalised to be between 1 and 0 m/s to prevent agents from taking too long
    // to exit the simulation. The relative velocity is scaled back down to provide an accurate drag force
    // while maintaining a faster simulation. The trade-off is the agents having a much higher momentum than
    // they realisticallly would in the real environment.
    Real3 vrel = (fluid_vel - velocity_)*Ub/1e+6;

    // Calculate Reynolds number
    real_t Re = d*mod(vrel)/visc;

    // Calculate relative velocity correlation
    real_t A = pow(1-alpha, 4.14);
    real_t B = pow(1-alpha, 2.65);
    if (alpha >= 0.15) {
        B = 0.8*pow(1-alpha, 1.28);
    }
    real_t Vs = 0.5*(A - (0.06*Re) + std::sqrt((0.06*Re)*(0.06*Re) +
                                        (0.12*Re*((2*B) - A)) +
                                        (A*A)));

    // Calculate drag coefficient. The drag coefficient equation has been expanded and multiplied by the relative
    // velocity in order to prevent infinities when the relative velocity approaches zero.
    real_t Cd = (0.3939*mod(vrel)) + (6.048*std::sqrt(Vs*visc*mod(vrel))/std::sqrt(d)) + (23.04*Vs*visc/d);

    // Drag force equation. This only contains one relative velocity as the other is included in the drag 
    // coefficient equation.
    Real3 F_drag = (Math::kPi/8.) * (d*d*rhof/(Vs*Vs)) * Cd * vrel;

    return F_drag;
  }



  // Calculate net force from neighbouring agents
  Real3 CalculateCellForce() {

    // Gather simulation pointers
    auto* sim = Simulation::GetActive();                   // Current simulation
    auto* ctxt = sim->GetExecutionContext();               // Execution context
    const auto* sparam = sim->GetParam()->Get<SimParam>(); // My simulation parameters

    // Trackers for total contact force and attraction force
    Real3 cell_contact = {0,0,0};
    Real3 cell_attraction = {0,0,0};

    // Retrieve agent information
    Real3 p1 = this->GetPosition();        // Position
    Real3 v1 = this->GetVelocity();        // Velocity
    real_t R1 = this->GetDiameter() / 2.0; // Radius
    real_t m1 = this->GetMass();           // Mass
    real_t t1 = this->GetCellType();       // Type (cell or EPS particle)

    // Retrieve spring stiffness constant and coefficient of restitution
    real_t kn = sparam->kn; 
    real_t e = sparam->e;

    // Retrieve agent-agent Hamaker constant and minimum separation distance
    real_t H = sparam->cell_H;
    real_t hmin = sparam->hmin;

    // Nearest neighbour function for calculating force between agent and neighbour
    auto calculate_neighbor_forces = L2F([&](Agent* other, real_t squared_distance) {
		
	    // Retrieve neighbour information    
            auto* neighbor = dynamic_cast<MyCell*>(other);
            Real3 p2 = neighbor->GetPosition();        // Position
            Real3 v2 = neighbor->GetVelocity();        // Velocity
            real_t R2 = neighbor->GetDiameter() / 2.0; // Radius
            real_t m2 = neighbor->GetMass();           // Mass
	    real_t t2 = neighbor->GetCellType();       // Type (cell or EPS particle)

	    // Calculate relative velocity between agent and neighbour
            Real3 vrel = v1 - v2;

	    // Calculate unit vector for distance between agent and neighbour
            Real3 r = p1 - p2;
            real_t rmag = mod(r);
            Real3 r_unit = r/rmag;

	    // Calculate overlap distance between agent and neighbour
            real_t delta = R1 + R2 - rmag;
	    // If agents overlap, calculate contact force
            if (delta > 0) {

		// Calculate normal relative velocity with respect to distance unit vector
        	Real3 vn = DotProduct(vrel, r_unit)*r_unit;

		// Calculate effective mass
        	real_t meff = (m1*m2)/(m1+m2);

		// Calculate spring damping constant
        	real_t gamman  = -2*std::log(e)/std::sqrt( ((Math::kPi*Math::kPi) + (std::log(e)*std::log(e))) * (meff/kn) );

		// Spring-daspot force equation
		//           spring force       Damping force
                Real3 fn = (kn*delta*r_unit) - (gamman*meff*vn);

		// Add contact force to total contact force
        	cell_contact += fn;
            }

	    // Calculate effective radius
            real_t Re = R1*R2/(R1+R2);

	    // If agent surfaces are closer than minimum separation distance, set overlap
	    // distance to minimum separation distance (this applies during agent contact).
	    if (delta > -hmin){
                delta = hmin;
            }
	    // Van der Waals attraction force between spheres
	    Real3 fa = - H * Re / (6 * delta*delta) * r_unit;


	    // Multiply the attraction force by 3 orders of magnitude if either the agent or
	    // the neighbour is an EPS particle.
	    // Agent and neighbour are EPS particles
	    if (t1 == 0 and t2 == 0) {
	        fa *= 3125;
	    }
	    // Agent is EPS particle
	    else if (t1 == 1 and t2 == 0) {
	        fa *= 1000;
	    }
	    // Neighbour is EPS particle
	    else if (t1 == 0 and t2 == 1) {
	        fa *= 1000;
	    }

	    // Add attraction force to total attraction force
            cell_attraction += fa;
    });

    // Run nearest neighbour funtion. This will update the total contact force and total
    // attraction force from all neighbouring particles.
    ctxt->ForEachNeighbor(calculate_neighbor_forces, *this, 2.56);

    // Calculate net force from all agents
    Real3 cell_force = cell_contact + cell_attraction;

    return cell_force;
  }


  // Calculate net force from the surface
  Real3 CalculateSurfaceForce() {

     // Gather simulation pointers
     auto* sim = Simulation::GetActive();                   // Current simulation
     const auto* sparam = sim->GetParam()->Get<SimParam>(); // My simulation parameters

     // Retrieve physical timestep size
     real_t dt = 1e-8;

     // Initialise contact force and friction force
     Real3 surface_contact = {0,0,0};
     Real3 surface_friction = {0,0,0};

     // Retrieve agent information
     Real3 pos = GetPosition(); // Position
     real_t d = GetDiameter();  // Diameter
     Real3 vel = GetVelocity(); // Velocity
     real_t mass = GetMass();   // Mass
     real_t t = GetCellType();  // Cell type (cell or EPS particle)

     // Retrieve spring stiffness constant and coefficient of restitution
     real_t kn = sparam->kn;
     real_t e = sparam->e;

     // Retrieve agent-agent Hamaker constant and minimum separation distance
     real_t H = sparam->surface_H;
     real_t hmin = sparam->hmin;

     // Calculate normal spring damping constant
     real_t gamman  = -2*std::log(e)/std::sqrt( ((Math::kPi*Math::kPi) + (std::log(e)*std::log(e))) * (mass/kn) );
     
     // Calculate tangential spring stiffness and damping constant
     real_t kt = 0.1*kn;
     real_t gammat = 0.5*gamman;


     // Tracker for if surface makes contact
     int contact = 0;

     // Unit distance between agent and surface (minimum z boundary)
     Real3 r_unit = {0, 0, 1};

     // Calculate overlap distance between agent and surface
     real_t delta = (d/2.) - pos[2];
     // If agent contacts surface, calculate contact and friction force
     if (delta > 0) {

	 // Calculate normal velocity of agent with respect to unit distance
         Real3 vn = DotProduct(vel, r_unit)*r_unit;
	 // Calculate tangential velocity of agent with respect to unit distance
         Real3 vt = vel - vn;

	 // spring-dashpot normal force
	 //                   contact force      damping force
         surface_contact += (kn*delta*r_unit) - (gamman*mass*vn);

	 contact = 1;

	 // Calculate tangential displacement between agent and surface
	 // If agent is already in contact with surface, update displacement
	 if (surface_contact_ == 1) {
	     surface_distance_ += vt*dt;
	 }
	 // If agent is making contact in this timestep, set displacement to zero
	 if (surface_contact_ == 0) {
	     surface_distance_ = {0,0,0};
	     surface_contact_ = 1;
	 }

	 // Zero vector
	 Real3 Zero = {0,0,0};

	 // Spring-dashpot tangential force
	 //                           Friction force        Damping force
         surface_friction += Zero-(kt*surface_distance_) - (gammat*mass*vt);
     }

     // If agent is not in contact with surface, reset tangential displacement
     if (contact == 0) {
         surface_distance_ = {0,0,0};
	 surface_contact_ = 0;
     }

     // If distance between agent surface and surface is closer than minimum separation
     // distance, set overlap distance to minimum separation distance (this applies
     // during surface contact).
     if (delta > -hmin){
         delta = hmin;
     }

     // Van der Waals attraction force between sphere and infinite plane
     Real3 surface_attraction = - H * (d/2.) / (6 * delta*delta) * r_unit;

     // Multiply the attraction force by 3 orders of magnitude if the agent is an EPS
     // particle.
     if (t == 0) {
         surface_attraction *= 1000;
     }

     // Calculate net force from surface
     Real3 surface_force = surface_contact + surface_friction + surface_attraction;

     return surface_force;
  }

 private:
  // Agent type
  int cell_type_ = 0;
  // Agent velocity
  Real3 velocity_ = {0,0,0};
  // Net force on agent
  Real3 fnet_ = {0,0,0};
  // Tracker for whether agent is in contact with surface
  int surface_contact_ = 0;
  // Tangential displacement between agent and surface
  Real3 surface_distance_ = {0,0,0};
  // Tracker for whether agent is outside the simulation space
  int outside_ = 0;
  // Tracker for whether agent is considered static
  int static_ = false;
  // Counter for when next check for static will occur (resets every 1000 physical
  // timesteps)
  int static_iter = 0;
  // Previous position used in check for static. Updates every 1000 timesteps
  Real3 pos_old = {0,0,0};
};

}; // namespace bdm

#endif

