#ifndef RELAX_AGENTS_H_
#define RELAX_AGENTS_H_

#include "core/operation/operation.h"
#include "core/resource_manager.h"
#include "cell_class/mycell.h"
#include "MyMath.h"
#include "Boundary.h"

using namespace bdm;

// Operation which relaxes the system in preparation for the particle-laden flow 
struct RelaxAgents : public StandaloneOperationImpl {
  BDM_OP_HEADER(RelaxAgents);

  void operator()() override {

    std::cout << "Biofilm grown" << std::endl;
  

    // Gather simulation pointers
    auto* sim = Simulation::GetActive();  // Current Simulation
    auto* rm = sim->GetResourceManager(); // Resource manager

    // Retrieve fluid flow grid
    auto* cfdgrid = CFDgrid::GetInstance();


    // Function to calculate net force on agent
    auto CalculateNetForce = L2F([&](Agent* agent) {
	// Retrieve agent information
        auto* cell = dynamic_cast<MyCell*>(agent);

	// Net force from neighbouring agents
        Real3 F_cells = cell->CalculateCellForce();

	// Net force from surface
        Real3 F_surface = cell->CalculateSurfaceForce();

	// Drag force
        Real3 F_drag = cell->CalculateDragForce();

	// Net force from simulation boundaries (excluding surface)
	Real3 F_boundary = BoundaryForce(cell);

	// Calculate net force
        Real3 F_net = F_cells + F_surface + F_drag + F_boundary;

	// Set net force for agent
        cell->SetFnet(F_net);
    });

    // Set physical timestep size and max velocity for agents
    real_t timestep = 1e-7;
    real_t max_speed = 1e+3;

    // Counter for number of iterations to reach relaxation
    int i = 0;
    // Tracker for finding steady state
    bool steady_state = false;

    // Relaxation loop
    while (steady_state == false) {

        // Reset execution context
        const auto& all_exec_ctxts = sim->GetAllExecCtxts();
        all_exec_ctxts[0]->TearDownIterationAll(all_exec_ctxts);
        sim->GetEnvironment()->ForcedUpdate();
        all_exec_ctxts[0]->SetupIterationAll(all_exec_ctxts);
        sim->GetEnvironment()->Update();

	// Track for maximum net force any agent is experiencing
	real_t max_force = 0;

	// Update solid fraction field every 100 iterations
	if (i%100 == 0) {
	    cfdgrid->UpdateAlpha();
	}

	// Calculate net force for each agent
        rm->ForEachAgentParallel(CalculateNetForce);

	// Apply net force for each agent and calculate net force
        rm->ForEachAgent([&](Agent* agent) {
	    // Retrieve agent information
            auto* cell = dynamic_cast<MyCell*>(agent);

	    // Reset static for all agent if its the start of the simulation
	    // loop
	    if (i == 0) {
	        cell->ResetStatic();
	    }

	    // Apply net force and move agents
	    Real3 f = cell->GetFnet();
	    cell->UpdateVelocity(f, timestep, max_speed);

	    // Label all agents outside the simulation space
            int outside = ExitSim(cell);
	    cell->SetOutside(outside);

	    // Calculate the maximum net force any agent is experiencing
	    if ( mod(f) > max_force) {
	        max_force = mod(f);
	    }
        });
        
	// If the maximum force is less than 0.5 nN, the system is considered relaxed
	if (max_force < 5e+14) {
	    steady_state = true;
	}

	// Simulation moves on if relaxation is not found after a long enough time
	if (i > 100000) {
            steady_state = true;
        }

	i++;
    }
    std::cout << "Agents relaxed after " << i << std::endl;
  }
};

#endif  // MOVE_H_

