#ifndef PHYSICAL_STEADY_STATE_H_
#define PHYSICAL_STEADY_STATE_H_

#include "core/operation/operation.h"
#include "core/resource_manager.h"
#include "cell_class/mycell.h"
#include "cell_class/HET.h"
#include "cell_class/EPS.h"
#include "MyMath.h"
#include "Boundary.h"

#include <chrono>

using namespace bdm;

inline void OutputAgentInfo(std::string file) {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    std::ofstream outdata;
    outdata.open(file);
    rm->ForEachAgent([&](Agent* agent) {
            auto* cell = dynamic_cast<MyCell*>(agent);
            int type = cell->GetCellType();
            Real3 pos = cell->GetPosition();
            Real3 vel = cell->GetVelocity();
            real_t d = cell->GetDiameter();
            outdata << pos[0] << " " << pos[1] << " " << pos[2] << " "
                    << vel[0] << " " << vel[1] << " " << vel[2] << " "
                    << d << " " << type << std::endl;
    });
    outdata.close();
}

struct PhysicalSteadyState : public StandaloneOperationImpl {
  BDM_OP_HEADER(PhysicalSteadyState);

  void operator()() override {
  
    // Gather simulation pointers
    auto* sim = Simulation::GetActive();                   // Current simulation
    auto* rm = sim->GetResourceManager();                  // Resource Manager
    const auto* sparam = sim->GetParam()->Get<SimParam>(); // My simulation parameters

    // Retrieve simulation timestep number
    int t = sim->GetScheduler()->GetSimulatedSteps();

    // Retrieve fluid flow grid
    auto* cfdgrid = CFDgrid::GetInstance();


    //std::string command = "mkdir output/steady_state_phy/"+std::to_string(t);
    //system(command.c_str());


    std::ofstream outdata;

    // If quorum-sensing is enabled, record all agents which switched phenotype and
    // removed them from the simulation.
    bool QS = sparam->QS;
    if (QS) {
        outdata.open("output/removed/QS/agent_data_"+std::to_string(t)+".csv");
        rm->ForEachAgent([&](Agent* agent) {
            auto* cell = dynamic_cast<MyCell*>(agent);
            int type = cell->GetCellType();
            if (type == 2) {
                Real3 pos = cell->GetPosition();
                Real3 vel = cell->GetVelocity();
                real_t d = cell->GetDiameter();
                outdata << pos[0] << " " << pos[1] << " " << pos[2] << " "
                        << vel[0] << " " << vel[1] << " " << vel[2] << " "
                        << d << " " << type << std::endl;

                cell->RemoveFromSimulation();
            }
        });
        outdata.close();
    }


    int i = 0;

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

    // Function to apply net force on agent and check for being outside simulation
    // space.
    auto ApplyNetForce = L2F([&](Agent* agent) {
	// Retrieve agent information
        auto* cell = dynamic_cast<MyCell*>(agent);

	// Update velocity and move agent
	Real3 f = cell->GetFnet();
        cell->UpdateVelocity(f, 1e-8, 1e+7);

	// Label all agents outside the simulation space
	int outside = ExitSim(cell);
        cell->SetOutside(outside);

    });

    // Tracker for finding steady state
    bool steady_state = false;

    // Particle-laden flow loop
    while (steady_state == false) {

	// Reset execution context
        const auto& all_exec_ctxts = sim->GetAllExecCtxts();
        all_exec_ctxts[0]->TearDownIterationAll(all_exec_ctxts);
        sim->GetEnvironment()->ForcedUpdate();
        all_exec_ctxts[0]->SetupIterationAll(all_exec_ctxts);
        sim->GetEnvironment()->Update();
        
	//if (i%1000 == 0) {
	//    OutputAgentInfo("output/steady_state_phy/"+std::to_string(t)+"/agent_data_"+std::to_string(i)+".csv");
	//}

	// Update solid fraction field every 100 timesteps
	if (i%100 == 0) {
	    cfdgrid->UpdateAlpha();
	}

	// Calculate net force for all agents
        rm->ForEachAgentParallel(CalculateNetForce);

	// Checking if each agent is static cannot be done in parallel, so every 1000 timesteps
	// the apply net force is run in series instead of parallel.
	if (i%1000 == 0) {
	    // Tracker for if the system is relaxed
	    int relaxed = 1;
            rm->ForEachAgent([&](Agent* agent) {
		// Retrieve agent information
                auto* cell = dynamic_cast<MyCell*>(agent);

		// Reset static for all agents at the start of the particle-laden flow
	        if (i == 0) {
                    cell->ResetStatic();
                }

		// Update velocity and move agents
	        Real3 f = cell->GetFnet();
	        cell->UpdateVelocity(f, 1e-8, 1e+6);

		// Label all agents outside the simulation space
		int outside = ExitSim(cell);
                cell->SetOutside(outside);

		// If all agents are either static or outside the simulation space, the system
		// is considered relaxed
	        if (cell->GetStatic() == 0 and outside==0) {
                    relaxed = 0;
                }

            });

	    // If the system is relaxed, the steady state has been found
	    if (relaxed == 1) {
                steady_state = true;
            }
	} else {
	    rm->ForEachAgentParallel(ApplyNetForce);
	}

	// The simulation moves on if the steady state is not found after some time.
	if (i > 150000) {
	    steady_state = true;
	}
	

	i++;
    }
    std::cout << "Mechanical teady state found after " << i << std::endl;

    // Record all agents that are outside the simulation space (or still moving inside 
    // the simulation space if the max iterations is reached) and remove them from the 
    // simulation.
    outdata.open("output/removed/agent_data_"+std::to_string(t)+".csv");
    rm->ForEachAgent([&](Agent* agent) {
            auto* cell = dynamic_cast<MyCell*>(agent);
            int type = cell->GetCellType();
            Real3 pos = cell->GetPosition();
            Real3 vel = cell->GetVelocity();
            real_t d = cell->GetDiameter();
	    int outside = cell->GetOutside();
	    if (outside == 1 or cell->GetStatic() == 0) {
            outdata << pos[0] << " " << pos[1] << " " << pos[2] << " "
                    << vel[0] << " " << vel[1] << " " << vel[2] << " "
                    << d << " " << type << std::endl;
	    }
	    if (outside == 1 or cell->GetStatic() == 0) {
	        cell->RemoveFromSimulation();
	    }
    });
    outdata.close();

    // Reset execution context
    const auto& all_exec_ctxts = sim->GetAllExecCtxts();
    all_exec_ctxts[0]->TearDownIterationAll(all_exec_ctxts);
    sim->GetEnvironment()->ForcedUpdate();
    all_exec_ctxts[0]->SetupIterationAll(all_exec_ctxts);
    sim->GetEnvironment()->Update();

    // Save agent data currently still inside the simulation
    OutputAgentInfo("output/agents/agent_data_"+std::to_string(t)+".csv");
  }
};

#endif  // PHYSICAL_STEADY_STATE_H_

