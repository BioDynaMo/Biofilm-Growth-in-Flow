// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef BIOFILM_FORMATION_H_
#define BIOFILM_FORMATION_H_

#include "biodynamo.h"

#include "cell_class/HET.h"
#include "cell_class/EPS.h"
#include "RelaxAgents.h"
#include "PhysicalSteadyState.h"
#include "CFDgrid.h"

#include "biology_functions/Growth.h"
#include "biology_functions/MakeEPS.h"
#include "biology_functions/sensing.h"

#include "core/environment/uniform_grid_environment.h"

namespace bdm {

inline int Simulate(int argc, const char** argv) {

  // Base simulation parameters needed to run
  auto set_param = [&](Param* param) {
          param->use_progress_bar = true;
          param->bound_space = Param::BoundSpaceMode::kOpen;
	  param->min_bound = 0;
	  param->max_bound = 1000;

          param->simulation_time_step = 1e-8;
  };
  Param::RegisterParamGroup(new SimParam());
  Simulation simulation(argc, argv, set_param);


  // Gather simulation pointers
  auto* rm = simulation.GetResourceManager();  // Resouce manager
  auto* scheduler = simulation.GetScheduler(); // Scheduler
  auto* param = simulation.GetParam();         // Base BioDynaMo parameters
  const auto* sparam = param->Get<SimParam>(); // My simulation parameters


  // Unschedule default operations that are not needed
  scheduler->UnscheduleOp(scheduler->GetOps("load balancing")[0]);
  scheduler->UnscheduleOp(scheduler->GetOps("bound space")[0]);
  scheduler->UnscheduleOp(scheduler->GetOps("continuum")[0]);
  // The RelaxAgents and PhysicalSteadyState operations effectively replace this
  scheduler->UnscheduleOp(scheduler->GetOps("mechanical forces")[0]);


  // Add RelaxAgents operation to schedule
  OperationRegistry::GetInstance()->AddOperationImpl(
      "RelaxAgents", OpComputeTarget::kCpu, new RelaxAgents());
  scheduler->ScheduleOp(NewOperation("RelaxAgents"), OpType::kSchedule);

  // Add PhysicalSteadyState operation to schedule
  OperationRegistry::GetInstance()->AddOperationImpl(
      "PhysicalSteadyState", OpComputeTarget::kCpu, new PhysicalSteadyState());
  scheduler->ScheduleOp(NewOperation("PhysicalSteadyState"), OpType::kSchedule);


  // Set uniform grid environment for neighbour searching
  auto* env = dynamic_cast<UniformGridEnvironment*>(simulation.GetEnvironment());
  env->SetBoxLength(3);
  env->SetDetermineSimSize(false);


  // Initialise fluid dynamics voxel grid
  auto* cfdgrid = new CFDgrid(0, 200, 100,
                              0, 50, 25,
                              0, 50, 25);
  // Setup fluid flow field
  cfdgrid->UpdateCFD();


  // Initialise starting bacteria and EPS particles
  // -------------------------------------------------------
  // This reads the file containing the information of all starting bacteria
  // and EPS particles
  std::string file = "initial_chunk.csv";
  int number_of_lines = 0;
  std::string line;
  std::ifstream myfile(file);
  while (std::getline(myfile, line))
      ++number_of_lines;

  std::ifstream inFile;
  inFile.open(file.c_str());
  real_t x, y, z, VC, VE, rho;
  int type;
  for (int i=0; i < number_of_lines; i++) {
    inFile >> type >> x >> y >> z >> VC >> VE >> rho;
    Real3 position = {x, y, z};
  // --------------------------------------------------------
    // Bacteria
    if (type == 1) {
        auto* cell = new HET(position);  // Initialise bacteria
	cell->SetVolume_cell(VC);        // Cell volume
	cell->SetVolume_eps(VE);         // EPS shell volume
	cell->SetDensity(1000);          // Cell density
	cell->SetVelocity({0,0,0});      // Start with no velocity
	cell->SetCellType(1);            // Label as bacteria for cell forces

	// EPS production behaviour. Added before growth and division behaviour
	// as both depend on the cell volume at the start of the timestep, which
	// is then modified in the growth and division behaviour.
	cell->AddBehavior(new MakeEPS(sparam->d_eps, sparam->f) );
	// Growth and division behaviour
        cell->AddBehavior(new Growth(sparam->th_g, sparam->f) );

	// Quorum-sensing behaviour
	bool QS = sparam->QS;
        if (QS) {
            cell->AddBehavior(new Sensing(sparam->r_sense, sparam->th_sense) );
        }

	rm->AddAgent(cell); // Add cell to simulation
    }
    // EPS
    else if (type == 0) {
        auto* cell = new EPS(position); // Initialise EPS particle
        cell->SetVolume(VE);            // EPS volume
        cell->SetDensity(1000);         // EPS density
        cell->SetVelocity({0,0,0});     // Start with no velocity
        cell->SetCellType(0);           // Label as EPS for cell forces
	rm->AddAgent(cell);             // Add EPS to simulation
    }
  }
  inFile.close();


  // Print list of operations for scheduler
  scheduler->PrintInfo(std::cout);

  // Run simulation
  simulation.GetScheduler()->Simulate(200);


  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

}  // namespace bdm

#endif  // PHYSICAL_STEADY_STATE_H_
