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
#ifndef MAKEEPS_H_
#define MAKEEPS_H_

#include <cmath>

#include "biodynamo.h"
#include "cell_class/HET.h"
#include "cell_class/EPS.h"
#include "sim_param.h"

namespace bdm {

struct MakeEPS : public Behavior {
  BDM_BEHAVIOR_HEADER(MakeEPS, Behavior, 1);

  public:
   MakeEPS(real_t th=1, real_t f=0.2) : th_(th), f_(f) { AlwaysCopyToNew(); }

  virtual ~MakeEPS() {}

  // Initialise behaviour such that each cell can be given different excretion thresholds
  // and EPS production fractions.
  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    auto* other = event.existing_behavior;
    if (MakeEPS* gdbm = dynamic_cast<MakeEPS*>(other)) {
      th_ = gdbm->th_;
      f_ = gdbm->f_;
    } else {
      Log::Fatal("MakeEPS::EventConstructor",
                 "other was not of type MakeEPS");
    }
  }

  void Run(Agent* so) override {

      // Gather simulation pointers
      auto* sim = Simulation::GetActive();                   // Current simulation
      auto* ctxt = sim->GetExecutionContext();               // Execution context
      auto* random = sim->GetRandom();                       // Random number generator
      const auto* sparam = sim->GetParam()->Get<SimParam>(); // My simulation parameters

      // Retrieve cell information     
      auto* cell = dynamic_cast<HET*>(so);
      
      // Set biological time step and growth rate
      real_t biology_time_step = sparam->biology_time_step;
      real_t mu = sparam->mu;

      // Increase EPS shell volume based on simulation parameters and current cell volume
      cell->ChangeVolume_eps( f_ * mu * cell->GetVolume_cell() * biology_time_step );

      // Retrieve ratio between cell + EPS shell diameter and cell diameter
      real_t ratio = cell->GetDiameter() / ( 2.0 * cbrt( (3.0 * cell->GetVolume_cell()) / (4.0 * Math::kPi) ) );
      
      // If ratio is above EPS excretion threshold
      if (ratio >= th_ ) {

	  // Initialise EPS particle
          auto* eps = new EPS({0,0,0});

	  // Entire shell is excreted, leaving the cell with no shell
          eps->SetVolume( cell->GetVolume_eps() );
          cell->SetVolume_eps(real_t(0.0));
	      
	  // The EPS particles is 'pushed' in a random direction in the 3D polar coordinate
          // axis. The cell is 'pushed' in the opposite direction so that the center
          // of mass of the cell and EPS particle is in the same position as the original cell.
          real_t phi = Math::kPi * random->Uniform(-1, 1);
          real_t theta = Math::kPi * random->Uniform(0, 1);
          real_t r = ( cell->GetDiameter() + eps->GetDiameter() ) * 0.15;

          real_t x_coord = std::sin(theta) * std::cos(phi);
          real_t y_coord = std::sin(theta) * std::sin(phi);
          real_t z_coord = std::cos(theta);
          Real3 coords = {x_coord, y_coord, z_coord};

	  eps->SetPosition(cell->GetPosition() + (coords * r) );
          cell->SetPosition(cell->GetPosition() - (coords * r) );

	  // Set the density of the EPS particle to the same as the cell and set the cell type.
          eps->SetDensity(cell->GetDensity());
	  eps->SetCellType(0);
	  // Set EPS particle velocity to zero.
	  eps->SetVelocity({0,0,0});

	  // Add the EPS particle to the simulation.
          ctxt->AddAgent(eps);
      }
  }

  private:
   real_t th_;
   real_t f_;
   real_t Y_;
};


}; // namespace bdm

#endif  // MAKEEPS_H_
