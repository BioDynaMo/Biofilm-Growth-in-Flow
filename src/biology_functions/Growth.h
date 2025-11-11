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
#ifndef GROWTH_H_
#define GROWTH_H_

#include "biodynamo.h"
#include "sim_param.h"
#include "cell_class/HET.h"
#include "biology_functions/MakeEPS.h"
#include "biology_functions/sensing.h"

namespace bdm {

struct Growth : public Behavior {
  BDM_BEHAVIOR_HEADER(Growth, Behavior, 1);

  public:
   Growth(real_t th=1.3, real_t f=0.5 ) : th_(th), f_(f) { AlwaysCopyToNew(); }

  virtual ~Growth() {}

  // Initialise behaviour such that each cell can be given different division thresholds
  // and EPS production fractions.
  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    auto* other = event.existing_behavior;
    if (Growth* gdbm = dynamic_cast<Growth*>(other)) {
      th_   = gdbm->th_;
      f_    = gdbm->f_;
    } else {
      Log::Fatal("Growth::EventConstructor",
                 "other was not of type Growth");
    }
  }
  

  void Run(Agent* so) override {
      
      // Gather simulation pointers
      auto* sim = Simulation::GetActive();                   // Current simulation
      auto* ctxt = sim->GetExecutionContext();               // Execution context
      auto* random = sim->GetRandom();                       // Random number generator
      const auto* sparam = sim->GetParam()->Get<SimParam>(); // My simulation parameters
    
      // Retrieve cell information
      auto* cell = static_cast<HET*>(so);

      // Set biological time step and growth rate
      real_t biology_time_step = sparam->biology_time_step;
      real_t mu = sparam->mu;

      // Increase cell volume based on simulation parameters and current cell volume
      // (note: this does not include EPS shell)
      cell->ChangeVolume_cell( (real_t(1.0) - f_) * mu * cell->GetVolume_cell() * biology_time_step );

      // Get cell diameter (note: this does not include EPS shell)
      real_t d = 2.0 * cbrt( (3.0 * cell->GetVolume_cell() ) / ( 4.0 * Math::kPi ) );
      
      // If diameter is greater than division threshold
      if (d >= th_) {

	  // Initialise the daughter cell.
	  auto* daughter = new HET({0,0,0});

	  // The volume is not split 50/50 between the mother and daughter cell. Instead
	  // a random volume fraction of the mother's volume is given to the daughter.
	  real_t volume_ratio = random->Uniform(real_t(0.45), real_t(0.55));

	  // This fraction applies only applies to cell volume and not EPS shell.
	  // The mother keeps the entire EPS shell
	  daughter->SetVolume_cell( cell->GetVolume_cell() * volume_ratio );
	  cell->SetVolume_cell( cell->GetVolume_cell() * (real_t(1.0)-volume_ratio) );


	  // The daughter is 'pushed' in a random direction in the 3D polar coordinate
	  // axis. The mother is 'pushed' in the opposite direction so that the center 
	  // of mass of the two cells is in the same position as the original mother
	  // cell.
          real_t phi = Math::kPi * random->Uniform(-1, 1);
          real_t theta = Math::kPi * random->Uniform(0, 1);
          real_t r = ( cell->GetDiameter() + daughter->GetDiameter() ) * 0.15;

          real_t x_coord = std::sin(theta) * std::cos(phi);
          real_t y_coord = std::sin(theta) * std::sin(phi);
          real_t z_coord = std::cos(theta);
          Real3 coords = {x_coord, y_coord, z_coord};

          daughter->SetPosition(cell->GetPosition() + (coords * r) );
          cell->SetPosition(cell->GetPosition() - (coords * r) );


	  // The daughter cell is given the same density and cell type as mother cell.
          daughter->SetDensity(cell->GetDensity());
	  daughter->SetCellType(1);
	  // Daughter velocity is set to zero.
	  daughter->SetVelocity({0,0,0});

	  // As the threshold for EPS particle excretion is based on the ratio between
	  // cell and EPS shell volume, all cells which start with no EPS shell will 
	  // eventually secrete an EPS particle in the same timestep, regardless of the
	  // initial cell size. To prevent this, the EPS production fraction is slightly
	  // varied around the EPS production fraction for each cell.
	  real_t var = random->Uniform(real_t(-0.03), real_t(0.03));
	  // Daughter is given EPS production and growth and division behaviour.
	  // EPS production behaviour is added before growth and division behaviour
          // as both depend on the cell volume at the start of the timestep, which
          // is then modified in the growth and division behaviour.
	  daughter->AddBehavior(new MakeEPS(sparam->d_eps, sparam->f+var) );
          daughter->AddBehavior(new Growth(1.3, sparam->f+var) );

	  // If quorum-sensing is switched on, daughter cell is also given quorum-sensing
	  // behaviour.
	  bool QS = sparam->QS;
	  if (QS) {
              daughter->AddBehavior(new Sensing(sparam->r_sense, sparam->th_sense) );
	  }

	  // Add the new daughter cell to the simulation.
	  ctxt->AddAgent(daughter);
      }
  }

  private:
   real_t th_;
   real_t f_;
};


}; // namespace bdm

#endif  // GROWTH_H_
