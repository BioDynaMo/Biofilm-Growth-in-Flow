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
#ifndef SENSING_H_
#define SENSING_H_

#include "biodynamo.h"
#include "cell_class/mycell.h"
#include "cell_class/HET.h"

namespace bdm {

struct Sensing : public Behavior {
  BDM_BEHAVIOR_HEADER(Sensing, Behavior, 1);

  public:
   Sensing(real_t r=5, int th=5 ) : r_(r), th_(th) { AlwaysCopyToNew(); }

  virtual ~Sensing() {}

  // Initialise behaviour such that each cell can be given different search radius
  // and cell count (although currently all cells are given the same parameters).
  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    auto* other = event.existing_behavior;
    if (Sensing* gdbm = dynamic_cast<Sensing*>(other)) {
      r_   = gdbm->r_;
      th_    = gdbm->th_;
    } else {
      Log::Fatal("Sensing::EventConstructor",
                 "other was not of type Sensing");
    }
  }
  

  void Run(Agent* so) override {

      // Gather simulation pointers
      auto* sim = Simulation::GetActive();     // Current simulation
      auto* ctxt = sim->GetExecutionContext(); // Execution context
    
      // Retrieve cell information
      auto* cell = static_cast<HET*>(so);

      // Retrieve cell position, QS search radius and QS threshold
      Real3 pos = cell->GetPosition();
      real_t search_radius = r_;
      int threshold = th_;

      // Pointer for number of cells found with search radius
      int tot = 0;

      // Nearest neighbour search function
      auto get_neighbor_number = L2F([&](Agent* other, real_t squared_distance) {
	    // Get neighbour information
            auto* neighbor = dynamic_cast<MyCell*>(other);

	    // Only cells are included in the count
	    int type = neighbor->GetCellType();
	    if (type != 0) {
	        tot += 1;
	    }
      });

      // Run nearest neighbour search function
      real_t search_squared = search_radius*search_radius;
      ctxt->ForEachNeighbor(get_neighbor_number, pos, search_squared);


      // If cell count is above threshold, change cell type to represent
      // phenotype switch.
      if (tot > threshold) {
          cell->SetCellType(2);
      }
      else {
          cell->SetCellType(1);
      }

  }

  private:
   real_t r_;
   int th_;
};


}; // namespace bdm

#endif  // GROWTH_H_
