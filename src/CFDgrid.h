#ifndef CFDGRID_H_
#define CFDGRID_H_

#include "biodynamo.h"
#include <cmath>
#include <string>

namespace bdm {

// Fluid flow field class
class CFDgrid {
 public:
  CFDgrid(real_t xmin, real_t xmax, size_t xres,
          real_t ymin, real_t ymax, size_t yres,
          real_t zmin, real_t zmax, size_t zres);
  ~CFDgrid();
  BDM_CLASS_DEF(CFDgrid, 1);

  // Retrieve field information
  static CFDgrid* GetInstance();

  // Update flow velocity values in field
  void UpdateCFD();

  // Update solid fraction values in field
  void UpdateAlpha();

  // Retrieve fluid velocity at given position
  Real3 GetFluidVelocity(Real3 position);

  // Retrieve solid fraction at given position
  real_t GetAlpha(Real3 position);

 private:
  // Surface instance
  static CFDgrid* grid_;
  // Flow velocity values
  ParallelResizeVector<Real3> cfd_ = {};
  // Solid fraction values
  ParallelResizeVector<real_t> alpha_ = {};
  // Upper bound for field grid (x, y, z)
  std::array<real_t, 3> low_bound_;
  // Lower bound for field grid (x, y, z)
  std::array<real_t, 3> high_bound_;
  // Total number of voxels in field grid
  size_t total_box_num_;
  // Resolution of field grid (x, y, z)
  std::array<size_t, 3> res_;
  // Voxel length in field grid (x, y, z)
  std::array<real_t, 3> box_length_;
};

}; // namespace bdm

#endif

