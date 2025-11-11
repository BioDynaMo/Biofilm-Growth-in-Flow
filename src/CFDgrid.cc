
#include "biodynamo.h"
#include "CFDgrid.h"

namespace bdm {

// Initialise grid
CFDgrid::CFDgrid(real_t xmin, real_t xmax, size_t xres,
		 real_t ymin, real_t ymax, size_t yres,
		 real_t zmin, real_t zmax, size_t zres){
 
  // Set lower and upper bound of grid
  low_bound_ = {xmin, ymin, zmin};
  high_bound_ = {xmax, ymax, zmax};

  // Set total number of voxels
  total_box_num_ = xres * yres * zres;

  // Set resolution
  res_ = {xres, yres, zres};

  // Set voxel length
  box_length_ = {(xmax-xmin)/xres,
                 (ymax-ymin)/yres,
                 (zmax-zmin)/zres};

  // Resize vectors containing flow velcoity and solid fraction values to fit
  // all voxels.
  cfd_.resize(total_box_num_);
  alpha_.resize(total_box_num_);

  // Output grid centers to csv to be read by OpenFoam
  std::ofstream outdata;
  outdata.open("cfd_grid.csv");
  for (size_t z=1; z<zres*2; z+=2) {
      for (size_t y=1; y<yres*2; y+=2) {
	  for (size_t x=1; x<xres*2; x+=2) {
               Real3 pos = {xmin + (box_length_[0]/2.)*x,
                            ymin + (box_length_[1]/2.)*y,
                            zmin + (box_length_[2]/2.)*z};
               outdata << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	  }
      }
  }

  // Allows the grid to be retrieved during simulation
  grid_ = this;

}

// Class destructor
CFDgrid::~CFDgrid(){
}
CFDgrid* CFDgrid::grid_ = nullptr;

// Getter for grid
CFDgrid* CFDgrid::GetInstance(){
  return grid_;
}

// Update flow velocity in grid using OpenFoam
void CFDgrid::UpdateCFD() {
  // Run OpenFoam command which returns flow velocity values
  std::string command = "updateGrid > /dev/null";
  system(command.c_str());

  // Read output values produced by OpenFoam into grid
  std::string inFileName = "cfdOutput.csv";
  std::ifstream inFile;
  inFile.open(inFileName.c_str());
  real_t Ux, Uy, Uz = 0;

  for (size_t i=0; i<total_box_num_; i++) {
      inFile >> Ux >> Uy >> Uz;
      Real3 U = {Ux, Uy, Uz};
      // OpenFoam values are normalised between 1 and 0, so all values
      // are scaled up to 1 m/s.
      cfd_[i] = U*1e+6;
  }
}

// Update solid fraction in grid
void CFDgrid::UpdateAlpha() {

  // Gather simulation pointers
  auto* sim = Simulation::GetActive();  // Current simulation
  auto* rm = sim->GetResourceManager(); // Resource manager

  // Set all solid fraction values to zero
  for (size_t i=0; i<total_box_num_; i++) {
      alpha_[i] = 0;
  }

  // Run through each agent to add solid fractions to each voxel
  rm->ForEachAgent([&](Agent* agent) {
      // Retrieve agent information
      auto* cell = dynamic_cast<Cell*>(agent);

      // Calculate agent volume and retrieve agent position
      real_t d = cell->GetDiameter();
      real_t Vp = (3./4.)*Math::kPi*(d/2.)*(d/2.)*(d/2.);
      Real3 position = cell->GetPosition();

      // Calculate ID of voxel the agent center resides in
      std::array<size_t, 3> box_coord;
      for (int i=0; i<3; i++) {
          box_coord[i] = std::floor((position[i] - low_bound_[i]) / box_length_[i]);
      }
      size_t idx = box_coord[2] * res_[0] * res_[1] +
               box_coord[1] * res_[0] + box_coord[0];
      
      // Calculate fraction of voxel taken up by agent and add to voxel solid fraction
      alpha_[idx] += Vp/(box_length_[0]*box_length_[1]*box_length_[2]);

      // Prevent solid fraction from exceeding 0.95 to prevent drag force approaching
      // infinity.
      if (alpha_[idx] > 0.95) {
          alpha_[idx] = 0.95;
      }
  });
}


// Retrieve fluid velocity with given position
Real3 CFDgrid::GetFluidVelocity(Real3 position) {
  std::array<size_t, 3> box_coord;
  for (int i=0; i<3; i++) {
      box_coord[i] = std::floor((position[i] - low_bound_[i]) / box_length_[i]);
  }
  size_t idx = box_coord[2] * res_[0] * res_[1] +
               box_coord[1] * res_[0] + box_coord[0];

  return cfd_[idx];
}

// Retrieve solid fraction with given position
real_t CFDgrid::GetAlpha(Real3 position) {
  std::array<size_t, 3> box_coord;
  for (int i=0; i<3; i++) {
      box_coord[i] = std::floor((position[i] - low_bound_[i]) / box_length_[i]);
  }
  size_t idx = box_coord[2] * res_[0] * res_[1] +
               box_coord[1] * res_[0] + box_coord[0];
  return alpha_[idx];
}

}; // namespace bdm
