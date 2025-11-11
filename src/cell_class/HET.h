#ifndef HET_H_
#define HET_H_

#include "biodynamo.h"
#include "core/agent/cell.h"
#include "cell_class/mycell.h"
#include <string>

namespace bdm {

// Bacteria class
class HET : public MyCell {
  BDM_AGENT_HEADER(HET, MyCell, 1);

 public:
  HET() {}
  explicit HET(const Real3& position) : Base(position) {}
  virtual ~HET() {}

  // Set cell volume.
  void SetVolume_cell(real_t Volume_cell) {
          Volume_cell_ = Volume_cell;
          UpdateVolume_total();
  }

  // Change cell volume
  void ChangeVolume_cell(real_t speed) {
    Volume_cell_ += speed;
    UpdateVolume_total();
  }

  // Get cell volume
  real_t GetVolume_cell() const { return Volume_cell_; }


  // Set EPS shell volume
  void SetVolume_eps(real_t Volume_eps) {
          Volume_eps_ = Volume_eps;
          UpdateVolume_total();
  }

  // Change EPS shell volume
  void ChangeVolume_eps(real_t speed) {
    Volume_eps_ += speed;
    UpdateVolume_total();
  }

  // Get EPS shell volume
  real_t GetVolume_eps() const { return Volume_eps_; }


  // Preserve base volume pointer for the cell class in BioDynaMo for
  // visualisation and other core functions.
  void UpdateVolume_total() {
          SetVolume(Volume_cell_ + Volume_eps_);
  }

 private:
  // A cell has two volumes which make up its total volume: 
  // Cell volume which is the biological cell itself
  real_t Volume_cell_ = 0.0;
  // EPS shell volume which is the shell of EPS which surrounds the 
  // cell.
  real_t Volume_eps_ = 0.0;
};

}; // namespace bdm

#endif
