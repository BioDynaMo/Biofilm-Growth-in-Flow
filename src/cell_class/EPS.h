#ifndef EPS_H_
#define EPS_H_

#include "biodynamo.h"
#include "core/agent/cell.h"
#include "cell_class/mycell.h"
#include <string>

namespace bdm {

// EPS particle class. EPS particles do not require any extra functions 
// than the MyCell class. This class is only defined to help visualise 
// bacteria and EPS particles separately in paraview.
class EPS : public MyCell {                              
  BDM_AGENT_HEADER(EPS, MyCell, 1);

 public:
  EPS() {}
  explicit EPS(const Real3& position) : Base(position) {}
  virtual ~EPS() {}
 
};

}; // namespace bdm

#endif
