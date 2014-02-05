#ifndef TriggerWeight_H
#define TriggerWeight_H

#include "SFrameTools/include/Objects.h"
#include "SFrameTools/include/BaseCycleContainer.h"

#include <algorithm>

#include "TH1.h"
#include "TFile.h"

class TriggerWeight{

 public:
  TriggerWeight(){};

  TriggerWeight(TString filename_mc, TString mode);

  ~TriggerWeight(){};

  double produceWeight(GenInfo* genInfo);
  
 private:

  TH1F* h_mc;
  TString m_mode;

};


#endif
