// Dear emacs, this is -*- c++ -*-
#ifndef HypothesisStatistics_H
#define HypothesisStatistics_H

#include <TString.h>
#include "SFrameTools/include/BaseCycleContainer.h"
#include "SFrameTools/include/ReconstructionHypothesis.h"
#include "core/include/SLogger.h"



/**
 *  @short calculate statistics for the comparison of two different reconstruction hypotheses
 *
 *  
 * 
 *  @author Thomas Peiffer
 */

class HypothesisStatistics{
 public:
  ///Default constructor
  HypothesisStatistics(std::string name);
  ///Default destructor
  ~HypothesisStatistics(){};

  ///set all hypothesis counters to 0
  void Reset();

  /// add two chosen hypotheses of the actual event; to be called in ExecuteEvent of a cycle
  void FillHyps(ReconstructionHypothesis *hyp1, ReconstructionHypothesis *hyp2);

  /// print the statistics for the comparison of the two hypothesis, to be called in EndInputData of a cycle
  void PrintStatistics();
  
 private:
  unsigned int m_ntotal;
  unsigned int m_identical;
  unsigned int m_toplep_identical;
  unsigned int m_tophad_identical;

  mutable SLogger m_logger;

};


#endif
