#ifndef BaseCycleContainer_H
#define BaseCycleContainer_H

#include "SFrameTools/include/Objects.h"
#include "SFrameTools/include/ReconstructionHypothesis.h"

/**
 *  @short container that contains all objects of the actual event
 * 
 *  @author Thomas Peiffer
 */


struct BaseCycleContainer{
  //input variables
  int run;
  int luminosityBlock;
  int event;

  float rho;

  bool isRealData;
  //bool HBHENoiseFilterResult;

  float beamspot_x0;
  float beamspot_y0;
  float beamspot_z0;
  
  std::vector< Electron >* electrons;  
  std::vector< Muon >* muons;
  std::vector< Tau >* taus;
  std::vector< Photon >* photons;
  std::vector< PrimaryVertex >* pvs;
  std::vector< Jet >* jets;
  std::vector< TopJet >* topjets;
  std::vector< GenTopJet >* topjetsgen;
  std::vector< TopJet >* prunedjets;
  std::vector< GenParticle >* genparticles;
  std::vector< PFParticle >* pfparticles;
  std::vector< Particle>* genjets;

  MET* met;
  
  GenInfo* genInfo;
 
  std::vector<std::string>* triggerNames;
  std::vector<bool>* triggerResults;

  //use this vector since triggerNames is only filled for first event of new run
  std::vector<std::string> triggerNames_actualrun;

  std::vector< ReconstructionHypothesis >* recoHyps;
  
};


#endif
