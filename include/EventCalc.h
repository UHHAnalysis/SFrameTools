// Dear emacs, this is -*- c++ -*-
#ifndef EventCalc_H
#define EventCalc_H

#include "core/include/SLogger.h"
#include "SFrameTools/include/Utils.h"
#include "SFrameTools/include/identifier.h"
#include "SFrameTools/include/TTbarGen.h"
#include "SFrameTools/include/LuminosityHandler.h"
#include "SFrameTools/include/BaseCycleContainer.h"

#include "TVector3.h"

/**
 *  @short class for the calculation of basic event variables
 *
 * 
 *
 *  @author Roman Kogler
 */

class EventCalc
{
 private:
  static EventCalc* m_instance;

 public:
  static EventCalc* Instance();

  void Reset();
  
  BaseCycleContainer* GetBaseCycleContainer();
  LuminosityHandler* GetLumiHandler();

  void SetLumiHandler(LuminosityHandler* lh);
  void SetBaseCycleContainer(BaseCycleContainer* bcc);

  int GetRunNum() {return (m_bcc ? m_bcc->run : -1);}
  int GetLumiBlock() {return (m_bcc ? m_bcc->luminosityBlock : -1);}
  int GetEventNum() {return (m_bcc ? m_bcc->event : -1);}
  bool IsRealData() {return (m_bcc ? m_bcc->isRealData : false);}
  //bool GetHBHENoiseFilterResult() {return (m_bcc ? m_bcc->HBHENoiseFilterResult : false);}

  float GetBeamSpotX0() {return (m_bcc ? m_bcc->beamspot_x0 : -999.);}
  float GetBeamSpotY0() {return (m_bcc ? m_bcc->beamspot_y0 : -999.);}
  float GetBeamSpotZ0() {return (m_bcc ? m_bcc->beamspot_z0 : -999.);}
  
  std::vector< Electron >* GetElectrons() {return (m_bcc ? m_bcc->electrons : NULL);}
  std::vector< Muon >* GetMuons() {return (m_bcc ? m_bcc->muons : NULL);}
  std::vector< Tau >* GetTaus() {return (m_bcc ? m_bcc->taus : NULL);}
  std::vector< Photon >* GetPhotons() {return (m_bcc ? m_bcc->photons : NULL);}
  std::vector< PrimaryVertex >* GetPrimaryVertices() {return (m_bcc ? m_bcc->pvs : NULL);}

  std::vector< Jet >* GetJets() { return (m_bcc ? m_bcc->jets : NULL);}
  std::vector< TopJet >* GetCAJets() {return (m_bcc ? m_bcc->topjets : NULL);}
  std::vector< TopJet >* GetPrunedCAJets() {return (m_bcc ? m_bcc->prunedjets : NULL);}
  std::vector< GenParticle >* GetGenParticles() {return (m_bcc ? m_bcc->genparticles : NULL);}
  std::vector< Particle >* GetGenJets() {return (m_bcc ? m_bcc->genjets : NULL);}
  /// returns all stored PF particles independently of their origin (to be used with care!)
  std::vector< PFParticle >* GetPFParticles() {return (m_bcc ? m_bcc->pfparticles : NULL);}

  /// PF particles clustered to all (top-)jets
  std::vector< PFParticle >* GetJetPFParticles();
  /// PF particles clustered to a specific jet
  std::vector< PFParticle > GetJetPFParticles( Jet* topjet);
  /// PF particles used to calculate lepton isolation
  std::vector< PFParticle >* GetIsoPFParticles();
  /// PF particles used to calculate pile-up contribution to lepton isolation
  std::vector< PFParticle >* GetPUIsoPFParticles();

  MET* GetMET() {return (m_bcc ? m_bcc->met : NULL);}

  GenInfo* GetGenInfo() {return (m_bcc ? m_bcc->genInfo : NULL);}
 
  std::vector<std::string>* GetTrigNames() {return (m_bcc ? m_bcc->triggerNames : NULL);}
  std::vector<bool>* GetTrigResults() {return (m_bcc ? m_bcc->triggerResults : NULL);}

  std::vector<std::string> GetTrigNamesActualRun() {return (m_bcc ? m_bcc->triggerNames_actualrun : std::vector<std::string>());}

  std::vector< ReconstructionHypothesis >* GetRecoHyps(){ return (m_bcc ? m_bcc->recoHyps : NULL);}
  
  /// scalar sum of the pt of all jets, leptons and missing transverse energy
  double GetHT();

  /// scalar sum of missing transverse energy and the pt of all leptons in the actual BaseCycleContainer
  double GetHTlep();

   /// scalar sum of the pt of all jets above certain pt-threshold within specified eta-range
  double GetHThad(double ptmin_jet, double etamax_jet);

  /// electron or muon in the actual BaseCycleContainer with largest transverse momentum
  Particle* GetPrimaryLepton();

  /// return ttbar generator information of the actual event
  TTbarGen* GetTTbarGen();

  /**
   * @short function to calculate the neutrino four-momentum from MET and charged lepton momenta
   *
   * Given the Decay:
   *
   * A -> B + Neutrino
   *
   * reconstruct the Neutrino pZ component of the Lorentz Vector given
   * it's Energy, px, and py are known.
   *
   * Calculation is carried with formula:
   *
   * P4A^2 = (P4B + P4Neutrino)^2
   *
   * assuming:
   *
   * P4Neutrino^2 = 0
   * P4B^2 = 0
   *
   * within the SM: m_neutrino = 0, mass of the second decay product
   * is neglected due to expected small mass, e.g. in case of the
   * electron: m_B = 0.5 MeV, for the muon: m_mu = 105 MeV and mass
   * of the W-boson: m_W = 80 GeV. m_B is used in formula in form:
   *
   * m_A^2 - m_B^2
   *
   * and therefore m_B can be neglected.
   *
   * The final equation is:
   *
   * (-pTlep^2) * x^2 + 2 * (mu * pZlep) * x + (mu^2 - Elep^2 * pTnu^2) = 0
   *
   * where
   * x is pz_nu
   * mu = mW^2 / 2 + pTlep * pTnu * cos(phi)
   * phi is angle between p_lepton and p_neutrino in transverse plane
   *
   */
  std::vector<LorentzVector> NeutrinoReconstruction(const LorentzVector lepton, const LorentzVector met);

  void FillHighMassTTbarHypotheses();

  /// print a list of all objects in the actual BaseCycleContainer
  void PrintEventContent();
  
  /// print a list of all generator level particles stored in the event
  /// it is possible to give the list a name, for example if the print
  /// function is called from different places in the code
  void PrintGenParticles(string name = "std");

  /// update the stored weight by multiplying it by the argument 'w'
  void ProduceWeight(double w);

  /// get the event weight
  double GetWeight();

  /// apply smearing of the tau energy to estimate uncertainty
  void ApplyTauEnergySmearing(double factor);

  /** @short integrated jet shape variable psi
   *
   * Psi = pT(rmin,rmax)/pT(rmin,R)
   *
   * R = jet radius, determined from jet type argument 
   * pT(r1,r2) = summed transverse momenta of all jet constituents between r1 and r2
   */
  double IntegratedJetShape(Jet* jet, double rmax, double rmin, E_JetType type );

  /// Jet charge, defined as summed charges of the PF constituents
  double JetCharge(Jet* jet);

  ///energy weighted jet charge, (1/Ejet * sum_j q_j*E_j^kappa), kappa to be chosen between ~0.2 and 1
  double EnergyWeightedJetCharge(Jet* jet, double kappa=1.);

  ///first jet moment <R^n>, where R_j is the distance of particle j to the jet axis
  double JetMoment(Jet* jet, int n=1);

  bool selection_passed(const identifier & selid) const{
      std::map<identifier, bool>::const_iterator it = selections.find(selid);
      if(it!=selections.end()) return it->second;
      throw std::runtime_error("did not find selection result for '" + selid.name() + "'");
  }
  
  void set_selection_passed(const identifier & selid, bool passed){
      selections[selid] = passed;
  }

 private:

  mutable SLogger m_logger;
  EventCalc();
  ~EventCalc();

  BaseCycleContainer* m_bcc;
  LuminosityHandler* m_lumi;

  // booleans to tell weather quantities have already been derived in an event
  bool b_Reconstruction;

  bool b_jetparticles;
  bool b_isoparticles;
  bool b_puisoparticles;

  // data members to store calculated results
  double m_TotalWeight;

  Particle* m_primlep;
  TTbarGen* m_ttgen;

  std::vector<PFParticle> m_jetparticles;
  std::vector<PFParticle> m_isoparticles;
  std::vector<PFParticle> m_puisoparticles;
  
  std::map<identifier, bool> selections;
};


#endif // EventCalc_H
