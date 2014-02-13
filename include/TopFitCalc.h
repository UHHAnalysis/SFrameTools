// Dear emacs, this is -*- c++ -*-
#ifndef TopFitCalc_H
#define TopFitCalc_H

#include "EventCalc.h"
#include "TVector3.h"
#include "Utils.h"
#include "TTbarGen.h"
#include "TMinuit.h"

/**
 * class for the calculation of Tops with hadronic Tags
 *
 * author Daniel Gonzalez
 */

class TopFitCalc
{
 private:
  static TopFitCalc* m_instance;

 public:
  static TopFitCalc* Instance();

  /*
   * Take the C/A-Jet with the highest pT and use it as hadronic Top. 
   * For the leptonic side the primary lepton is taken. A list is generated 
   * with all possible combinations of one ak 0.5 Jets. Later one ak 0.5 Jet is selected  
   * via chi_square.
   * ak 0.5 Jets with  
   */
  void CalculateSelection();

  /*
   * Take the C/A-Jet with the TopTag and use it as hadronic Top.
   * For the leptonic side the primary lepton is taken. A list is generated 
   * with all possible combinations of single ak5 jets. Later one hypothesis is selected 
   * via chi_square. Make sure a reasonable selection is run. There is hardcoded a minimum
   * distance of lepton and TopTag of 0.8 and between ak5 jet and TopTag of 1.3 (measure in delta R)
   *
   */
  void CalculateTopTag();

  /*
   * @short Reconstruction of ttbar system at high mass
   *
   * Reconstructs the ttbar system by constructing a chi2 function
   * taking into account different combinations of jets with the lepton and MET. 
   * The solution which minimizes chi2 in an event is used, the hadronic side
   * is allowed to consist of 1,2 or 3 jets to account for very boosted topologies.
   */
  void FillHighMassTTbarHypotheses();


  void Reset();

  /*
   * @short Extension of EventCalc::NeutrinoReconstruction
   *
   * Extension of the neutrino reconstruction in EventCalc. This function 
   * performs a fit if the sqrt() part of the solution becomes imaginary.
   *
   * Basic idea is that for the case the Neutrino becomes an Im-Part you do a fit
   * on the px and py components. If you look into the kinematics there is a way of 
   * doing this with a W-mass constraint. As the name indicates polarcoordinates 
   * are used in the computation.
   * 
   */
  std::vector<LorentzVector> NeutrinoFitPolar(const LorentzVector lepton, const LorentzVector met);

 private:

  mutable SLogger m_logger;
  TopFitCalc();
  ~TopFitCalc();


  TMinuit* positiv;// For the Neutrino-fit

  // data members to store calculated results
  bool b_Reconstruction;

};



#endif // TopFitCalc
