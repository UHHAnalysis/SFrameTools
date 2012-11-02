#ifndef MCDataScaleFactors_H
#define MCDataScaleFactors_H

#include "include/Utils.h"
#include "include/EventCalc.h"

/**
 *  @short module to apply data-MC lepton scale factors for trigger and ID
 *
 *  
 */
class LeptonScaleFactors{
 public:
 /**
  * constructor
  *
  * first argument: list of corrections to be applied together with weight factors for luminosity, example: "MuonRunA 1.5 MuonRunB 2.6 MuonRunC 7.8"
  *
  * second argument: systematic shift
  * @see E_SystShift
  */
  LeptonScaleFactors(std::vector<std::string> correctionlist, E_SystShift syst_shift=e_Default);
  ///Default destructor
  ~LeptonScaleFactors(){};

  ///return the weighted correction factor
  double GetWeight();

 private:
  E_SystShift m_syst_shift;
  std::vector<std::pair<std::string, double> > m_correctionlist;

};


#endif
