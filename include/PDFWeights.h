#ifndef PDFWeights_H
#define PDFWeights_H

#include "SFrameTools/include/EventCalc.h"

#include "LHAPDF/LHAPDF.h"
#include "TSystem.h"

/**
 *  @short pdf re-weighting tool
 *
 *
 */
class PDFWeights {
public:

  ///Default constructor
  PDFWeights(E_SystShift syst_shift, TString pdfname= "cteq66", TString pdfweightdir="");
  ///Default destructor
  ~PDFWeights() {};

  ///returns quadratically summed pdf uncertainties
  double GetWeight(unsigned int index);

  ///returns list of pdf uncertainties from eigenvectors
  std::vector<double> GetWeightList();

  unsigned int GetNWeights(){return m_N_unc;}
  
private:

  mutable SLogger m_logger;

  bool m_libvalid;
  E_SystShift m_syst_shift;

  unsigned int m_N_unc;

  bool m_normalize_to_total_sum;
  std::vector<double> m_sumofweights;
  int m_N_tot;
  

};

#endif
