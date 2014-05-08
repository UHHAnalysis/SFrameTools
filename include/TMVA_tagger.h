#ifndef TMVA_tagger_H
#define TMVA_tagger_H

#include "SFrameTools/include/Objects.h"
#include "SFrameTools/include/fwd.h"
#include "SFrameTools/include/boost_includes.h" // for shared_array
#include <TMVA/Reader.h>
#include "TVector3.h"
#include <limits>
#include <algorithm>
#include <memory>
#include <TF1.h>
#include "Utils.h"
#include "EventCalc.h"


//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVectorXYZE;
class TMVA_tagger{
 private:
  static TMVA_tagger* m_instance;
  mutable SLogger m_logger;
 TMVA::Reader *reader=NULL;

 Float_t Nsubjets; //done
 Float_t Topjetmass; //done
 Float_t Subjet12mass; //done
 Float_t Subjet13mass; //done
 Float_t Subjet23mass; //done
 Float_t t3t2; //done
 Float_t t2t1; //done
 Float_t Q_valatility; //done
  Float_t npv; 
  Float_t Subjet12mass_pruned; //done
  Float_t Subjet13mass_pruned; //done
  Float_t Subjet23mass_pruned; //done
  Float_t t3t2p;
  Float_t t2t1p;
  Float_t subjet1pt; //sjpt1 done
  Float_t subjet2pt;//done	 
  Float_t subjet3pt; //done
  Float_t tau1; //done
  Float_t tau2; //done
  Float_t tau3; //done
  Float_t tau4;
  Float_t Topjetmass_pruned; //done
  Float_t tau1pruned; //done t1p
  Float_t tau2pruned; //done t2p
  Float_t tau3pruned; //done t3p
  Float_t HelAng12; //done HA12
  Float_t HelAng13; //done HA13
  Float_t HelAng23; //done HA23
  Float_t P02;
  Float_t P04;
  Float_t P06;
  Float_t P08;
  Float_t P10;
  Float_t w;
  Float_t st1;
  Float_t st2;
  Float_t st3;
 public:
TMVA_tagger();
  ~TMVA_tagger();

 
static TMVA_tagger* Instance();


 void Set_Reader(TString object);
 void push_variables(TopJet topjet);
 void reset_variables();
 bool IsTobiasTagged(TopJet topjet);
 bool IsTagged(TString taggername, TopJet topjet, double cutvalue, double& mva_value);
 Double_t GetMVA_value();
};

#endif
