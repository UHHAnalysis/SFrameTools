#include "include/TriggerWeight.h"
#include "SFrameTools/include/BaseCycleContainer.h"
#include "SFrameTools/include/EventCalc.h"

TriggerWeight::TriggerWeight(TString filename_mc, TString mode){

  TFile *file_mc = new TFile(filename_mc);

  h_mc=(TH1F*) file_mc->Get("triggerSF");
  
  m_mode=mode;

}

double TriggerWeight::produceWeight(BaseCycleContainer* bcc){

  double HTSubJets = 0.;

  std::vector<Particle> subjets_top;

  for (unsigned int itj = 0; itj < bcc->topjets->size(); itj++){
   if(bcc->topjets->at(itj).pt()<150.) continue;
   TopJet topjet = bcc->topjets->at(itj);
   subjets_top=topjet.subjets();
   for (unsigned int subj = 0; subj < subjets_top.size(); subj++){
     HTSubJets += subjets_top.at(subj).pt();
   }
  }
 
  if(HTSubJets>800.){
    return 1.;
  }
  else{
    int bin=h_mc->GetXaxis()->FindBin(HTSubJets);
    double weight=h_mc->GetBinContent(bin);
    if(m_mode=="up"){
      weight=weight*1.5;
    }
    if(m_mode=="down"){
      weight=weight*0.5;
    }

    return weight;
  }

}
