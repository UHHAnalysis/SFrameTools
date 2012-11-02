#include "include/MCDataScaleFactors.h"

LeptonScaleFactors::LeptonScaleFactors(std::vector<std::string> correctionlist, E_SystShift syst_shift){

  m_syst_shift = syst_shift;

  if(correctionlist.size()%2!=0){
    std::cerr<< "not a valid list of correction factors given to LeptonScaleFactors" <<std::endl;
  }

  for(unsigned int i=0; i< correctionlist.size()/2; ++i){
    std::pair<std::string, double> correction (correctionlist[2*i], atof(correctionlist[2*i+1].c_str()));
    //std::cout << "Apply correction " << correction.first << " with factor " << correction.second <<std::endl;
    m_correctionlist.push_back(correction);
  }
  

}


double LeptonScaleFactors::GetWeight(){

  double triggerfactor = 0;
  double IDfactor = 0;
  double isofactor=0;
  double sumofweights=0;

  double Mu40triggerRunA[3] = {0.9799, 0.9621, 0.9851};
  double Mu40triggerRunB[3] = {0.9773, 0.9573, 0.9754};
  double Mu40triggerRunC[3] = {0.9817, 0.9640, 0.9973};

  double TightID_RunAB[3] = {0.9941, 0.9917, 0.9982};
  double TightID_RunC[3] = {0.9934, 0.9903, 0.9979};

  EventCalc* calc = EventCalc::Instance();

  //Take lepton with highest transverse momentum as reference. Make sure that routine is called after selection of exactly one good lepton.
  Particle* primlep = calc->GetPrimaryLepton();

  double eta = fabs(primlep->eta());

  int etabin=0;
  if(eta<0.9) etabin=0;
  else if(eta>=0.9 && eta<1.2) etabin=1;
  else if(eta>=1.2) etabin=2;


  for(unsigned int i=0; i<m_correctionlist.size(); ++i){

    double weight = m_correctionlist[i].second;

    //non isolated muons
    if(m_correctionlist[i].first == "MuonRunA"){
      triggerfactor += weight*Mu40triggerRunA[etabin];
      isofactor+=weight;
      IDfactor+= weight*TightID_RunAB[etabin];
    }
    else if (m_correctionlist[i].first == "MuonRunB"){
      triggerfactor += weight*Mu40triggerRunB[etabin];
      isofactor+=weight;
      IDfactor+= weight*TightID_RunAB[etabin];
    }
    else if (m_correctionlist[i].first == "MuonRunC"){
      triggerfactor += weight*Mu40triggerRunC[etabin];
      isofactor+=weight;
      IDfactor+= weight*TightID_RunC[etabin];
    }
    //isolated muons
    else if(m_correctionlist[i].first == "IsoMuonRunA"){
          
    }
    else if (m_correctionlist[i].first == "IsoMuonRunB"){

    }
    else if (m_correctionlist[i].first == "IsoMuonRunC"){

    }

    else{
      std::cerr<< "No information found for lepton correction named " << m_correctionlist[i].first <<std::endl;
    }

    sumofweights += m_correctionlist[i].second;

  }

  triggerfactor/=sumofweights;
  IDfactor/=sumofweights;
  isofactor/=sumofweights;

  return triggerfactor*IDfactor*isofactor;
}
