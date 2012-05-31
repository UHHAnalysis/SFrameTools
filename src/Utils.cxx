#include "../include/Utils.h"


bool TopTag(TopJet topjet,  double &mjet, int &nsubjets, double &mmin){
 
  nsubjets=topjet.numberOfDaughters();

  LorentzVector allsubjets(0,0,0,0);
  
  for(int j=0; j<topjet.numberOfDaughters(); ++j){
    allsubjets += topjet.subjets()[j].v4();
  }
  if(!allsubjets.isTimelike()){
    mjet=0;
    mmin=0;
    return false;
  }

  mjet = allsubjets.M();
  
  if(nsubjets>=3) {

    std::vector<Particle> subjets = topjet.subjets();
    sort(subjets.begin(), subjets.end(), HigherPt());

    double m01 = 0;
    if( (subjets[0].v4()+subjets[1].v4()).isTimelike()) 
      m01=(subjets[0].v4()+subjets[1].v4()).M();
    double m02 = 0;
    if( (subjets[0].v4()+subjets[2].v4()).isTimelike() ) 
      m02=(subjets[0].v4()+subjets[2].v4()).M();
    double m12 = 0;
    if( (subjets[1].v4()+subjets[2].v4()).isTimelike()  )
      m12 = (subjets[1].v4()+subjets[2].v4()).M();
    
    //minimum pairwise mass > 50 GeV/c^2
    mmin = std::min(m01,std::min(m02,m12));
  }

  //at least 3 sub-jets
  if(nsubjets<3) return false;
  //minimum pairwise mass > 50 GeV/c^2
  if(mmin<50) return false;
  //jet mass between 140 and 250 GeV/c^2
  if(mjet<140 || mjet>250) return false;
  
  return true;
}

double HTlep(const BaseCycleContainer *bcc){

  double htlep=0;

  if(bcc->electrons){
    for(unsigned int i=0; i<bcc->electrons->size(); ++i){
      htlep += bcc->electrons->at(i).pt();
    }
  }
  if(bcc->muons){
    for(unsigned int i=0; i<bcc->muons->size(); ++i){
      htlep += bcc->muons->at(i).pt();
    }
  }
  if(bcc->met) htlep += bcc->met->pt();

  return htlep;

}

Jet* nextJet(const Particle *p, std::vector<Jet> *jets){

  double deltarmin = double_infinity();
  Jet* nextjet=0;
  for(unsigned int i=0; i<jets->size();++i){
    if(jets->at(i).deltaR(*p) < deltarmin){
      deltarmin = jets->at(i).deltaR(*p);
      nextjet = &jets->at(i);
    }
  }

  return nextjet;
}

double pTrel(const Particle *p, std::vector<Jet> *jets){

  double ptrel=0;

  Jet* nextjet =  nextJet(p,jets);

  TVector3 p3(p->v4().Px(),p->v4().Py(),p->v4().Pz());
  TVector3 jet3(nextjet->v4().Px(),nextjet->v4().Py(),nextjet->v4().Pz());

  if(p3.Mag()!=0 && jet3.Mag()!=0){
    double sin_alpha = (p3.Cross(jet3)).Mag()/p3.Mag()/jet3.Mag();
    ptrel = p3.Mag()*sin_alpha;
  }
  else{
    std::cout << "something strange happend in the ptrel calculation: either lepton or jet momentum is 0" <<std::endl;
  }

  return ptrel;
}

double deltaRmin(const Particle *p, std::vector<Jet> *jets){
  return nextJet(p,jets)->deltaR(*p);
}


double double_infinity(){
  return std::numeric_limits<double>::infinity() ;
}

int int_infinity(){
  return std::numeric_limits<int>::max() ;
}
