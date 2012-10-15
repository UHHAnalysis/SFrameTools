#include "../include/Utils.h"


//variable HEP Tagger from Rebekka

bool variableHepTopTag(TopJet topjet, double ptJetMin, double massWindowLower, double massWindowUpper, double cutCondition2, double cutCondition3){
 
  double mjet;
  double ptjet;
  int nsubjets;

  double topmass=172.3;
  double wmass=80.4;

  nsubjets=topjet.numberOfDaughters();

  LorentzVector allsubjets(0,0,0,0);
  
  for(int j=0; j<topjet.numberOfDaughters(); ++j){
    allsubjets += topjet.subjets()[j].v4();
  }
  if(!allsubjets.isTimelike()){
    mjet=0;
    return false;
  }

  mjet = allsubjets.M();
  ptjet= allsubjets.Pt();

  double m12, m13, m23;

  //The subjetcs have to be three
  if(nsubjets==3) {

    std::vector<Particle> subjets = topjet.subjets();
    sort(subjets.begin(), subjets.end(), HigherPt());

    m12 = 0;
    if( (subjets[0].v4()+subjets[1].v4()).isTimelike()) 
      m12=(subjets[0].v4()+subjets[1].v4()).M();
    m13 = 0;
    if( (subjets[0].v4()+subjets[2].v4()).isTimelike() ) 
      m13=(subjets[0].v4()+subjets[2].v4()).M();
    m23 = 0;
    if( (subjets[1].v4()+subjets[2].v4()).isTimelike()  )
      m23 = (subjets[1].v4()+subjets[2].v4()).M();
   
  }
  else{
    return false;
  }

  double rmin=massWindowLower*wmass/topmass;
  double rmax=massWindowUpper*wmass/topmass;

  int keep=0;

  //Conditions on the subjects: at least one has to be true
  //1 condition
  if(atan(m13/m12)>0.2 && atan(m13/m12)<1.3 && m23/mjet>rmin && m23/mjet<rmax) keep=1;
  
  //2 condition
  double cond2left=pow(rmin,2)*(1+pow((m13/m12),2));
  double cond2cent=1-pow(m23/mjet,2);
  double cond2right=pow(rmax,2)*(1+pow(m13/m12,2));

  if(cond2left<cond2cent && cond2cent<cond2right && m23/mjet>cutCondition2) keep=1;

  //3 condition
  double cond3left=pow(rmin,2)*(1+pow((m12/m13),2));
  double cond3cent=1-pow(m23/mjet,2);
  double cond3right=pow(rmax,2)*(1+pow(m12/m13,2));

  if(cond3left<cond3cent && cond3cent<cond3right && m23/mjet>cutCondition3) keep=1;

  //Final requirement: at least one of the three subjets conditions and total pt
  if(keep==1 && ptjet>ptJetMin){
    return true;
  }
  else{
    return false;
  }

}






//HEP Tagger from Ivan

bool HepTopTag(TopJet topjet){
 
  double mjet;
  double ptjet;
  int nsubjets;

  double topmass=172.3;
  double wmass=80.4;

  nsubjets=topjet.numberOfDaughters();

  LorentzVector allsubjets(0,0,0,0);
  
  for(int j=0; j<topjet.numberOfDaughters(); ++j){
    allsubjets += topjet.subjets()[j].v4();
  }
  if(!allsubjets.isTimelike()){
    mjet=0;
    return false;
  }

  mjet = allsubjets.M();
  ptjet= allsubjets.Pt();

  double m12, m13, m23;

  //The subjetcs have to be three
  if(nsubjets==3) {

    std::vector<Particle> subjets = topjet.subjets();
    sort(subjets.begin(), subjets.end(), HigherPt());

    m12 = 0;
    if( (subjets[0].v4()+subjets[1].v4()).isTimelike()) 
      m12=(subjets[0].v4()+subjets[1].v4()).M();
    m13 = 0;
    if( (subjets[0].v4()+subjets[2].v4()).isTimelike() ) 
      m13=(subjets[0].v4()+subjets[2].v4()).M();
    m23 = 0;
    if( (subjets[1].v4()+subjets[2].v4()).isTimelike()  )
      m23 = (subjets[1].v4()+subjets[2].v4()).M();
   
  }
  else{
    return false;
  }

  double rmin=0.85*wmass/topmass;
  double rmax=1.15*wmass/topmass;

  int keep=0;

  //Conditions on the subjects: at least one has to be true
  //1 condition
  if(atan(m13/m12)>0.2 && atan(m13/m12)<1.3 && m23/mjet>rmin && m23/mjet<rmax) keep=1;
  
  //2 condition
  double cond2left=pow(rmin,2)*(1+pow((m13/m12),2));
  double cond2cent=1-pow(m23/mjet,2);
  double cond2right=pow(rmax,2)*(1+pow(m13/m12,2));

  if(cond2left<cond2cent && cond2cent<cond2right && m23/mjet>0.35) keep=1;

  //3 condition
  double cond3left=pow(rmin,2)*(1+pow((m12/m13),2));
  double cond3cent=1-pow(m23/mjet,2);
  double cond3right=pow(rmax,2)*(1+pow(m12/m13,2));

  if(cond3left<cond3cent && cond3cent<cond3right && m23/mjet>0.35) keep=1;

  //Final requirement: at least one of the three subjets conditions and total pt
  if(keep==1 && ptjet>200){
    return true;
  }
  else{
    return false;
  }

}


//default values (mminLower=50., mjetLower=140, mjetUpper=250.) defined in Utils.h
 bool variableTopTag(TopJet topjet, double &mjet, int &nsubjets, double &mmin, double mminLower, double mjetLower, double mjetUpper){
 
  nsubjets=topjet.numberOfDaughters();
 
  LorentzVector allsubjets(0,0,0,0);
 
  for(int j=0; j<topjet.numberOfDaughters(); ++j){
  allsubjets += topjet.subjets()[j].v4();
  }
  if(!allsubjets.isTimelike()){
  mjet=0;
  mmin=0;
 //  mminLower=50;
//   mjetLower=140;
//   mjetUpper=250;
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
  if( (subjets[1].v4()+subjets[2].v4()).isTimelike() )
  m12 = (subjets[1].v4()+subjets[2].v4()).M();
 
  //minimum pairwise mass
  mmin = std::min(m01,std::min(m02,m12));
  }
 
  //at least 3 sub-jets
  if(nsubjets<3) return false;
  //minimum pairwise mass > 50 GeV/c^2
  if(mmin<mminLower) return false;
  //jet mass between 140 and 250 GeV/c^2
  if(mjet<mjetLower || mjet>mjetUpper) return false;
 
  return true;
 }


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
    
    //minimum pairwise mass 
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

bool WTag(TopJet prunedjet,  double& mjet, int &nsubjets, double& massdrop){

  nsubjets=prunedjet.numberOfDaughters();

  mjet = 0;
  if(prunedjet.v4().isTimelike())
    mjet = prunedjet.v4().M(); 

  //calculate mass drop for first sub-jet ordered in pt
  massdrop = 0;
  if(nsubjets>=1 && mjet>0){

    std::vector< Particle > subjets = prunedjet.subjets();
    sort(subjets.begin(), subjets.end(), HigherPt());

    double m1 = 0;
    if(subjets[0].v4().isTimelike())
      m1 = subjets[0].v4().M();

    massdrop = m1/mjet;
  }
  
  //at least 2 sub-jets
  if(nsubjets<2) return false;
  //60 GeV < pruned jet mass < 100 GeV
  if(mjet <= 60 || mjet >= 100) return false;
  //mass drop < 0.4
  if(massdrop>=0.4) return false;

  return true;

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

TVector3 toVector(LorentzVector v4){

  TVector3 v3(0,0,0);
  v3.SetX(v4.X());
  v3.SetY(v4.Y());
  v3.SetZ(v4.Z());
  return v3;
}

TVector3 toVector(LorentzVectorXYZE v4){

  TVector3 v3(0,0,0);
  v3.SetX(v4.X());
  v3.SetY(v4.Y());
  v3.SetZ(v4.Z());
  return v3;
}

LorentzVectorXYZE toXYZ(LorentzVector v4){

  LorentzVectorXYZE v4_new(0,0,0,0);
  v4_new.SetPx(v4.Px());
  v4_new.SetPy(v4.Py());
  v4_new.SetPz(v4.Pz());
  v4_new.SetE(v4.E());
  return v4_new;
}

LorentzVector toPtEtaPhi(LorentzVectorXYZE v4){

  LorentzVector v4_new(0,0,0,0);
  v4_new.SetPt(v4.Pt());
  v4_new.SetEta(v4.Eta());
  v4_new.SetPhi(v4.Phi());
  v4_new.SetE(v4.E());
  return v4_new;
}

double deltaR(LorentzVector v1, LorentzVector v2){

  Particle p1;
  p1.set_v4(v1);
  Particle p2;
  p2.set_v4(v2);
  return p1.deltaR(p2);
}

double double_infinity(){
  return std::numeric_limits<double>::infinity() ;
}

int int_infinity(){
  return std::numeric_limits<int>::max() ;
}

int myPow(int x, unsigned int p) {
  int i = 1;
  for (unsigned int j = 1; j <= p; j++)  i *= x;
  return i;
}
