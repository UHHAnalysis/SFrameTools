#include "include/Utils.h"
#include "NtupleWriter/include/JetProps.h"

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
namespace external {
#include "include/HEPTopTagger.h"
}



//subjet b-tagger, returns number of b-tagged subjets

int subJetBTag(TopJet topjet, E_BtagType type){
  int nBTagsSub = 0;
  float discriminator_cut;
  if(type==e_CSVL) discriminator_cut = 0.244;
  if(type==e_CSVM) discriminator_cut = 0.679;
  if(type==e_CSVT) discriminator_cut = 0.898;

  //Create a vector of subjets
  std::vector<Particle> subjets_top;
  //Create a float vector of the subjets discriminators
  std::vector<float> btagsub_combinedSecondaryVertex_top;

  //Fill the vector of subjets with the subjets of a topjet
  subjets_top=topjet.subjets();
  //Fill the vector of discriminators with the discriminators of the subjets of a certain topjet
  btagsub_combinedSecondaryVertex_top=topjet.btagsub_combinedSecondaryVertex();

  //Looping over subjets and checking if they are b-tagged
  for(unsigned int i=0; i < btagsub_combinedSecondaryVertex_top.size(); ++i){
    float test=btagsub_combinedSecondaryVertex_top[i];
    if(test>discriminator_cut){
      nBTagsSub += 1;
      //This means it is b-tagged
    }
  }
  return nBTagsSub;
}


bool HiggsTag(TopJet topjet, E_BtagType type1, E_BtagType type2){
  int nBTagsSub1 = 0;
  int nBTagsSub2 = 0;
  float discriminator_cut1;
  float discriminator_cut2;
  if(type1==e_CSVL) discriminator_cut1 = 0.244;
  if(type1==e_CSVM) discriminator_cut1 = 0.679;
  if(type1==e_CSVT) discriminator_cut1 = 0.898;
  if(type2==e_CSVL) discriminator_cut2 = 0.244;
  if(type2==e_CSVM) discriminator_cut2 = 0.679;
  if(type2==e_CSVT) discriminator_cut2 = 0.898;
  // std::cout << "discriminator_cut1: " <<  discriminator_cut1 << " discriminator_cut2: "<< discriminator_cut1 << std::endl;
  //Create a vector of subjets
  std::vector<Particle> subjets_top;
  //Create a float vector of the subjets discriminators
  std::vector<float> btagsub_combinedSecondaryVertex_top;

  //Fill the vector of subjets with the subjets of a topjet
  subjets_top=topjet.subjets();
  //Fill the vector of discriminators with the discriminators of the subjets of a certain topjet
  btagsub_combinedSecondaryVertex_top=topjet.btagsub_combinedSecondaryVertex();

  //Looping over subjets and checking if they are b-tagged
  for(unsigned int i=0; i < btagsub_combinedSecondaryVertex_top.size(); ++i){
    float test=btagsub_combinedSecondaryVertex_top[i];
    if (nBTagsSub1 != 0 && test>discriminator_cut2) nBTagsSub2 =+ 1;
    if(test>discriminator_cut1){
      if(test>discriminator_cut2 && nBTagsSub2==0) nBTagsSub2+=1;
      else nBTagsSub1 += 1;      //This means it is b-tagged  
 
    }
  }
  if (nBTagsSub1!=0 && nBTagsSub2!=0) return true;
  else return false;
}


bool HepTopTagFull(TopJet topjet, std::vector<PFParticle>* allparts){

  //Transform the SFrame TopJet object in a fastjet::PseudoJet

  if (!allparts) return false;

  fastjet::ClusterSequence* JetFinder;
  fastjet::JetDefinition* JetDef ;

  JetProps jp(&topjet, allparts);
  std::vector<fastjet::PseudoJet> jetpart = jp.GetJetConstituents();
 
  //Clustering definition
  double conesize=3;
  JetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm,conesize); 

  JetFinder = new fastjet::ClusterSequence(jetpart, *JetDef);

  std::vector<fastjet::PseudoJet> tops = JetFinder->inclusive_jets(10.);

  if (tops.size() != 1){
    std::cout << "Problem! Doesn't give exactly one jet!!!!!!!!!!!!!!Gives " << tops.size() << " jets" << std::endl;
    delete JetFinder;
    delete JetDef;
    return false;
  }

  std::vector<fastjet::PseudoJet> SortedJets = sorted_by_pt(tops);


  //Run the HEPTopTagger
  external::HEPTopTagger tagger(*JetFinder, SortedJets[0]);

  //Mass window to be applied in a second step
  tagger.set_top_range(0.0, 10000.0);
  tagger.set_mass_drop_threshold(0.8);
  tagger.set_max_subjet_mass(30);

  tagger.run_tagger();

  delete JetFinder;
  delete JetDef;
 
  if (tagger.is_masscut_passed()) return true;
  else return false;

  return true;
 
}



// global function to define a tagged jet

bool IsTagged(Jet & jet, E_BtagType type)
{
    if(type==e_CSVL && jet.btag_combinedSecondaryVertex()>0.244) return true;
    if(type==e_CSVM && jet.btag_combinedSecondaryVertex()>0.679) return true;
    if(type==e_CSVT && jet.btag_combinedSecondaryVertex()>0.898) return true;
    if(type==e_JPL && jet.btag_jetProbability()>0.275) return true;
    if(type==e_JPM && jet.btag_jetProbability()>0.545) return true;
    if(type==e_JPT && jet.btag_jetProbability()>0.790) return true;       
    
    return false;
}

//variable HEP Tagger from Rebekka

bool variableHepTopTag(TopJet topjet, double ptJetMin, double massWindowLower, double massWindowUpper, double cutCondition2, double cutCondition3)
{
    double mjet;
    double ptjet;
    int nsubjets;

    double topmass=172.3;
    double wmass=80.4;

    nsubjets=topjet.numberOfDaughters();

    LorentzVector allsubjets(0,0,0,0);

    for(int j=0; j<topjet.numberOfDaughters(); ++j) {
        allsubjets += topjet.subjets()[j].v4();
    }
    if(!allsubjets.isTimelike()) {
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

    } else {
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
    if(keep==1 && ptjet>ptJetMin) {
        return true;
    } else {
        return false;
    }

}


//HEP Tagger from Ivan

bool HepTopTag(TopJet topjet)
{
  //call variable tagger with default parameters
  return variableHepTopTag(topjet);

}

//default values (mminLower=50., mjetLower=140, mjetUpper=250.) defined in Utils.h
bool variableTopTag(TopJet topjet, double &mjet, int &nsubjets, double &mmin, double mminLower, double mjetLower, double mjetUpper)
{

    nsubjets=topjet.numberOfDaughters();

    LorentzVector allsubjets(0,0,0,0);

    for(int j=0; j<topjet.numberOfDaughters(); ++j) {
        allsubjets += topjet.subjets()[j].v4();
    }
    if(!allsubjets.isTimelike()) {
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


bool TopTag(TopJet topjet,  double &mjet, int &nsubjets, double &mmin)
{
  //call variable tagger with default parameters
  return variableTopTag(topjet, mjet, nsubjets, mmin);

}
 
Jet* nextJet(const Particle *p, std::vector<Jet> *jets)
{

    double deltarmin = double_infinity();
    Jet* nextjet=0;
    for(unsigned int i=0; i<jets->size(); ++i) {
        Jet ji = jets->at(i);
	if (fabs(p->pt() - ji.pt())<1e-8) continue; // skip identical particle
        if(jets->at(i).deltaR(*p) < deltarmin) {
            deltarmin = jets->at(i).deltaR(*p);
            nextjet = &jets->at(i);
        }
    }

    return nextjet;
}

bool WTag(TopJet prunedjet,  double& mjet, int &nsubjets, double& massdrop)
{

    nsubjets=prunedjet.numberOfDaughters();

    mjet = 0;
    if(prunedjet.v4().isTimelike())
        mjet = prunedjet.v4().M();

    //calculate mass drop for first sub-jet ordered in pt
    massdrop = 0;
    if(nsubjets>=1 && mjet>0) {

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

float relIsoMuon( Muon mu, float deltaR ){
  float chargedHadronIso=0;
  float neutralHadronIso=0;
  float photonIso=0;
  float puiso=0;

  EventCalc* calc = EventCalc::Instance();

  unsigned int npfp=calc->GetIsoPFParticles()->size();
  for(unsigned int j=0; j<npfp; ++j){
    PFParticle pfp = calc->GetIsoPFParticles()->at(j);
    if(pfp.deltaR(mu)<deltaR ){
      if(pfp.particleID() == PFParticle::eH && pfp.pt()>0.0 && pfp.deltaR(mu)>0.0001 ) chargedHadronIso += pfp.pt();
      if(pfp.particleID() == PFParticle::eH0 && pfp.pt()>0.5 && pfp.deltaR(mu)>0.01) neutralHadronIso += pfp.pt();
      if(pfp.particleID() == PFParticle::eGamma && pfp.pt()>0.5 && pfp.deltaR(mu)>0.01) photonIso += pfp.pt();
    }
  }
  
  unsigned int npfppu=calc->GetPUIsoPFParticles()->size();
  for(unsigned int j=0; j<npfppu; ++j){
    PFParticle pfp = calc->GetPUIsoPFParticles()->at(j);
    if(pfp.deltaR(mu)<deltaR ){
      if(pfp.particleID() == PFParticle::eH && pfp.pt()>0.5 && pfp.deltaR(mu)>0.01 ) puiso += pfp.pt();
    }
  }
  
  return (chargedHadronIso + std::max( 0.0, neutralHadronIso + photonIso - 0.5*puiso))/mu.pt();
}

double pTrel(const Particle *p, std::vector<Jet> *jets)
{

    double ptrel=0;
    Jet* nextjet =  nextJet(p,jets);
    if (!nextjet) return ptrel;

    TVector3 p3(p->v4().Px(),p->v4().Py(),p->v4().Pz());
    TVector3 jet3(nextjet->v4().Px(),nextjet->v4().Py(),nextjet->v4().Pz());

    if(p3.Mag()!=0 && jet3.Mag()!=0) {
        double sin_alpha = (p3.Cross(jet3)).Mag()/p3.Mag()/jet3.Mag();
        ptrel = p3.Mag()*sin_alpha;
    } else {
        std::cout << "something strange happend in the ptrel calculation: either lepton or jet momentum is 0" <<std::endl;
    }

    return ptrel;
}

double deltaRmin(const Particle *p, std::vector<Jet> *jets)
{
    Jet* j = nextJet(p,jets);
    double dr = 999.;
    if (j) dr = j->deltaR(*p);
    return dr;
}

TVector3 toVector(LorentzVector v4)
{

    TVector3 v3(0,0,0);
    v3.SetX(v4.X());
    v3.SetY(v4.Y());
    v3.SetZ(v4.Z());
    return v3;
}

TVector3 toVector(LorentzVectorXYZE v4)
{

    TVector3 v3(0,0,0);
    v3.SetX(v4.X());
    v3.SetY(v4.Y());
    v3.SetZ(v4.Z());
    return v3;
}

LorentzVectorXYZE toXYZ(LorentzVector v4)
{

    LorentzVectorXYZE v4_new(0,0,0,0);
    v4_new.SetPx(v4.Px());
    v4_new.SetPy(v4.Py());
    v4_new.SetPz(v4.Pz());
    v4_new.SetE(v4.E());
    return v4_new;
}

LorentzVector toPtEtaPhi(LorentzVectorXYZE v4)
{

    LorentzVector v4_new(0,0,0,0);
    v4_new.SetPt(v4.Pt());
    v4_new.SetEta(v4.Eta());
    v4_new.SetPhi(v4.Phi());
    v4_new.SetE(v4.E());
    return v4_new;
}

double deltaR(LorentzVector v1, LorentzVector v2)
{

    Particle p1;
    p1.set_v4(v1);
    Particle p2;
    p2.set_v4(v2);
    return p1.deltaR(p2);
}

double double_infinity()
{
    return std::numeric_limits<double>::infinity() ;
}

int int_infinity()
{
    return std::numeric_limits<int>::max() ;
}

int myPow(int x, unsigned int p)
{
    int i = 1;
    for (unsigned int j = 1; j <= p; j++)  i *= x;
    return i;
}

