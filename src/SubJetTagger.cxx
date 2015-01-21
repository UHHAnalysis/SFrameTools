#include "SubJetTagger.h"
#include "include/Utils.h"

using namespace std;

SubJetTagger::~SubJetTagger(){}

CMSTopTagger::CMSTopTagger(double mminLower, double mjetLower, double mjetUpper, double tau32Upper) {
  m_mjet = 0;
  m_nsubjets = 0;
  m_mmin = 0;
  m_mminLower = mminLower;
  m_mjetLower = mjetLower; 
  m_mjetUpper= mjetUpper;
  m_tau32Upper = tau32Upper;
} 

CMSTopTagger::~CMSTopTagger(){}


bool CMSTopTagger::Tag(const TopJet& topjet){
  m_nsubjets=topjet.numberOfDaughters();
  m_mjet = m_mmin = 0.0;
  LorentzVector allsubjets(0,0,0,0);
  for(int j=0; j<topjet.numberOfDaughters(); ++j) {
    allsubjets += topjet.subjets()[j].v4();
  }
  if(!allsubjets.isTimelike()) {
    return false;
  }
  m_mjet = allsubjets.M();
  
  if(m_nsubjets>=3) {
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
    m_mmin = std::min(m01,std::min(m02,m12));
  }
  //at least 3 sub-jets
  if(m_nsubjets<3) return false;
  //minimum pairwise mass > 50 GeV/c^2 or what you diceded
  if(m_mmin<m_mminLower) return false;
  //jet mass between 140 and 250 GeV/c^2 -""-
  if(m_mjet<m_mjetLower || m_mjet>m_mjetUpper) return false;
  if(std::isfinite(m_tau32Upper)){
      if(topjet.tau2() == 0) return false;
      if(topjet.tau3() / topjet.tau2() > m_tau32Upper) return false;
  }
  return true;
}  

void CMSTopTagger::SetTau32Cut(double cut){
  m_tau32Upper=cut;
}


map<string,double> CMSTopTagger::TagVar(){
  map<string,double> mymap;
  mymap.insert(pair<string,double>("mjet",m_mjet));
  mymap.insert(pair<string,double>("nsubjets",m_nsubjets));
  mymap.insert(pair<string,double>("mmin",m_mmin));
  return mymap;
}
