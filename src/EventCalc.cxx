// Dear emacs, this is -*- c++ -*-

#include "include/EventCalc.h"
#include <iostream>

using namespace std;

EventCalc* EventCalc::m_instance = NULL;

EventCalc* EventCalc::Instance()
{
  // Get a pointer to the object handler.
  // This is the only way to access this class, 
  // since it's a singleton. This method is accessible
  // from everywhere.

  if (m_instance == NULL){
    m_instance = new EventCalc();
  }
  return m_instance;    
}

EventCalc::EventCalc() : m_logger( "EventCalc" )
{
  // constructor: initialise all variables
  m_logger << DEBUG << "Constructor called." << SLogger::endmsg;
  m_bcc = NULL;
  m_lumi = NULL;
  m_primlep = NULL;
  m_ttgen = NULL;
}

EventCalc::~EventCalc()
{
  // default destructor
}

void EventCalc::Reset()
{
  // reset: set all booleans to false
  // this has to be done at the beginning of each event
  // after a call to reset all quantities will be re-calculated 
  // when they are accessed
  
  // (re-)set the pointers using the ObjectHandler

  //m_bcc = ?
  //m_lumi = ?

  // reset booleans
  b_Reconstruction = false;
  b_jetparticles = false;
  b_isoparticles = false;
  b_puisoparticles = false;
  selections.clear();


  m_TotalWeight = 1.;
  
  m_primlep = NULL;
  delete m_ttgen;
  m_ttgen = NULL;
}

void EventCalc::SetBaseCycleContainer(BaseCycleContainer* bcc)
{
  // set the internal pointer to the container with all objects
  m_bcc = bcc;
  m_logger << DEBUG << "Pointer to BaseCycleContainer set." << SLogger::endmsg;
}

void EventCalc::SetLumiHandler(LuminosityHandler* lh)
{
  // set the internal pointer to the container with all objects
  m_lumi = lh;
  m_logger << DEBUG << "Pointer to LumiHandler set." << SLogger::endmsg;
}

BaseCycleContainer* EventCalc::GetBaseCycleContainer()
{
  // return the pointer to the container with all objects
  if (!m_bcc){
    m_logger << WARNING << "Pointer to BaseCycleContainer is NULL." << SLogger::endmsg;
  }
  return m_bcc;
}

LuminosityHandler* EventCalc::GetLumiHandler()
{
  // return the pointer to the container with all objects
  if (!m_lumi){
    m_logger << DEBUG << "Pointer to LumiHandler is NULL." << SLogger::endmsg;
  }
  return m_lumi;
}

double EventCalc::GetHT()
{
    // calculate HT, which is defined as the scalar sum of all
    // jets, leptons and missing transverse momentum in the event
    double m_HT = GetHTlep();
    // sum over pt of all jets
    if(m_bcc->jets){
        for(vector<Jet>::const_iterator jet = m_bcc->jets->begin(); jet != m_bcc->jets->end(); ++jet){
            m_HT += jet->pt();
        }
    }
    return m_HT;
}

double EventCalc::GetHTlep()
{
    // calculate HT_lep, which is defined as the scalar sum of all
    // leptons and missing transverse momentum in the event
    double m_HTlep=0;

    // sum over pt of all electrons
    if(m_bcc->electrons){
        for(vector<Electron>::const_iterator ele = m_bcc->electrons->begin(); ele != m_bcc->electrons->end(); ++ele){
            m_HTlep += ele->pt();
        }
    }

    // sum over pt of all muons
    if(m_bcc->muons){
        for(vector<Muon>::const_iterator mu = m_bcc->muons->begin(); mu != m_bcc->muons->end(); ++mu){
            m_HTlep += mu->pt();
        }
    }

    // sum over pt of all taus
    if(m_bcc->taus){
        for(vector<Tau>::const_iterator tau = m_bcc->taus->begin(); tau != m_bcc->taus->end(); ++tau){
            m_HTlep += tau->pt();
        }
    }

    if(m_bcc->met) m_HTlep += m_bcc->met->pt();   
    return m_HTlep;
}

Particle* EventCalc::GetPrimaryLepton(){

  if(!m_primlep){
    double ptmax = -999;
    for(unsigned int i=0; i<m_bcc->electrons->size(); ++i){
      if(m_bcc->electrons->at(i).pt()>ptmax){
	ptmax = m_bcc->electrons->at(i).pt();
	m_primlep = &m_bcc->electrons->at(i);
      }
    }
    for(unsigned int i=0; i<m_bcc->muons->size(); ++i){
      if(m_bcc->muons->at(i).pt()>ptmax){
	ptmax = m_bcc->muons->at(i).pt();
	m_primlep = &m_bcc->muons->at(i);
      }
    }
  }

  return m_primlep;
}

TTbarGen* EventCalc::GetTTbarGen(){

  if(!m_ttgen){
    m_ttgen = new TTbarGen(m_bcc);
  }

  return m_ttgen;
}

void EventCalc::ApplyTauEnergySmearing(double factor)
{
  using namespace std;
  for(unsigned int i=0; i<m_bcc->taus->size(); ++i){
    LorentzVector vec = m_bcc->taus->at(i).v4();    
    LorentzVector svec;
    svec.SetPxPyPzE(vec.px()*factor, vec.py()*factor, vec.pz()*factor, vec.e()*factor);
    m_bcc->taus->at(i).set_v4(svec);
  }

}

std::vector<LorentzVector> EventCalc::NeutrinoReconstruction(const LorentzVector lepton, const LorentzVector met){

  
  TVector3 lepton_pT = toVector(lepton);
  lepton_pT.SetZ(0);
  
  TVector3 neutrino_pT = toVector(met);
  neutrino_pT.SetZ(0);
  
  const float mass_w = 80.399;
  float mu = mass_w * mass_w / 2 + lepton_pT * neutrino_pT;
  
  float A = - (lepton_pT * lepton_pT);
  float B = mu * lepton.pz();
  float C = mu * mu - lepton.e() * lepton.e() * (neutrino_pT * neutrino_pT);
  
  float discriminant = B * B - A * C;
  
  std::vector<LorentzVector> solutions;
  
  if (0 >= discriminant)
    {
      // Take only real part of the solution
      //
      LorentzVectorXYZE solution (0,0,0,0);
      solution.SetPx(met.Px());
      solution.SetPy(met.Py());
      solution.SetPz(-B / A);
      solution.SetE(toVector(solution).Mag());
      
      solutions.push_back(toPtEtaPhi(solution));
      
      //_solutions = 0 > discriminant ? 0 : 1;
    }
  else
    {
      discriminant = sqrt(discriminant);
      
      LorentzVectorXYZE solution (0,0,0,0);
      solution.SetPx(met.Px());
      solution.SetPy(met.Py());
      solution.SetPz((-B - discriminant) / A);
      solution.SetE(toVector(solution).Mag());
      
      solutions.push_back(toPtEtaPhi(solution));
      
      LorentzVectorXYZE solution2 (0,0,0,0);
      solution2.SetPx(met.Px());
      solution2.SetPy(met.Py());
      solution2.SetPz((-B + discriminant) / A);
      solution2.SetE(toVector(solution2).Mag());
      
      solutions.push_back(toPtEtaPhi(solution2));
      
      //_solutions = 2;
    }
  
  return solutions;
}

void EventCalc::FillHighMassTTbarHypotheses(){

  if(b_Reconstruction) return;
  b_Reconstruction=true;

  //clear hypothesis list
  m_bcc->recoHyps->clear();

  //find primary charged lepton
  Particle* lepton = GetPrimaryLepton();

  //reconstruct neutrino
  std::vector<LorentzVector> neutrinos = NeutrinoReconstruction( lepton->v4(), m_bcc->met->v4());
  
  ReconstructionHypothesis hyp;

  hyp.set_lepton(*lepton);

  //loop over neutrino solutions and jet assignments to fill hyotheses
  for(unsigned int i=0; i< neutrinos.size();++i){

    hyp.set_neutrino_v4(neutrinos[i]);
    LorentzVector wlep_v4 = lepton->v4()+neutrinos[i];

    unsigned int n_jets = m_bcc->jets->size();
    if(n_jets>10) n_jets=10; //avoid crashes in events with many jets
    unsigned int max_j = myPow(3, n_jets);
    for (unsigned int j=0; j < max_j; j++) {
      LorentzVector tophad_v4(0,0,0,0);
      LorentzVector toplep_v4 = wlep_v4;
      int hadjets=0;
      int lepjets=0;
      int num = j;
      hyp.clear_jetindices();
      for (unsigned int k=0; k<n_jets; k++) {
	// num is the k-th digit of j if you
	// write j in a base-3 system. According
	// to the value of this digit (which takes
	// values from 0 to 2,
	// in all possible combinations with the other digits),
	// decide how to treat the jet.
	
	if(num%3==0) {
	  tophad_v4 = tophad_v4 + m_bcc->jets->at(k).v4();
	  hyp.add_tophad_jet_index(k);
	  hadjets++;
	}
	
	if(num%3==1) {
	  toplep_v4 = toplep_v4 + m_bcc->jets->at(k).v4();
	  hyp.add_toplep_jet_index(k);
	  lepjets++;
	}
	//if(num%3==2); //do not take this jet
	
	//shift the trigits of num to the right:
	num /= 3;
      }
      
      //search jet with highest pt assigned to leptonic top
      float maxpt=-999;
      int maxind=-1;
      for(unsigned int i=0; i<hyp.toplep_jets_indices().size(); ++i){
	float pt = m_bcc->jets->at(hyp.toplep_jets_indices().at(i)).pt();
	if(pt>maxpt){
	  maxpt=pt;
	  maxind=hyp.toplep_jets_indices().at(i);
	}
      }
      hyp.set_blep_index(maxind);

      
      //fill only hypotheses with at least one jet assigned to each top quark
      if(hadjets>0 && lepjets>0){
	hyp.set_tophad_v4(tophad_v4);
	hyp.set_toplep_v4(toplep_v4);
	m_bcc->recoHyps->push_back(hyp);
      }
    }
  }
  
}


std::vector< PFParticle > EventCalc::GetJetPFParticles(Jet* jet){

  std::vector<PFParticle>* pfps = GetPFParticles();
  std::vector<unsigned int> jet_parts_ind = jet->pfconstituents_indices();

  std::vector<PFParticle> jetparticles;
  for(unsigned int ic=0; ic<jet_parts_ind.size(); ++ic){
    PFParticle p = pfps->at(jet_parts_ind[ic]);
    jetparticles.push_back(p);
  }
  return jetparticles;
}

std::vector< PFParticle >* EventCalc::GetJetPFParticles(){
  if(b_jetparticles) return &m_jetparticles;

  std::vector<PFParticle>* pfps = GetPFParticles();
  m_jetparticles.clear();

  for(unsigned int i=0; i< pfps->size(); ++i){
    if(pfps->at(i).isJetParticle()) m_jetparticles.push_back(pfps->at(i));
  }
  b_jetparticles=true;
  return &m_jetparticles;
}
 
std::vector< PFParticle >* EventCalc::GetIsoPFParticles(){
  if(b_isoparticles) return &m_isoparticles;

  std::vector<PFParticle>* pfps = GetPFParticles();
  m_isoparticles.clear();

  for(unsigned int i=0; i< pfps->size(); ++i){
    if(pfps->at(i).isIsoParticle()) m_isoparticles.push_back(pfps->at(i));
  }
  b_isoparticles=true;
  return &m_isoparticles;
}
 
std::vector< PFParticle >* EventCalc::GetPUIsoPFParticles(){
  if(b_puisoparticles) return &m_puisoparticles;

  std::vector<PFParticle>* pfps = GetPFParticles();
  m_puisoparticles.clear();

  for(unsigned int i=0; i< pfps->size(); ++i){
    if(pfps->at(i).isPUIsoParticle()) m_puisoparticles.push_back(pfps->at(i));
  }
  b_puisoparticles=true;
  return &m_puisoparticles;
}

double EventCalc::IntegratedJetShape(Jet* jet, double rmax, double rmin, E_JetType type ){

  double R = 0;
  if (type == e_CA8) R=0.8;
  else if (type == e_CA15) R=1.5; 
  else if (type == e_AK5) R=0.5;
  else {
    m_logger << ERROR << "EventCalc::IntegratedJetShape: can not determine radius for this type of jet, return Psi=-1" << SLogger::endmsg;
    return -1;
  }

  std::vector<PFParticle> pfps = GetJetPFParticles(jet);

  double pt_tot=0;
  double pt_frac=0;

  for(unsigned int i=0; i< pfps.size(); i++){
    if( pfps[i].deltaR(*jet) >= rmin && pfps[i].deltaR(*jet) <= rmax) pt_frac += pfps[i].pt();
    if( pfps[i].deltaR(*jet) >= rmin && pfps[i].deltaR(*jet) <= R) pt_tot += pfps[i].pt(); 
  }
  
  double psi=0;
  if(pt_tot>0) psi=pt_frac/pt_tot;
  return psi;
}

double EventCalc::JetCharge(Jet* jet){

  double charge = 0;

  std::vector<PFParticle> pfps = GetJetPFParticles(jet);
  for(unsigned int i=0; i< pfps.size(); i++){
    charge += pfps[i].charge();
  }
  return charge;
}

double EventCalc::EnergyWeightedJetCharge(Jet* jet, double kappa){

  double Q = 0;
  std::vector<PFParticle> pfps = GetJetPFParticles(jet);
  for(unsigned int i=0; i< pfps.size(); i++){
    Q += pfps[i].charge() * pow(pfps[i].energy(),kappa);
  }
  
  LorentzVector jet_v4_raw = jet->v4()*jet->JEC_factor_raw();
  double E_jet = jet_v4_raw.E();
  Q*= 1/pow(E_jet, kappa);
  return Q;
}

double EventCalc::JetMoment(Jet* jet, int n){
  
  double moment = 0;
  std::vector<PFParticle> pfps = GetJetPFParticles(jet);
  for(unsigned int i=0; i< pfps.size(); i++){
    moment += pow(jet->deltaR(pfps[i]),n);
  }
  
  if(pfps.size()!=0) moment /= 1.*pfps.size();
  return moment;
}

void EventCalc::PrintEventContent(){

  m_logger << INFO << "----------------- event content -----------------" << SLogger::endmsg;
  m_logger << INFO << "run: " << m_bcc->run <<  "   lumi block:" << m_bcc->luminosityBlock << "   event: " << m_bcc->event << SLogger::endmsg;
  m_logger << INFO << "MET = " << m_bcc->met->pt() << "   METphi = " << m_bcc->met->phi() <<  "   HTlep = " << GetHTlep() << SLogger::endmsg;
  if(m_bcc->electrons){m_logger << INFO << "Electrons:" << SLogger::endmsg;
    for(unsigned int i=0; i<m_bcc->electrons->size(); ++i){
      m_logger << INFO << "     " << i+1 << " pt = " << m_bcc->electrons->at(i).pt() <<"   eta = " << m_bcc->electrons->at(i).eta() 
	       << "   supercluster eta = " <<m_bcc->electrons->at(i).supercluster_eta() << "   IP wrt bsp = " << fabs(m_bcc->electrons->at(i).gsfTrack_dxy_vertex(m_bcc->pvs->at(0).x(), m_bcc->pvs->at(0).y()))
	       << "   pass conv veto = " << m_bcc->electrons->at(i).passconversionveto()  << "   mvaTrigV0 = " <<  m_bcc->electrons->at(i).mvaTrigV0()
	       << "   dEtaIn = " << m_bcc->electrons->at(i).dEtaIn() << "   sigmaIEtaIEta = " << m_bcc->electrons->at(i).sigmaIEtaIEta()
	       << "   HoverE = " << m_bcc->electrons->at(i).HoverE() << "   EcalEnergy = " << m_bcc->electrons->at(i).EcalEnergy() 
	       << "   EoverPIn = " << m_bcc->electrons->at(i).EoverPIn() << "   trackMomentumAtVtx = " << m_bcc->electrons->at(i).EcalEnergy()/m_bcc->electrons->at(i).EoverPIn()
	       << SLogger::endmsg;
    }
  }
  if(m_bcc->muons){m_logger << INFO << "Muons:" << SLogger::endmsg;
    for(unsigned int i=0; i<m_bcc->muons->size(); ++i){
      m_logger << INFO << "     " << i+1 << " pt = " << m_bcc->muons->at(i).pt() <<"   eta = " << m_bcc->muons->at(i).eta() << SLogger::endmsg;
    }
  }
  if(m_bcc->taus){m_logger << INFO << "Taus:" << SLogger::endmsg;
    for(unsigned int i=0; i<m_bcc->taus->size(); ++i){
      m_logger << INFO << "     " << i+1 << " pt = " << m_bcc->taus->at(i).pt() <<"   eta = " << m_bcc->taus->at(i).eta() << SLogger::endmsg;
    }
  }
  if(m_bcc->jets){m_logger << INFO << "Jets:" << SLogger::endmsg;
    for(unsigned int i=0; i<m_bcc->jets->size(); ++i){
      m_logger << INFO << "     " << i+1 << " pt = " << m_bcc->jets->at(i).pt() <<"   eta = " << m_bcc->jets->at(i).eta() << "   genjet_pt = " << m_bcc->jets->at(i).genjet_pt() << "   genjet_eta = " <<   m_bcc->jets->at(i).genjet_eta() <<SLogger::endmsg;
    }
  }
  if(m_bcc->topjets){m_logger << INFO << "TopJets:" << SLogger::endmsg;
    for(unsigned int i=0; i<m_bcc->topjets->size(); ++i){
      m_logger << INFO << "     " << i+1 << " pt = " << m_bcc->topjets->at(i).pt() <<"   eta = " << m_bcc->topjets->at(i).eta() << SLogger::endmsg;
    }
  }
  if(m_bcc->photons){m_logger << INFO << "Photons:" << SLogger::endmsg;
    for(unsigned int i=0; i<m_bcc->photons->size(); ++i){
      m_logger << INFO << "     " << i+1 << " pt = " << m_bcc->photons->at(i).pt() <<"   eta = " << m_bcc->photons->at(i).eta() << SLogger::endmsg;
    }
  }
}


void EventCalc::ProduceWeight(double weight)
{
  m_TotalWeight = m_TotalWeight * weight;
}

double EventCalc::GetWeight()
{

return m_TotalWeight;
}

void EventCalc::PrintGenParticles(string name)
{
  
  m_logger << INFO << "Printing generator particles for event " << GetEventNum() << " (name = " << name << ")" << SLogger::endmsg;
  if (GetGenParticles()->size()>0){

    std::cout << setw(10) << "index" << '|';
    std::cout << setw(10) << "pdgId" << '|';
    std::cout << setw(10) << "status" << '|';
    std::cout << setw(20) << "mother1" << '|';
    std::cout << setw(20) << "mother2" << '|';
    std::cout << setw(20) << "daughter1" << '|';
    std::cout << setw(20) << "daughter2" << '|';
    std::cout << setw(10) << "Px" << '|';
    std::cout << setw(10) << "Py" << '|';
    std::cout << setw(10) << "Pz"<< '|';
    std::cout << setw(10) << "energy" << '|';
    std::cout << setw(10) << "Pt" << '|';
    std::cout << setw(10) << "M" << std::endl; 

    std::cout.fill('-');
    std::cout << setw(11) << "+";
    std::cout << setw(11) << "+";
    std::cout << setw(11) << "+";
    std::cout << setw(21) << "+";
    std::cout << setw(21) << "+";
    std::cout << setw(21) << "+";
    std::cout << setw(21) << "+";
    std::cout << setw(11) << "+";
    std::cout << setw(11) << "+";
    std::cout << setw(11) << "+";
    std::cout << setw(11) << "+";
    std::cout << setw(11) << "+";
    std::cout << setw(11) << "" << std::endl; 
    std::cout.fill(' ');

    for (unsigned int i=0; i<GetGenParticles()->size(); ++i){
      GenParticle gp = GetGenParticles()->at(i);
      gp.Print(GetGenParticles());
    }
    std::cout << std::endl;
  } else {
    m_logger << INFO << "No generator particles found." << SLogger::endmsg;
  }

}

