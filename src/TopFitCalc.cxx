// Dear emacs, this is -*- c++ -*-

#include "include/TopFitCalc.h"
#include "include/SubJetTagger.h"
#include "include/Utils.h"
#include "TSystem.h"
#include <stdio.h>
#include <iostream>
#include <stdexcept>

TopFitCalc* TopFitCalc::m_instance = NULL;

TopFitCalc* TopFitCalc::Instance()
{
  // Get a pointer to the object handler.
  // This is the only way to access this class,
  // since it's a singleton. This method is accessible
  // from everywhere.

  if (m_instance == NULL){
    m_instance = new TopFitCalc();
  }
  return m_instance;
}



TopFitCalc::TopFitCalc() : m_logger("TopFitCalc")
{
  // constructor: initialise all variables

  m_logger << DEBUG << "Constructor of TopFitCalc called." << SLogger::endmsg;

  positiv = new TMinuit(5);

  b_Reconstruction = false;

  //gSystem->Load("libMinuit");

  //Reset();

}

void TopFitCalc::Reset()
{
  b_Reconstruction = false;
}


TopFitCalc::~TopFitCalc()
{
  // default destructor

}

void TopFitCalc::FillHighMassTTbarHypotheses(){

  if(b_Reconstruction) return;
  b_Reconstruction=true;

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* m_bcc = calc->GetBaseCycleContainer();

  //clear hypothesis list
  m_bcc->recoHyps->clear();

  //find primary charged lepton
  Particle* lepton = calc->GetPrimaryLepton();

  //reconstruct neutrino
  std::vector<LorentzVector> neutrinos = calc->NeutrinoReconstruction( lepton->v4(), m_bcc->met->v4());
  //std::vector<LorentzVector> neutrinos = NeutrinoFitPolar( lepton->v4(), m_bcc->met->v4() );

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


double DeltaPolarNeutrino(double PhiN, double metPx, double metPy, double PhiLep, double PtLep)
{
  using namespace std;

  double PyN;
  double PxN;
  const double mass_w = 80.399;

  double num = 10.e-7;

  if(1-cos(fabs(PhiLep-PhiN)> PI ? 2*PI-fabs(PhiLep-PhiN) : fabs(PhiLep-PhiN )) < num){
    PyN = 0.5*mass_w*mass_w* sin(PhiN)/(PtLep*num);
    PxN = 0.5*mass_w*mass_w* cos(PhiN)/(PtLep*num);
  }
  else{
    PyN = 0.5*mass_w*mass_w* sin(PhiN)/(PtLep*(1-cos(fabs(PhiLep-PhiN)> PI ? 2*PI-fabs(PhiLep-PhiN) : fabs(PhiLep-PhiN ))));
    PxN = 0.5*mass_w*mass_w* cos(PhiN)/(PtLep*(1-cos(fabs(PhiLep-PhiN)> PI ? 2*PI-fabs(PhiLep-PhiN) : fabs(PhiLep-PhiN ))));
  }

  return pow(PxN-metPx,2)+pow(PyN-metPy,2);

}

static void polarminuitfunc(int& nDim, double* gout, double& result, double par[], int flg){
  result = DeltaPolarNeutrino(par[0],par[1],par[2],par[3],par[4]);

}

std::vector<LorentzVector> TopFitCalc::NeutrinoFitPolar(const LorentzVector lepton, const LorentzVector met)
{

  TVector3 lepton_pT = toVector(lepton);
  lepton_pT.SetZ(0);

  TVector3 neutrino_pT = toVector(met);
  neutrino_pT.SetZ(0);

  const double mass_w = 80.399;


  double min = -2*PI;
  double max = 2*PI;
  double start = met.phi();
  double step = 10.e-5;

  double mu = mass_w * mass_w / 2 + lepton_pT * neutrino_pT;
  double A = - (lepton_pT * lepton_pT);
  double B = mu * lepton.pz();
  double C = mu * mu - lepton.e() * lepton.e() * (neutrino_pT * neutrino_pT);

  double discriminant = B * B - A * C;

  std::vector<LorentzVector> solutions;


  if (0 >= discriminant)
    {

    double resultPhi = 0;
    double error = 0;
    int ierflg;
    double* arg = new double[1];

    positiv->SetPrintLevel(-1); // -1 quiet, 0 normal, 1 verbose; Preset 0

    positiv->SetFCN(polarminuitfunc);

    positiv->DefineParameter(0,"PhiN",start, step,  min, max);
    positiv->DefineParameter(1,"metPx",met.px(),0,0,0);
    positiv->DefineParameter(2,"metPy",met.py(),0,0,0);
    positiv->DefineParameter(3,"PhiLep",lepton.phi(),0,0,0);
    positiv->DefineParameter(4,"PtLep",lepton.pt(),0,0,0);

    positiv->FixParameter(1);
    positiv->FixParameter(2);
    positiv->FixParameter(3);
    positiv->FixParameter(4);

    positiv->SetMaxIterations(500);


    arg[0]= 2;
    positiv->mnexcm("SET STR",arg,1,ierflg);

    positiv->Migrad();

    positiv->GetParameter(0,resultPhi,error);

    delete[] arg;

    if(resultPhi != resultPhi)
      {
	std::cerr << "neutrino phi is NAN " << std::endl;
      }

    if(resultPhi > PI) resultPhi = resultPhi-2*PI;
    if(resultPhi < PI) resultPhi = resultPhi+2*PI;


    double PyN;
    double PxN;

    double num = 10.e-7;

    if(1-cos(deltaPhiAbs(lepton.phi(), resultPhi)) < num){
      PyN = 0.5*mass_w*mass_w* sin(resultPhi)/(lepton.pt()*num);
      PxN = 0.5*mass_w*mass_w* cos(resultPhi)/(lepton.pt()*num);
    }
    else{
      PyN = 0.5*mass_w*mass_w* sin(resultPhi)/(lepton.pt()*(1-cos(deltaPhiAbs(lepton.phi(), resultPhi))));
      PxN = 0.5*mass_w*mass_w* cos(resultPhi)/(lepton.pt()*(1-cos(deltaPhiAbs(lepton.phi(), resultPhi))));
    }

    LorentzVectorXYZE neutrino_result(0,0,0,0);
    neutrino_result.SetPx(PxN);
    neutrino_result.SetPy(PyN);

    double pzfit =  lepton.pz()*neutrino_result.pt()/lepton.pt();

    LorentzVectorXYZE solution (0,0,0,0);
    solution.SetPx(PxN);
    solution.SetPy(PyN);
    solution.SetPz(pzfit);
    solution.SetE(toVector(solution).Mag());


    solutions.push_back(toPtEtaPhi(solution));

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

   }

  return solutions;
}


void TopFitCalc::CalculateSelection()
{

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* m_bcc = calc->GetBaseCycleContainer();
  double deltaR_Lep_Tophad = .8;
  double deltaR_Jet_Tophad = 1.3;
  double mjet=0;
  //clear hypothesis list
  //m_bcc->recoHyps->clear();

  //find primary charged lepton
  Particle* lepton = calc->GetPrimaryLepton();

  std::vector<Jet>* antikjets = m_bcc->jets;
  std::vector<TopJet>* cajets = m_bcc->topjets;

    std::vector<LorentzVector> neutrinos = calc->NeutrinoReconstruction(lepton->v4(),m_bcc->met->v4());
    // std::vector<LorentzVector> neutrinos = NeutrinoFitPolar(lepton->v4(),m_bcc->met->v4());

  ReconstructionHypothesis hyp;

  hyp.set_lepton(*lepton);
  double cajets_pt;
  if(cajets->size()>0) {
   cajets_pt   =  cajets->at(0).pt();
  int caposi_pt = 0;

  for(unsigned int m = 0; m<cajets->size(); ++m)
    {
      TopJet cajet = cajets->at(m);

      if(cajets_pt < cajet.pt())
	{
	  cajets_pt = cajet.pt();
	  caposi_pt = m;
	}

    }
   LorentzVector top_had;
  if(cajets->size()>0) top_had = cajets->at(caposi_pt).v4();
  mjet=top_had.M();
     for(unsigned int m = 0; m<cajets->size(); ++m)
    {
       TopJet cajet = cajets->at(m);
       if(deltaR_Lep_Tophad < deltaR(cajet.v4(),lepton->v4())){

	int n_jets = antikjets->size();
	if(n_jets>6) n_jets=6;
	int max_j = myPow(3, n_jets);


	for(unsigned int i = 0; i < neutrinos.size();++i){

	  Particle wboson_lep;
	  wboson_lep.set_v4(lepton->v4()+neutrinos.at(i));


	  for(int j=0; j<max_j; ++j){
	    LorentzVector top_lep(0,0,0,0);
	    LorentzVector b_lep(0,0,0,0);
	    int num = j;
	    hyp.clear_jetindices();
	    for(unsigned int p=0; p<antikjets->size(); ++p){
	      if(deltaR(top_had,antikjets->at(p).v4())> deltaR_Jet_Tophad && num%3==0){
		b_lep = b_lep + antikjets->at(p).v4();
		top_lep = wboson_lep.v4() + b_lep;
		hyp.set_blep_index(p);
		hyp.set_blep_v4(b_lep);
		hyp.add_toplep_jet_index(p);
		hyp.add_tophad_jet_index(caposi_pt);

		hyp.set_neutrino_v4(neutrinos[i]);

		double egroomed = sqrt(cajet.v4().P2()+mjet*mjet);
		top_had.SetPxPyPzE(cajet.v4().Px(),cajet.v4().Py(),cajet.v4().Pz(),egroomed);
		hyp.set_tophad_v4(top_had);

		hyp.set_toplep_v4(top_lep);

		m_bcc->recoHyps->push_back(hyp);
	      }
	      num/=3;
	    }

	  }}}}

  }}


void TopFitCalc::CalculateTopTag()
{

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* m_bcc = calc->GetBaseCycleContainer();
  m_bcc->recoHyps->clear();

  //find primary charged lepton
  Particle* lepton = calc->GetPrimaryLepton();

  std::vector<Jet>* antikjets = m_bcc->jets;
  std::vector<TopJet>* cajets = m_bcc->topjets;


  std::vector<LorentzVector> neutrinos = calc->NeutrinoReconstruction(lepton->v4(),m_bcc->met->v4());
  //std::vector<LorentzVector> neutrinos = NeutrinoFitPolar(lepton->v4(),m_bcc->met->v4());

  ReconstructionHypothesis hyp;

  hyp.set_lepton(*lepton);

  int best_neutrino = 0;


  int NCAJets = cajets->size();

  int caposi=-1;
  double deltaR_Lep_Tophad = .8;
  double deltaR_Jet_Tophad = 1.3;

  CMSTopTagger toptag;
  toptag.SetTau32Cut();

  for(unsigned int m = 0; m<cajets->size(); ++m)
    {
      TopJet cajet = cajets->at(m);

      if(toptag.Tag(cajet) && deltaR_Lep_Tophad < deltaR(cajet.v4(),lepton->v4())){
	caposi = m;
	//deltaR_Lep_Tophad = delR(cajet.v4(),lepton->v4());

	//if(caposi==-1) return;
	LorentzVector top_had = cajet.v4();
        LorentzVector allsubjets(0,0,0,0);

        for(int j=0; j<cajet.numberOfDaughters(); ++j) {
            allsubjets += cajet.subjets()[j].v4();
        }
        double mgroomed = allsubjets.M();

	int n_jets = antikjets->size();
	if(n_jets>6) n_jets=6;
	int max_j = myPow(3, n_jets);


	for(unsigned int i = 0; i < neutrinos.size();++i){

	  Particle wboson_lep;
	  wboson_lep.set_v4(lepton->v4()+neutrinos.at(i));


	  for(int j=0; j<max_j; ++j){
	    LorentzVector top_lep(0,0,0,0);
	    LorentzVector b_lep(0,0,0,0);
	    int num = j;
	    hyp.clear_jetindices();
	    for(unsigned int p=0; p<antikjets->size(); ++p){
	      if(deltaR(top_had,antikjets->at(p).v4())> deltaR_Jet_Tophad && num%3==0){
		b_lep = b_lep + antikjets->at(p).v4();
		top_lep = wboson_lep.v4() + b_lep;
		hyp.set_blep_index(p);
		hyp.set_blep_v4(b_lep);
		hyp.add_toplep_jet_index(p);
		hyp.add_tophad_jet_index(caposi);

		hyp.set_neutrino_v4(neutrinos[i]);

		double egroomed = sqrt(cajet.v4().P2()+mgroomed*mgroomed);
		top_had.SetPxPyPzE(cajet.v4().Px(),cajet.v4().Py(),cajet.v4().Pz(),egroomed);
		hyp.set_tophad_v4(top_had);

		hyp.set_toplep_v4(top_lep);

		m_bcc->recoHyps->push_back(hyp);
	      }
	      num/=3;
	    }
	  }
	}
      }
    }
}

