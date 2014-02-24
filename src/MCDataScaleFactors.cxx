#include "include/MCDataScaleFactors.h"
#include "include/Utils.h"
#include <TMath.h>

LeptonScaleFactors::LeptonScaleFactors(std::vector<std::string> correctionlist)
{
    m_syst_shift = e_Default;
    m_muon_unc = false;
    m_ele_unc = false;
    m_tau_unc = false;
    m_tauele_unc = false;
    m_tau_eff_unc = false;
    
    if(correctionlist.size()%2!=0) {
        std::cerr<< "not a valid list of correction factors given to LeptonScaleFactors" <<std::endl;
        std::cerr<< "usage: \"<name> <weight> <name> <weight> <name> <weight> ...\" " <<std::endl;
	exit(EXIT_FAILURE);
    }

    for(unsigned int i=0; i< correctionlist.size()/2; ++i) {
        std::pair<std::string, double> correction (correctionlist[2*i], atof(correctionlist[2*i+1].c_str()));
        //std::cout << "Apply correction " << correction.first << " with factor " << correction.second <<std::endl;
        m_correctionlist.push_back(correction);
	m_apply = true;
    }

    if (correctionlist.size()==0){
      m_apply = false;
    }

    m_current_run = 0;
}

bool LeptonScaleFactors::IsUpToDate()
{
  EventCalc* calc = EventCalc::Instance();

  // no run dependency forseen in MC
  if (!calc->IsRealData()){
    if (m_current_run==0){
      m_current_run = 1;
      return false;
    } else {
      return true;
    }
  }
  // check run number in data
  if (m_current_run == calc->GetRunNum()){ 
    return true;
  } else {
    m_current_run = calc->GetRunNum();
    return false;
  }
}

void LeptonScaleFactors::FillWeights()
{
    if (!m_apply) return;

    // clear old scale factors
    m_mu_id.clear();
    m_mu_trig.clear();
    m_mu_iso.clear();
    m_ele_trig.clear();
    m_weights.clear();
    
    // initialize pointer to SF 2D-histogram
    m_ele_mva = NULL;

    // initialise scale factors to 1.0
    m_ele_trig.push_back(1.0);	
    m_ele_trig.push_back(0.0);	
    m_ele_trig.push_back(1.0);	

    // initialise arrays for eta bins
    for (int i=0; i<3; ++i){
      m_mu_id.push_back(std::vector<TGraphAsymmErrors*>());
      m_mu_trig.push_back(std::vector<TGraphAsymmErrors*>());
      m_mu_iso.push_back(std::vector<TGraphAsymmErrors*>());
    }
    
    // parameters of 1-d correction function for Ele30 trigger
    double Ele30_Trig_RelIso_Par[2] = { 0.9791, -0.5403 };

    // open file with scale factors 
    TFile* file = new TFile("$SFRAME_DIR/SFrameTools/efficiencies/muon_effs_2012_53x.root", "READ");
    if (!file->IsOpen()){
      file = new TFile("muon_effs_2012_53x.root", "READ");
    }
    if (!file->IsOpen()){
      std::cerr << "Could not find file with muon scale factors. Filename = muon_effs_2012_53x.root" << std::endl;
      std::cerr << "I looked in directory $SFRAME_DIR/SFrameTools/efficiencies/ and ./" << std::endl;
      std::cerr << "Please make sure the file is available." << std::endl;
      exit(EXIT_FAILURE);
    }

    TFile* egm_file = new TFile("$SFRAME_DIR/SFrameTools/efficiencies/EgammaPOG_TrigMVA_SF_2013.root", "READ");
    if (!egm_file->IsOpen()){
      std::cerr << "Could not find file with electron scale factors ";
      std::cerr << "in directory $SFRAME_DIR/SFrameTools/efficiencies/." << std::endl;
      std::cerr << " Please make sure the file is available." << std::endl;
      exit(EXIT_FAILURE);
    }

    double sum_mu_weights = 0;    
    // parse arguments of configuration
    for(unsigned int i=0; i<m_correctionlist.size(); ++i) {

        double weight = m_correctionlist[i].second;
	bool isok = false;
	
	TString name = m_correctionlist[i].first;
	if (name.Contains("Muon")){
	  sum_mu_weights += weight;
	}
	m_weights.push_back(weight);

	// ----------------- muons ------------------
	// muons: loop over eta bins
	for (int etabin=0; etabin<3; ++etabin){
	  
	  TString eta_name = TString::Format("eta%d", etabin);
	  
 	  //non isolated muons
 	  if(m_correctionlist[i].first == "MuonRunABCD") {
 	    m_mu_trig[etabin].push_back((TGraphAsymmErrors*) file->Get("2012ABCD/SF_" + eta_name + "_TRIG_Mu40"));
 	    m_mu_id[etabin].push_back((TGraphAsymmErrors*) file->Get("2012ABCD/SF_" + eta_name + "_ID_tight"));
 	    isok = true;
 	  }
	  
	   //isolated muons
	  else if(m_correctionlist[i].first == "IsoMuonRunABCD") {
	    m_mu_trig[etabin].push_back((TGraphAsymmErrors*) file->Get("2012ABCD/SF_" + eta_name + "_TRIG_IsoMu24"));
	    m_mu_iso[etabin].push_back((TGraphAsymmErrors*) file->Get("2012ABCD/SF_" + eta_name + "_ISO_tight"));
	    m_mu_id[etabin].push_back((TGraphAsymmErrors*) file->Get("2012ABCD/SF_" + eta_name + "_ID_tight"));
	    isok = true;
 	  }
	}
	
	// ----------------- electrons ------------------
	
	//trigger efficiency for electron trigger HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v*
	if (m_correctionlist[i].first == "HLT_Ele30") {
	  m_ele_trig[0] = Ele30_Trig_RelIso_Par[0];
	  m_ele_trig[1] = Ele30_Trig_RelIso_Par[1];
	  m_ele_trig[2] = weight;
	  isok = true;
	}	  

	// Scale Factor for Electron-ID based on Triggering MVA (Egamma-POG)
	if (m_correctionlist[i].first == "EGMTrigMVA") {

          m_ele_mva = (TH2F*) egm_file->Get("electronsDATAMCratio_FO_ID");
	  isok = true;
	}

	if (!isok){
	  std::cerr<< "No information found for lepton correction named " << m_correctionlist[i].first <<std::endl;
        }

    }
					       
    if (sum_mu_weights > 0.){
      for (unsigned int ii=0; ii<m_weights.size(); ++ii){
	m_weights[ii] = m_weights[ii]/sum_mu_weights;
      }
    }
    return;
}

int LeptonScaleFactors::GetMuonEtaBin(double eta) 
{
  int etabin=0;
  if(fabs(eta)<0.9) etabin=0;
  else if(fabs(eta)>=0.9 && fabs(eta)<1.2) etabin=1;
  else if(fabs(eta)>=1.2) etabin=2;
  return etabin;

}

int LeptonScaleFactors::GetBin(double xval, TGraphAsymmErrors* graph) 
{
  int bin = -1;
  double currx = -100;
  int i = 0;
  while (currx<xval){
    if (i >= graph->GetN()){
      bin = i-1;
      break;
    }
    double x, y;
    graph->GetPoint(i, x, y);
    double xlow = x - graph->GetErrorXlow(i);
    double xup = x + graph->GetErrorXhigh(i);
    bin = i-1;
    currx = xlow;
    ++i;
  }
  return bin;
}

double LeptonScaleFactors::GetMuonIDWeight()
{
  if (!m_apply) return 1.;
  if (!IsUpToDate()){
    FillWeights();
  }
  static EventCalc* calc = EventCalc::Instance();
  if (calc->GetMuons()->size()==0){
    return 1.;
  }
  Muon mu = calc->GetMuons()->at(0);
  int etabin = GetMuonEtaBin(mu.eta());
  double weight = 0.;
  if (m_mu_id[etabin].size()==0){
    return 1.;
  }
  for (unsigned int i=0; i<m_mu_id[etabin].size(); ++i){
    int ptbin = GetBin(mu.pt(), m_mu_id[etabin][i]);
    if (ptbin<0){
      weight += m_weights[i];
      continue;
    }
    double w = m_mu_id[etabin][i]->GetY()[ptbin];
    
    if (m_muon_unc){
      double sys = 0.005 * w; // systematic error on muon ID: 0.5%
      if (m_syst_shift==e_Down){
	double stat = m_mu_id[etabin][i]->GetEYlow()[ptbin];
	double err = TMath::Sqrt(stat*stat + sys*sys);
	w -= err;
      }
      if (m_syst_shift==e_Up){
	double stat = m_mu_id[etabin][i]->GetEYhigh()[ptbin];
	double err = TMath::Sqrt(stat*stat + sys*sys);
	w += err;
      }
    }
    weight += m_weights[i] * w;
  }

  return weight;
}

double LeptonScaleFactors::GetMuonTrigWeight()
{
  if (!m_apply) return 1.;
  if (!IsUpToDate()){
    FillWeights();
  }
  static EventCalc* calc = EventCalc::Instance();
  if (calc->GetMuons()->size()==0){
    return 1.;
  }
  Muon mu = calc->GetMuons()->at(0);
  int etabin = GetMuonEtaBin(mu.eta());
  double weight = 0.;
  if (m_mu_trig[etabin].size()==0){
    return 1.;
  }
  for (unsigned int i=0; i<m_mu_trig[etabin].size(); ++i){
    int ptbin = GetBin(mu.pt(), m_mu_trig[etabin][i]);
    if (ptbin<0){
      weight += m_weights[i];
      continue;
    }
    double w = m_mu_trig[etabin][i]->GetY()[ptbin];

    if (m_muon_unc){
      double sys = 0.002 * w; // systematic error on single muon trigger: 0.2%
      if (m_syst_shift==e_Down){
	double stat = m_mu_trig[etabin][i]->GetEYlow()[ptbin];
	double err = TMath::Sqrt(stat*stat + sys*sys);
	w -= err;
      }
      if (m_syst_shift==e_Up){
	double stat = m_mu_trig[etabin][i]->GetEYhigh()[ptbin];
	double err = TMath::Sqrt(stat*stat + sys*sys);
	w += err;
      }
    }
    weight += m_weights[i] * w;
  }
  return weight;
}

double LeptonScaleFactors::GetMuonIsoWeight()
{
  if (!m_apply) return 1.;
  if (!IsUpToDate()){
    FillWeights();
  }
  static EventCalc* calc = EventCalc::Instance();
  if (calc->GetMuons()->size()==0){
    return 1.;
  }
  Muon mu = calc->GetMuons()->at(0);
  int etabin = GetMuonEtaBin(mu.eta());
  double weight = 0.;
  if (m_mu_iso[etabin].size()==0){
    return 1.;
  }
  for (unsigned int i=0; i<m_mu_iso[etabin].size(); ++i){
    int ptbin = GetBin(mu.pt(), m_mu_iso[etabin][i]);
    if (ptbin<0){
      weight += m_weights[i];
      continue;
    }
    double w = m_mu_iso[etabin][i]->GetY()[ptbin];

    if (m_muon_unc){
      double sys = 0.002 * w; // systematic error on single muon isolation: 0.2%
      if (m_syst_shift==e_Down){
	double stat = m_mu_iso[etabin][i]->GetEYlow()[ptbin];
	double err = TMath::Sqrt(stat*stat + sys*sys);      
	w -= err;
      }
      if (m_syst_shift==e_Up){
	double stat = m_mu_iso[etabin][i]->GetEYhigh()[ptbin];
	double err = TMath::Sqrt(stat*stat + sys*sys);
	w += err;
      }
    }

    weight += m_weights[i] * w;
  }
  return weight;
}

double LeptonScaleFactors::GetMuonWeight()
{
  if (!m_apply) return 1.;
  if (!IsUpToDate()){
    FillWeights();
  }
  double id = GetMuonIDWeight();
  double trig = GetMuonTrigWeight();
  double iso = GetMuonIsoWeight();
  return id*trig*iso;
}

double LeptonScaleFactors::GetTauWeight()
{
   static EventCalc* calc = EventCalc::Instance();
   BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
   
   double weight = 1.;
   
   std::vector<Tau> fake_taus;
   std::vector<Tau> ele_fake_taus;
   bool fake = true;
   for(unsigned int j=0; j<bcc->taus->size(); ++j)
      {
         fake = true;
         Tau tau = bcc->taus->at(j);
         for(unsigned int i=0; i<bcc->genparticles->size(); ++i)
            {
               GenParticle genp = bcc->genparticles->at(i);
               double deltaR = genp.deltaR(tau);
               if (deltaR < 0.5 && abs(genp.pdgId())==15) fake =false; 
            }
         if(!fake) continue;
         bool fakeEle = false;
         for(unsigned int i=0; i<bcc->genparticles->size(); ++i)
            {
               GenParticle genp = bcc->genparticles->at(i);
               double deltaR = genp.deltaR(tau);
               if (deltaR < 0.5 && abs(genp.pdgId())==11 && genp.status()==3) fakeEle =true; 
            }
         if (fakeEle) {
            ele_fake_taus.push_back(tau);
            continue;
         }
         bool fakeMuon = false;
         for(unsigned int i=0; i<bcc->genparticles->size(); ++i)
            {
               GenParticle genp = bcc->genparticles->at(i);
               double deltaR = genp.deltaR(tau);
               if (deltaR < 0.5 && abs(genp.pdgId())==13 && genp.status()==3) fakeMuon =true; 
            }
         if (fakeMuon) continue;
         
         if (fake) fake_taus.push_back(tau);    
      }
   
   for(unsigned int i=0; i<ele_fake_taus.size(); ++i)
      {
         Tau tau = ele_fake_taus[i];
         if (!m_tauele_unc)
            {
               if (fabs(tau.eta()) <= 2.1) weight = weight*0.85;
               if (fabs(tau.eta()) > 2.1) weight = weight*0.65;
            } else 
            {
               if (m_syst_shift==e_Down)
                  {
                     if (fabs(tau.eta()) <= 2.1) weight = weight*0.85 - weight*0.85*0.2;
                     if (fabs(tau.eta()) > 2.1) weight = weight*0.65 - weight*0.65*0.25;	 
                  }
               if (m_syst_shift==e_Up)
                  {
                     if (fabs(tau.eta()) <= 2.1) weight = weight*0.85 + weight*0.85*0.2;
                     if (fabs(tau.eta()) > 2.1) weight = weight*0.65 + weight*0.65*0.25;	 
                  }
            }
      }
   for(unsigned int i=0; i<fake_taus.size(); ++i)
      {
         Tau tau = fake_taus[i];	      
         if (!m_tau_unc)
            {
               if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*1.0235;
               if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.7719;
               if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.4929;
               if (tau.pt() > 200) weight = weight*1.0813;
               
            } else 
            {
               
               if (m_syst_shift==e_Down)
                  {
                     if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*0.9567;
                     if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.6874;
                     if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.2870;
                     if (tau.pt() > 200) weight = weight*0.3583;
                  }
               if (m_syst_shift==e_Up)
                  {
                     if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*1.0898;
                     if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.8566;
                     if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.7064;
                     if (tau.pt() > 200) weight = weight*1.8769;
                  }
            }
      }
   return weight;
}
   

double LeptonScaleFactors::GetDecayModeFindingWeight()
{
  static EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  
  double weight = 1.;

  std::vector<Tau> fake_taus;
  bool fake = true;
   for(unsigned int j=0; j<bcc->taus->size(); ++j)
	{
	  fake = true;
	  Tau tau = bcc->taus->at(j);
	  for(unsigned int i=0; i<bcc->genparticles->size(); ++i)
	    {
	      GenParticle genp = bcc->genparticles->at(i);
	      double deltaR = genp.deltaR(tau);
	      if (deltaR < 0.5 && abs(genp.pdgId())==15) fake =false; 
	    }
	  if(!fake) continue;
	  bool fakeEle = false;
	  for(unsigned int i=0; i<bcc->genparticles->size(); ++i)
	    {
	      GenParticle genp = bcc->genparticles->at(i);
	      double deltaR = genp.deltaR(tau);
	      if (deltaR < 0.5 && abs(genp.pdgId())==11 && genp.status()==3) fakeEle = true; 
	    }
	  if (fakeEle)continue;
    
	  bool fakeMuon = false;
	  for(unsigned int i=0; i<bcc->genparticles->size(); ++i)
	    {
	      GenParticle genp = bcc->genparticles->at(i);
	      double deltaR = genp.deltaR(tau);
	      if (deltaR < 0.5 && abs(genp.pdgId())==13 && genp.status()==3) fakeMuon = true; 
	    }
	  if (fakeMuon) continue;

	  if (fake) fake_taus.push_back(tau);    
	}
   
   for(unsigned int i=0; i<fake_taus.size(); ++i)
     {
       Tau tau = fake_taus[i];	      
       weight = weight* 0.9354;
     }
   return weight;
}







double LeptonScaleFactors::GetTauEffUnc()
{
  static EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  
  
  double weight = 1;

  std::vector<Tau> real_taus;

   for(unsigned int j=0; j<bcc->taus->size(); ++j)
	{
	  bool real = false;
	  Tau tau = bcc->taus->at(j);
	  for(unsigned int i=0; i<bcc->genparticles->size(); ++i)
	    {
	      GenParticle genp = bcc->genparticles->at(i);
	      double deltaR = genp.deltaR(tau);
	      if (deltaR < 0.5 && abs(genp.pdgId())==15) real = true;
	    }
	  if (real) real_taus.push_back(tau);    
	}
    for(unsigned int i=0; i<real_taus.size(); ++i)
	{
	  Tau tau = real_taus[i];
	  if (!m_tau_eff_unc)
	    {
	      weight = 1.;
	    }
	  else
	    {
	      if (m_syst_shift==e_Down)
		{
		  weight = weight*0.94;	      
		}
	      if (m_syst_shift==e_Up)
		{
		  weight = weight*1.06;
		}
	    }
	}
    return weight;
}




double LeptonScaleFactors::GetElectronWeight()
{
  if (!m_apply) return 1.;
  if (!IsUpToDate()){
    FillWeights();
  }
  //double trig = GetElectronTrigWeight();
  double mva = GetElectronMVAIDWeight();
  return mva;
}

double LeptonScaleFactors::GetElectronTrigWeight()
{
  if (!m_apply) return 1.;
  if (!IsUpToDate()){
    FillWeights();
  }
  static EventCalc* calc = EventCalc::Instance();
  if (calc->GetElectrons()->size()==0){
    return 1.0;
  }
  Electron ele = calc->GetElectrons()->at(0);
  double iso = ele.relIso();
  double w = m_ele_trig[0] + m_ele_trig[1]*iso;

  // uncertainty	
  if (m_ele_unc){	
    double unc[] = {0.00829184, 0.187612};
    if (m_syst_shift==e_Down){
      w = (m_ele_trig[0]-unc[0]) + (m_ele_trig[1]-unc[1])*iso;
    } 
    if (m_syst_shift==e_Up){
      w = (m_ele_trig[0]+unc[0]) + (m_ele_trig[1]+unc[1])*iso;
    }
  }
  
  w *= m_ele_trig[2];
  if (w>1. || w<0.){ // sanity check
    w = 1.;
  }
  return w;
}

double LeptonScaleFactors::GetElectronMVAIDWeight()
{
  /* Data/MC scale factor for Electron-ID based on Triggering-MVA
   * provided by Egamma-POG for 22Jan2013-ReReco
   */

  if(!m_apply) return 1.;
  if(!IsUpToDate()) FillWeights();

  if(!m_ele_mva) return 1.;

  static EventCalc* calc = EventCalc::Instance();

  if(calc->GetElectrons()->size()==0) return 1.;

  double aScEta = std::abs(calc->GetElectrons()->at(0).supercluster_eta());
  double pt = calc->GetElectrons()->at(0).pt();
  if(pt>=m_ele_mva->GetYaxis()->GetXmax()) pt = m_ele_mva->GetYaxis()->GetXmax()-0.01;

  double w = m_ele_mva->GetBinContent(m_ele_mva->FindBin(aScEta,pt));

  // uncertainty	
  if (m_ele_unc){

    double unc = m_ele_mva->GetBinError(m_ele_mva->FindBin(aScEta,pt));
    if(m_syst_shift==e_Down) w -= unc;
    if(m_syst_shift==e_Up) w += unc;
  }

  if(w<=0.) w = 1.;// sanity check

  return w;
}

double LeptonScaleFactors::GetElectronORJetTrigWeight(const std::string& sys)
{
  /*  Data/MC scale factor for the "Ele30 OR PFJet320" trigger,
   *  measured as a function of the pT of the highest-pT jet;
   *  currently assigning a 1% flat systematic on the SF;
   *  Ref: https://indico.cern.ch/getFile.py/access?contribId=3&resId=0&materialId=slides&confId=290680
   */

  if(!m_apply) return 1.;

  EventCalc* calc = EventCalc::Instance();
  std::vector< Jet >* jets = calc->GetJets();
  if(!jets) return 1.;

  TFile* file = TFile::Open("$SFRAME_DIR/SFrameTools/efficiencies/Ele30_OR_PFJet320_trigSF.root");
  if (!file->IsOpen()){
    std::cerr << "Could not find file with Ele30_OR_PFJet320 trigger SF";
    std::cerr << " in $SFRAME_DIR/SFrameTools/efficiencies/.\n";
    std::cerr << "Please make sure the file is available.\n";
    exit(EXIT_FAILURE);
  }

  TF1* SFfit = (TF1*) file->Get("fit");
  if (!SFfit){
    std::cerr << "TF1 obj not found in Ele30_OR_PFJet320_trigSF.root.\n";
    std::cerr << "Please make sure the 'fit' function is available.\n";
    exit(EXIT_FAILURE);
  }

  float arg = jets->at(0).pt();
  if(arg < SFfit->GetXmin()){ arg = SFfit->GetXmin(); }
  else if(arg > SFfit->GetXmax()){ arg = SFfit->GetXmax(); }

  double w(1.);
  if(sys == "none"){ w = SFfit->Eval(arg); }
  else if(sys == "UP")  { w = SFfit->Eval(arg) * (1 + 0.01); }
  else if(sys == "DOWN"){ w = SFfit->Eval(arg) * (1 - 0.01); }
  else{
    std::cerr << "Incorrect argument for LeptonScaleFactors::GetElectronORJetTrigWeight().\n";
    std::cerr << "Must be either 'none', 'UP' or 'DOWN'. Exiting.\n";
    exit(EXIT_FAILURE);
  }

  file->Close();

  return w;
}

double LeptonScaleFactors::GetWeight()
{
    if (!m_apply) return 1.;
    double mu_weight = GetMuonWeight();
    double ele_weight = GetElectronWeight();
    return mu_weight * ele_weight;
}


BTaggingScaleFactors::BTaggingScaleFactors(
    E_BtagType btagtype, E_LeptonSelection lepton_selection, E_SystShift sys_bjets, E_SystShift sys_ljets
)
{
    m_sys_bjets = sys_bjets;
    m_sys_ljets = sys_ljets;
    
    m_btagtype = btagtype;
    m_lepton_selection = lepton_selection;

    _scale_btag=new BtagScale(btagtype);
    _eff_btag=new BtagEfficiency(btagtype, lepton_selection);

    _scale_ctag=new CtagScale(btagtype);
    _eff_ctag=new CtagEfficiency(btagtype, lepton_selection);

    _scale_light=new LtagScale(btagtype);
    _eff_light=new LtagEfficiency(btagtype, lepton_selection);
}


double BTaggingScaleFactors::GetWeight()
{
    EventCalc* calc = EventCalc::Instance();

    std::vector< Jet > *jets =  calc->GetJets();
    if(!jets) return 1.0;

    double scale_factor = 1.;

    for(unsigned int i=0; i<jets->size(); ++i) {

        Jet jet = jets->at(i);

        bool result = IsTagged(jet, m_btagtype);
        double scale_jet = 1.0;
        float jet_pt = jet.pt();
        float jet_eta = fabs(jet.eta());
	//do not consider jets outside b-tagging range
	if(jet_eta>=2.4){
	  continue;
	}

        switch(abs(jet.flavor())) {
        case 5: // b-quark
	    scale_jet = scale(result, jet_pt, jet_eta,
                              _scale_btag, _eff_btag,
                              m_sys_bjets);
            /*std::cout << "b jet pt: " << jet_pt << " is tagged: " << result << " scale: ";
            if (m_sys_bjets == e_Default)
                std::cout << _scale_btag->value(jet_pt) << " eff: " << _eff_btag->value(jet_pt);
            else if (m_sys_bjets == e_Up)
                std::cout << _scale_btag->value_plus(jet_pt) << " eff: " << _eff_btag->value_plus(jet_pt);
            else
                std::cout << _scale_btag->value_minus(jet_pt) << " eff: " << _eff_btag->value_minus(jet_pt);
            std::cout << " weight: " << scale_jet << std::endl;*/
            break;

        case 4: // c-quark
            scale_jet = scale(result, jet_pt, jet_eta,
                              _scale_ctag, _eff_ctag,
                              m_sys_bjets);
            /*std::cout << "b jet pt: " << jet_pt << " is tagged: " << result << " scale: ";
            if (m_sys_bjets == e_Default)
                std::cout << _scale_btag->value(jet_pt) << " eff: " << _eff_btag->value(jet_pt);
            else if (m_sys_bjets == e_Up)
                std::cout << _scale_btag->value_plus(jet_pt) << " eff: " << _eff_btag->value_plus(jet_pt);
            else
                std::cout << _scale_btag->value_minus(jet_pt) << " eff: " << _eff_btag->value_minus(jet_pt);
            std::cout << " weight: " << scale_jet << std::endl;*/
            break;

        case 3: // s-quark
        case 2: // d-quark
        case 1: // u-quark
        case 21: // gluon
            scale_jet = scale(result, jet_pt, jet_eta,
                              _scale_light, _eff_light,
                              m_sys_ljets);
            /*std::cout << "l jet pt: " << jet_pt << " is tagged: " << result << " scale: ";
            if (m_sys_ljets == e_Default)
                std::cout << _scale_light->value(jet_pt) << " eff: " << _eff_light->value(jet_pt);
            else if (m_sys_ljets == e_Up)  
                std::cout << _scale_light->value_plus(jet_pt) << " eff: " << _eff_light->value_plus(jet_pt);
            else
                std::cout << _scale_light->value_minus(jet_pt) << " eff: " << _eff_light->value_minus(jet_pt);
            std::cout << " weight: " << scale_jet << std::endl;*/
            break;

        default:
            break;
        }
        scale_factor *= scale_jet;
    }

    return scale_factor;
}


// Private
//
float BTaggingScaleFactors::scale(const bool &is_tagged,
                                  const float &jet_pt,
                                  const float &jet_eta,
                                  const BtagFunction* sf,
                                  const BtagFunction* eff,
                                  const E_SystShift &systematic)
{
    switch(systematic) {
    case e_Default:
        return is_tagged ?
	  sf->value(jet_pt,jet_eta) :
	  (1 - sf->value(jet_pt,jet_eta) * eff->value(jet_pt,jet_eta)) /
	  (1 - eff->value(jet_pt,jet_eta));
        break;

    case e_Up:
        return is_tagged ?
               sf->value_plus(jet_pt,jet_eta) :
               (1 - sf->value_plus(jet_pt,jet_eta) * eff->value_plus(jet_pt,jet_eta)) /
               (1 - eff->value_plus(jet_pt,jet_eta));
        break;

    case e_Down:
        return is_tagged ?
               sf->value_minus(jet_pt,jet_eta) :
               (1 - sf->value_minus(jet_pt,jet_eta) * eff->value_minus(jet_pt,jet_eta)) /
               (1 - eff->value_minus(jet_pt,jet_eta));
        break;

    default:
        std::cerr <<  "unsupported systematic" <<std::endl;
        break;
    }
    return 1.;
}


// Btag Scale
//
BtagScale::BtagScale(E_BtagType btagtype) : BtagFunction(btagtype) 
{
    // EPS13 prescription
    const float bins[] = {
        20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800
    };

    const float CSVTErrors[] = {
      0.0515703,
      0.0264008,
      0.0272757,
      0.0275565,
      0.0248745,
      0.0218456,
      0.0253845,
      0.0239588,
      0.0271791,
      0.0273912,
      0.0379822,
      0.0411624,
      0.0786307,
      0.0866832,
      0.0942053,
      0.102403 
    };
    const float CSVLErrors[] = {
      0.033299,
      0.0146768,
      0.013803,
      0.0170145,
      0.0166976,
      0.0137879,
      0.0149072,
      0.0153068,
      0.0133077,
      0.0123737,
      0.0157152,
      0.0175161,
      0.0209241,
      0.0278605,
      0.0346928,
      0.0350099 
    };
    switch(btagtype) {
    case e_CSVT: // EPS13 prescription
        _scale = new TF1("csvtb", "(0.927563+(1.55479e-05*x))+(-1.90666e-07*(x*x))", 20.0, 800.0);
        _bins.assign(bins, bins + 17);
        _errors.assign(CSVTErrors, CSVTErrors + 16);
        break;
    case e_CSVL:
        _scale = new TF1("csvlb", "0.997942*((1.+(0.00923753*x))/(1.+(0.0096119*x)))", 20.0, 800.0);
        _bins.assign(bins, bins + 17);
        _errors.assign(CSVLErrors, CSVLErrors + 16);
        break;
    default:
        std::cerr <<  "unsupported b-tagging operating point" <<std::endl;
        break;
    }
}


const unsigned int BtagScale::find_bin(const float &jet_pt, const float &jet_eta) const
{
    if (jet_pt < _bins.front())
        return 0;

    if (jet_pt > _bins.back())
        return _bins.size() - 2;

    unsigned int bin = 0;
    for(std::vector<float>::const_iterator bin_pt = _bins.begin();
            _bins.end() != ++bin_pt && *bin_pt < jet_pt;
            ++bin);

    return bin;
}


float BtagScale::value(const float &jet_pt, const float &jet_eta) const
{
    if (_scale->GetXmin() > jet_pt)
      return value(_scale->GetXmin(), jet_eta);

    if (_scale->GetXmax() < jet_pt)
        return value(_scale->GetXmax(), jet_eta);

    return _scale->Eval(jet_pt, jet_eta);
}


float BtagScale::error(const float &jet_pt, const float &jet_eta) const
{
    if (_scale->GetXmin() > jet_pt)
        return 2.0 * error(_scale->GetXmin(), jet_eta);

    if (_scale->GetXmax() < jet_pt)
        return 2.0 * error(_scale->GetXmax(), jet_eta);

    return _errors.at(find_bin(jet_pt, jet_eta));
}


// Ctag scale
//
float CtagScale::error(const float &jet_pt, const float &jet_eta) const
{
  return 2 * BtagScale::error(jet_pt, jet_eta);
}

// Light-tag scale
//
LtagScale::LtagScale(E_BtagType btagtype) : BtagFunction(btagtype)
{
    switch(btagtype) {
    case e_CSVT: // EPS13 prescription
         _scale = new TF2("csvtl","((1.00462+(0.00325971*x))+(-7.79184e-06*(x*x)))+(5.22506e-09*(x*(x*x)))",20.,800.0,0,2.4);
        _scale_plus = new TF2("cvstl_plus","((1.16361+(0.00464695*x))+(-1.09467e-05*(x*x)))+(7.21896e-09*(x*(x*x)))",20.,800.0,0,2.4);
        _scale_minus = new TF2("cvstl_minus","((0.845757+(0.00186422*x))+(-4.6133e-06*(x*x)))+(3.21723e-09*(x*(x*x)))",20.,800.0,0,2.4);
        break;
    case e_CSVL:
        _scale = new TF2("csvll","(((1.01177+(0.0023066*x))+(-4.56052e-06*(x*x)))+(2.57917e-09*(x*(x*x))))*(0.5>y)+(((0.975966+(0.00196354*x))+(-3.83768e-06*(x*x)))+(2.17466e-09*(x*(x*x))))*(y>=0.5)*(1.0>y)+(((0.93821+(0.00180935*x))+(-3.86937e-06*(x*x)))+(2.43222e-09*(x*(x*x))))*(y>=1.0)*(1.5>y)+(((1.00022+(0.0010998*x))+(-3.10672e-06*(x*x)))+(2.35006e-09*(x*(x*x))))*(y>=1.5)" ,20.,800.0,0,2.4);
        _scale_plus = new TF2("cvsll_plus","(((1.04582+(0.00290226*x))+(-5.89124e-06*(x*x)))+(3.37128e-09*(x*(x*x))))*(0.5>y)+(((1.00683+(0.00246404*x))+(-4.96729e-06*(x*x)))+(2.85697e-09*(x*(x*x))))*(y>=0.5)*(1.0>y)+(((0.964787+(0.00219574*x))+(-4.85552e-06*(x*x)))+(3.09457e-09*(x*(x*x))))*(y>=1.0)*(1.5>y)+(((1.03039+(0.0013358*x))+(-3.89284e-06*(x*x)))+(3.01155e-09*(x*(x*x))))*(y>=1.5)",20.,800.0,0,2.4);
        _scale_minus = new TF2("cvsll_minus","(((0.977761+(0.00170704*x))+(-3.2197e-06*(x*x)))+(1.78139e-09*(x*(x*x))))*(0.5>y)+(((0.945135+(0.00146006*x))+(-2.70048e-06*(x*x)))+(1.4883e-09*(x*(x*x))))*(y>=0.5)*(1.0>y)+(((0.911657+(0.00142008*x))+(-2.87569e-06*(x*x)))+(1.76619e-09*(x*(x*x))))*(y>=1.0)*(1.5>y)+(((0.970045+(0.000862284*x))+(-2.31714e-06*(x*x)))+(1.68866e-09*(x*(x*x))))*(y>=1.5)",20.,800.0,0,2.4);
        break;
    default:
        std::cerr <<  "unsupported b-tagging operating point" <<std::endl;
        break;
    }
}




float LtagScale::value(const float &jet_pt, const float &jet_eta) const
{
    if (_scale->GetXmin() > jet_pt)
        return value(_scale->GetXmin(), jet_eta);

    if (_scale->GetXmax() < jet_pt)
        return value(_scale->GetXmax(), jet_eta);

    return _scale->Eval(jet_pt, jet_eta);
}


float LtagScale::value_plus(const float &jet_pt, const float &jet_eta) const
{
    if (_scale_plus->GetXmin() > jet_pt)
    {
      double error = 2.0*(value_plus(_scale_plus->GetXmin(), jet_eta) - value(_scale_plus->GetXmin(), jet_eta));
        return value(_scale_plus->GetXmin(), jet_eta) + error;
    }

    if (_scale_plus->GetXmax() < jet_pt)
    {
      double error = 2.0*(value_plus(_scale_plus->GetXmax(), jet_eta) - value(_scale_plus->GetXmax(), jet_eta));
        return value(_scale_plus->GetXmax(), jet_eta) + error;
    }

    return _scale_plus->Eval(jet_pt, jet_eta);
}


float LtagScale::value_minus(const float &jet_pt, const float &jet_eta) const
{
    if (_scale_minus->GetXmin() > jet_pt)
    {
        double error = 2.0*(value(_scale->GetXmin(), jet_eta) - value_minus(_scale_minus->GetXmin(), jet_eta));
        double scale = value(_scale_minus->GetXmin(), jet_eta) - error;
        if ( scale >= 0.0 ) return scale;
        return 0.0;
    }

    if (_scale_minus->GetXmax() < jet_pt)
    {
        double error = 2.0*(value(_scale->GetXmax(), jet_eta) - value_minus(_scale_minus->GetXmax(), jet_eta));
        double scale = value(_scale_minus->GetXmax(), jet_eta) - error;
        if ( scale >= 0.0 ) return scale;
        return 0.0;
    }

    return _scale_minus->Eval(jet_pt, jet_eta);
}


// Btag Efficiency
//
BtagEfficiency::BtagEfficiency(E_BtagType btagtype, E_LeptonSelection leptonsel) : BtagFunction(btagtype)
{
    const float bins[] = {
        20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 120.0, 160.0, 210.0, 260.0, 320.0, 400.0, 500.0, 600.0, 800.0, 1600.0
    };



    const float CSVTEfficiencies[] = {
        0.45995393113468314, 0.46970557780785482, 0.47945722448102651, 0.4892088711541982, 
        0.49896051782736989, 0.52736480487115456, 0.51718656298800869, 0.51778381299772003, 
        0.49273414735703341, 0.40596266795544039, 0.34807215987045181, 0.31600602509673009,
        0.27177222600495071, 0.17082550051964149, 0.10729452077597022, 0.088296058189995086, 0.11093069992211117
    };

    const float CSVTEfficiencies_mu[] = {
      0, 0, 0, 0.527496 , 0.538724 , 0.555376 , 0.540102 , 0.528351 , 0.50122 , 0.41426 , 0.371213 , 0.311879 , 0.264041 , 0.215617 , 0.208605 , 0.150797 , 0.130312
    };
    const float CSVLEfficiencies_mu[] = {
      0, 0, 0, 0.843019 , 0.836389 , 0.847278 , 0.839844 , 0.833074 , 0.835171 , 0.820583 , 0.804997 , 0.762067 , 0.730621 , 0.707983 , 0.655583 , 0.607146 , 0.495463
    };

    if (btagtype == e_CSVT && leptonsel == e_Electron) { 
        _bins.assign(bins, bins + 18);
        _values.assign(CSVTEfficiencies, CSVTEfficiencies + 17);
    } 
    else if (btagtype == e_CSVT && leptonsel == e_Muon) { 
      _bins.assign(bins, bins + 18);
      _values.assign(CSVTEfficiencies_mu, CSVTEfficiencies_mu + 17);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Muon) { 
      _bins.assign(bins, bins + 18);
      _values.assign(CSVLEfficiencies_mu, CSVLEfficiencies_mu + 17);
    }
    else {
        std::cerr <<  "unsupported b-tagging operating point and lepton selection" <<std::endl;
    }
}


const unsigned int BtagEfficiency::find_bin(const float &jet_pt, const float &jet_eta ) const
{
    if (jet_pt < _bins.front())
        return 0;

    if (jet_pt > _bins.back())
        return _bins.size() - 2;

    unsigned int bin = 0;
    for(std::vector<float>::const_iterator bin_pt = _bins.begin();
            _bins.end() != ++bin_pt && *bin_pt < jet_pt;
            ++bin);

    return bin;
}


float BtagEfficiency::value(const float &jet_pt, const float &jet_eta) const
{
  return _values.at(find_bin(jet_pt, jet_eta));
}


CtagEfficiency::CtagEfficiency(E_BtagType btagtype, E_LeptonSelection leptonsel) : BtagEfficiency(btagtype, leptonsel)
{
    const float bins[] = {
        20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 120.0, 160.0, 210.0, 260.0, 320.0, 400.0, 1600.0
    };

    const float CSVTEfficiencies[] = {
        0.078277168566833921, 0.074915742421118273, 0.071554316275402624, 0.068192890129686975,
        0.064831463983971327, 0.056055860000873536, 0.053423145350989264, 0.049437287248529672, 
        0.051175126071014654, 0.031147321304156712, 0.028672334664972543, 0.017483901927108567, 
        0.012445728161565342, 0.013059366026755743
    };

    const float CSVTEfficiencies_mu[] = {
      0, 0, 0 , 0.0746531 , 0.0648892 , 0.0661782 , 0.0593803 , 0.0483994 , 0.0424287 , 0.0284464 , 0.0307131 , 0.0275881 , 0.0173931 , 0.0147546
    };
    const float CSVLEfficiencies_mu[] = {
      0, 0, 0 , 0.438554 , 0.409188 , 0.417346 , 0.417047 , 0.385195 , 0.376154 , 0.372287 , 0.342093 , 0.324201 , 0.312492 , 0.252123
    };

    if (btagtype == e_CSVT && leptonsel == e_Electron) {
        _bins.assign(bins, bins + 15);
        _values.assign(CSVTEfficiencies, CSVTEfficiencies + 14);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Muon) { 
      _bins.assign(bins, bins + 15);
      _values.assign(CSVTEfficiencies_mu, CSVTEfficiencies_mu + 14);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Muon) { 
      _bins.assign(bins, bins + 15);
      _values.assign(CSVLEfficiencies_mu, CSVLEfficiencies_mu + 14);
    }
    else {
        std::cerr <<  "unsupported b-tagging operating point and lepton selection" <<std::endl;
    }
}


// Light Efficiency
//
LtagEfficiency::LtagEfficiency(E_BtagType btagtype, E_LeptonSelection leptonsel) : BtagEfficiency(btagtype, leptonsel)
{
    // BtagFunction::BtagFunction(btagtype);

    const float bins[] = {
        20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 120.0, 160.0, 210.0, 260.0, 320.0, 400.0, 500.0, 600.0, 1600.0
    };

    const float CSVTEfficiencies[] = {
        0.0041961524551977604, 0.0044655765042000087, 0.0047350005532022571, 0.0050044246022045054,
        0.0052738486512067537, 0.0045889898119126932, 0.0050790777581588833, 0.0059039415342432922,
        0.0060358742501952448, 0.0058222771697384801, 0.0060753523061352734, 0.0042071929862444917, 
        0.0061857466264373896, 0.0046817449333456688, 0.016176856044657056, 0.0065362936670525645
    };

    const float CSVTEfficiencies_mu[] = {
      0, 0, 0, 0.00569484 , 0.00548976 , 0.00547509 , 0.00712279 , 0.00675918 , 0.00670482 , 0.00444024 , 0.0045339 , 0.00446893 , 0.00537538 , 0.00461259 , 0.00792868 , 0.00387412
    };

    const float CSVLEfficiencies_mu[] = {
      0, 0, 0, 0.13864 , 0.112884 , 0.110701 , 0.109742 , 0.0995396 , 0.099064 , 0.104746 , 0.101926 , 0.100509 , 0.104804 , 0.108763 , 0.104505 , 0.101038
    };

    if (btagtype == e_CSVT && leptonsel == e_Electron) {
        _bins.assign(bins, bins + 17);
        _values.assign(CSVTEfficiencies, CSVTEfficiencies + 16);
    } 
    else if (btagtype == e_CSVT && leptonsel == e_Muon) { 
      _bins.assign(bins, bins + 17);
      _values.assign(CSVTEfficiencies_mu, CSVTEfficiencies_mu + 16);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Muon) { 
      _bins.assign(bins, bins + 17);
      _values.assign(CSVLEfficiencies_mu, CSVLEfficiencies_mu + 16);
    }
    else {
        std::cerr <<  "unsupported b-tagging operating point and lepton selection" <<std::endl;
    }
}



JetpTReweightingInWJets::JetpTReweightingInWJets()
{
    m_syst_shift = e_Default;
    m_jetpTreweigting_unc = false;
}


double JetpTReweightingInWJets::GetWeight()
{
   EventCalc* calc = EventCalc::Instance();
   BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
   double weight = 1.;
   
   if (!m_jetpTreweigting_unc)
      {
         if(bcc->genjets->size() > 0) {
            float genpt = bcc->genjets->at(0).pt();
            weight = TMath::Exp(-0.690-0.009 *genpt) +  0.879;
         }
      } else
      {
         if (m_syst_shift==e_Down)
            {
               weight = 1.;
            }
         if (m_syst_shift==e_Up)
            {
               if(bcc->genjets->size() > 0) {
                  float genpt = bcc->genjets->at(0).pt();
                  weight = TMath::Exp(-0.690-0.009 *genpt) +  0.879;
                  weight = weight *weight;
               }
            }
      }
   
   return weight;
}



