#include "include/MCDataScaleFactors.h"
#include "include/Utils.h"
#include <TMath.h>



LeptonScaleFactors::LeptonScaleFactors(std::vector<std::string> correctionlist, TString channel)
{
    m_syst_shift = e_Default;
    m_muon_unc = false;
    m_ele_unc = false;
    m_tau_unc = false;
    m_tauele_unc = false;
    m_tau_eff_unc = false;
    m_channel = channel;

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
         if (m_channel == "inclusive_LQ")
            {
               if (!m_tau_unc)
                  {
                     if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*1.04765;
                     if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.857812;
                     if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.59419;
                     if (tau.pt() > 200) weight = weight*1.06272;
                  } else
                  {

                     if (m_syst_shift==e_Down)
                        {
                           if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*0.9352;
                           if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.6769;
                           if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.287;
                           if (tau.pt() > 200) weight = weight*0.3583;
                        }
                     if (m_syst_shift==e_Up)
                        {
                           if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*1.1113;
                           if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.8671;
                           if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.7064;
                           if (tau.pt() > 200) weight = weight*1.8769;
                        }
                  }
            }
         if (m_channel == "combination_LQ")
            {
               if (!m_tau_unc)
                  {
                     if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*1.04765;
                     if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.857812;
                     if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.59419;
                     if (tau.pt() > 200) weight = weight*1.06272;
                  } else
                  {

                     if (m_syst_shift==e_Down)
                        {
                           if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*0.9142;
                           if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.6652;
                           if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.281;
                           if (tau.pt() > 200) weight = weight*0.3525;
                        }
                     if (m_syst_shift==e_Up)
                        {
                           if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*1.1324;
                           if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.8788;
                           if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.7104;
                           if (tau.pt() > 200) weight = weight*1.8824;
                        }
                  }
            }
         if (m_channel == "inclusive_ClosureTest")
            {
              if (!m_tau_unc)
                  {
                     if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*1.04765;
                     if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.857812;
                     if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.59419;
                     if (tau.pt() > 200) weight = weight*1.06272;
                  } else
                  {

                     if (m_syst_shift==e_Down)
                        {
                           if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*0.8873;
                           if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.6488;
                           if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.2792;
                           if (tau.pt() > 200) weight = weight*0.3535;
                        }
                     if (m_syst_shift==e_Up)
                        {
                           if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*1.1593;
                           if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.8951;
                           if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.7139;
                           if (tau.pt() > 200) weight = weight*1.8867;
                        }
                  }
            }
         if (m_channel == "SS_ClosureTest")
            {
              if (!m_tau_unc)
                  {
                     if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*1.04765;
                     if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.857812;
                     if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.59419;
                     if (tau.pt() > 200) weight = weight*1.06272;
                  } else
                  {

                     if (m_syst_shift==e_Down)
                        {
                           if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*0.8053;
                           if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.5939;
                           if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.264;
                           if (tau.pt() > 200) weight = weight*0.3257;
                        }
                     if (m_syst_shift==e_Up)
                        {
                           if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*1.2415;
                           if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.95;
                           if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.7287;
                           if (tau.pt() > 200) weight = weight*1.9066;
                        }
                  }
            }
         if (m_channel == "OS_ClosureTest")
            {
               if (!m_tau_unc)
                  {
                     if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*1.04765;
                     if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.857812;
                     if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.59419;
                     if (tau.pt() > 200) weight = weight*1.06272;
                  } else
                  {

                     if (m_syst_shift==e_Down)
                        {
                           if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*0.7563;
                           if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.5593;
                           if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.2522;
                           if (tau.pt() > 200) weight = weight*0.3535;
                        }
                     if (m_syst_shift==e_Up)
                        {
                           if (tau.pt() > 20 && tau.pt() <= 60) weight = weight*1.2906;
                           if (tau.pt() > 60 && tau.pt() <= 120) weight = weight*0.9846;
                           if (tau.pt() > 120 && tau.pt() <= 200) weight = weight*0.7401;
                           if (tau.pt() > 200) weight = weight*1.9225;
                        }
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

TopTaggingScaleFactors::TopTaggingScaleFactors(E_SystShift sys_toptag, E_SystShift sys_mistag)
{

    m_sys_toptag = sys_toptag;
    m_sys_mistag = sys_mistag;


    _scale_toptag = new ToptagScale();
    _eff_toptag = new ToptagEfficiency();
    _scale_topmistag = new TopMistagScale();
    _eff_topmistag = new TopMistagEfficiency();

}

double TopTaggingScaleFactors::GetWeight()
{
    static EventCalc* calc = EventCalc::Instance();
    BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

    std::vector< TopJet > *jets =  calc->GetCAJets();
    if(!jets) return 1.0;

    std::vector< GenParticle > *genparticles = calc->GetGenParticles();
    if(!genparticles) return 1.0;

    double scale_factor = 1.;

    bool toptagevent = false;
    for(unsigned int i=0; i<jets->size(); ++i) {

        TopJet jet = jets->at(i);

        double scale_jet = 1.0;
        double mmin=0;
        double mjet=0;
        int nsubjets=0;
        bool result = TopTag(jet,mjet,nsubjets,mmin);
        if(result)
            toptagevent = true;
        float jet_pt = jet.pt();
        float jet_eta = fabs(jet.eta());

        //Only apply corrections to high pT TopJets
        //Mistag measurement starts at 200, Efficiency starts at 400, but TopTag algorithm is rated for 350
        if(jet_pt < 400.0)
            continue;

        bool truetop = false;
        for(unsigned int j=0; j<genparticles->size(); ++j) {
            GenParticle p = genparticles->at(j);
            double deltaR = jet.deltaR(p);
            if (deltaR < 0.5 && abs(p.pdgId())==6) {
                bool leptonic_decay = false;
                const GenParticle* d1 = p.daughter(genparticles,1);
                const GenParticle* d2 = p.daughter(genparticles,2);
                const GenParticle* d11 = 0;
                const GenParticle* d12 = 0;
                if(abs(d1->pdgId())==24) {
                    d11 = d1->daughter(genparticles,1);
                    d12 = d1->daughter(genparticles,2);
                } else if(abs(d2->pdgId())==24) {
                    d11 = d2->daughter(genparticles,1);
                    d12 = d2->daughter(genparticles,2);
                }
                if(d11 && d12) {
                    if(abs(d11->pdgId())==11)
                        leptonic_decay = true;
                    else if(abs(d11->pdgId())==12)
                        leptonic_decay = true;
                    else if(abs(d11->pdgId())==13)
                        leptonic_decay = true;
                    else if(abs(d11->pdgId())==14)
                        leptonic_decay = true;
                    else if(abs(d11->pdgId())==15)
                        leptonic_decay = true;
                    else if(abs(d11->pdgId())==16)
                        leptonic_decay = true;
                    else if(abs(d12->pdgId())==11)
                        leptonic_decay = true;
                    else if(abs(d12->pdgId())==12)
                        leptonic_decay = true;
                    else if(abs(d12->pdgId())==13)
                        leptonic_decay = true;
                    else if(abs(d12->pdgId())==14)
                        leptonic_decay = true;
                    else if(abs(d12->pdgId())==15)
                        leptonic_decay = true;
                    else if(abs(d12->pdgId())==16)
                        leptonic_decay = true;
                }
                if(!leptonic_decay) {
                        truetop = true;
                }
            }
        }

        if(truetop) { // Found a top quark that could be parent of jet, so it is not possible to be mistagged. Apply toptag efficiency scale factor
            /*
            if(result) {
                std::cout << "Top-Tag event!" << std::endl;
            } else {
                std::cout << "Top-NoTag event!" << std::endl;
            }
            */
	    scale_jet = scale(result, jet_pt, jet_eta,
                              _scale_toptag, _eff_toptag,
                              m_sys_toptag);
        } else { // No true top quarks found near this jet, apply mistag scale factor
            /*
            if(result) {
                std::cout << "Mistag event!" << std::endl;
            } else {
                std::cout << "NoTop-NoTag event!" << std::endl;
            }
            */
	    scale_jet = scale(result, jet_pt, jet_eta,
                              _scale_topmistag, _eff_topmistag,
                              m_sys_mistag);
        }

        scale_factor *= scale_jet;
    }
    /*
    if(toptagevent)
        std::cout << "TopTag weight: " << scale_factor << std::endl;
    else
        std::cout << "NoTTag weight: " << scale_factor << std::endl;

    */

    return scale_factor;
}


float TopTaggingScaleFactors::scale(const bool &is_tagged,
                                  const float &jet_pt,
                                  const float &jet_eta,
                                  const ToptagFunction* sf,
                                  const ToptagFunction* eff,
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


BTaggingScaleFactors::BTaggingScaleFactors(
   E_BtagType btagtype, E_LeptonSelection lepton_selection, E_SystShift sys_bjets, E_SystShift sys_ljets, bool use_subjet_btags
)
{
    m_sys_bjets = sys_bjets;
    m_sys_ljets = sys_ljets;

    m_btagtype = btagtype;
    m_lepton_selection = lepton_selection;

    _scale_btag=new BtagScale(btagtype);
    _eff_btag=new BtagEfficiency(btagtype, lepton_selection);
    _eff_btag_subj=new BtagEfficiency(btagtype, lepton_selection, true, false);
    _eff_btag_topj=new BtagEfficiency(btagtype, lepton_selection, false, true);

    _scale_ctag=new CtagScale(btagtype);
    _eff_ctag=new CtagEfficiency(btagtype, lepton_selection);
    _eff_ctag_subj=new CtagEfficiency(btagtype, lepton_selection, true, false);
    _eff_ctag_topj=new CtagEfficiency(btagtype, lepton_selection, false, true);

    _scale_light=new LtagScale(btagtype);
    _eff_light=new LtagEfficiency(btagtype, lepton_selection);
    _eff_light_subj=new LtagEfficiency(btagtype, lepton_selection, true, false);
    _eff_light_topj=new LtagEfficiency(btagtype, lepton_selection, false, true);

    m_use_subjet_btags =use_subjet_btags;
}


double BTaggingScaleFactors::GetWeight()
{
    EventCalc* calc = EventCalc::Instance();

    std::vector< Jet > *jets =  calc->GetJets();
    if(!jets) return 1.0;

    double scale_factor = 1.;

    for(unsigned int i=0; i<jets->size(); ++i) {

        Jet jet = jets->at(i);

	//skip jets that overlap with top-jets
	bool overlap_with_topjet = false;
	for(unsigned int m = 0; m< calc->GetCAJets()->size();++m){
	  TopJet topjet = calc->GetCAJets()->at(m);
	  if(topjet.pt()>250. && topjet.deltaR(jets->at(i))<1.3) overlap_with_topjet = true;
	}
	if(m_use_subjet_btags && overlap_with_topjet) continue;

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

    if(m_use_subjet_btags){

      for(unsigned int m = 0; m< calc->GetCAJets()->size();++m){

	TopJet topjet = calc->GetCAJets()->at(m);
	if(topjet.pt() < 250.) continue;

	unsigned int nsubjets = topjet.subjets().size();
	double min_dr = 5.;
	if(nsubjets > 1){
	  // loop over subjets, find minimum DeltaR between them
	  for(unsigned int g = 0; g<nsubjets-1; ++g){

	    Particle subjetg = topjet.subjets().at(g);
	    for(unsigned int k = g+1; k<nsubjets; ++k){

	      double dr = subjetg.deltaR(topjet.subjets().at(k));
	      if(dr < min_dr) min_dr = dr;
	    }
	  }
	}


	//apply scale factor for individual subjets if min_dr>0.4
	if(min_dr>0.4){
	  double discriminator_cut= 0.0;
	  if(m_btagtype==e_CSVL) discriminator_cut = 0.244;
	  if(m_btagtype==e_CSVM) discriminator_cut = 0.679;
	  if(m_btagtype==e_CSVT) discriminator_cut = 0.898;

	  std::vector<float> btagsub_combinedSecondaryVertex_top;
	  std::vector<int> flavorsub_top;

	  btagsub_combinedSecondaryVertex_top=topjet.btagsub_combinedSecondaryVertex();
	  flavorsub_top=topjet.flavorsub();

	  for(unsigned int j = 0; j<nsubjets; ++j){
	    bool result = btagsub_combinedSecondaryVertex_top[j]>discriminator_cut;
	    double scale_jet = 1.0;
	    float jet_pt = topjet.subjets().at(j).pt();
	    float jet_eta = fabs(topjet.subjets().at(j).eta());

	    //do not consider jets outside b-tagging range
	    if(jet_eta>=2.4){
	      continue;
	    }
	    switch(abs(flavorsub_top[j])) {
	    case 5: // b-quark
	      scale_jet = scale(result, jet_pt, jet_eta,
				_scale_btag, _eff_btag_subj,
				m_sys_bjets);
            break;

	    case 4: // c-quark
	      scale_jet = scale(result, jet_pt, jet_eta,
				_scale_ctag, _eff_ctag_subj,
				m_sys_bjets);
	    break;

	    case 3: // s-quark
	    case 2: // d-quark
	    case 1: // u-quark
	    case 21: // gluon
	      scale_jet = scale(result, jet_pt, jet_eta,
				_scale_light, _eff_light_subj,
				m_sys_ljets);
            break;

	    default:
	      break;

	    }

	    scale_factor *= scale_jet;
	  }
	}
	//apply scale factor for entire topjet if min_dr<0.4
	else{
	  bool result = IsTagged(topjet, m_btagtype);
	  double scale_jet = 1.0;
	  float jet_pt = topjet.pt();
	  float jet_eta = fabs(topjet.eta());

	  //b-tagging for fat jets with eta>2.4 possible? -> use scale factors for eta=2.4
	  if(jet_eta>=2.4){
	    jet_eta=2.4;
	  }
	  switch(abs(topjet.flavor())) {
	  case 5: // b-quark
	    scale_jet = scale(result, jet_pt, jet_eta,
			      _scale_btag, _eff_btag_topj,
			      m_sys_bjets);
            break;

	  case 4: // c-quark
	    scale_jet = scale(result, jet_pt, jet_eta,
				_scale_ctag, _eff_ctag_topj,
			      m_sys_bjets);
	    break;

	  case 3: // s-quark
	  case 2: // d-quark
	  case 1: // u-quark
	  case 21: // gluon
	    scale_jet = scale(result, jet_pt, jet_eta,
			      _scale_light, _eff_light_topj,
			      m_sys_ljets);
            break;

	  default:
	    break;

	  }

	  scale_factor *= scale_jet;
	}
      }
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

// Toptag Scale
//
ToptagScale::ToptagScale() : ToptagFunction() {}


float ToptagScale::value(const float &jet_pt, const float &jet_eta) const
{
    //Flat value for fit in theta
    return 1.0;
    //From
    // https://cds.cern.ch/record/1647419/files/JME-13-007-pas.pdf
    //if(abs(jet_eta) < 1.0)
    //    return 0.985;
    //else
    //    return 0.644;
}


float ToptagScale::error(const float &jet_pt, const float &jet_eta) const
{
    //Flat error for fit in theta
    return 0.2;
    //From
    // https://cds.cern.ch/record/1647419/files/JME-13-007-pas.pdf
    //if(abs(jet_eta) < 1.0)
    //    return 0.073;
    //else
    //    return 0.100;
}



// TopMistag Scale
//
TopMistagScale::TopMistagScale() : ToptagFunction() {;}


float TopMistagScale::value(const float &jet_pt, const float &jet_eta) const
{
    //From Mistag Study
    return 0.83;
}


float TopMistagScale::error(const float &jet_pt, const float &jet_eta) const
{
    //From Mistag Study
    return 0.21;
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
    const float CSVMErrors[] = {
      0.0415707,
      0.0204209,
      0.0223227,
      0.0206655,
      0.0199325,
      0.0174121,
      0.0202332,
      0.0182446,
      0.0159777,
      0.0218531,
      0.0204688,
      0.0265191,
      0.0313175,
      0.0415417,
      0.0740446,
      0.0596716
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
    case e_CSVM:
        _scale = new TF1("csvmb", "(0.938887+(0.00017124*x))+(-2.76366e-07*(x*x))", 20.0, 800.0);
        _bins.assign(bins, bins + 17);
        _errors.assign(CSVMErrors, CSVMErrors + 16);
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
        _scale = new TF2("csvtl","((1.00462+(0.00325971*x))+(-7.79184e-06*(x*x)))+(5.22506e-09*(x*(x*x)))+y-y",20.,800.0,0,2.4);
        _scale_plus = new TF2("cvstl_plus","((1.16361+(0.00464695*x))+(-1.09467e-05*(x*x)))+(7.21896e-09*(x*(x*x)))+y-y",20.,800.0,0,2.4);
        _scale_minus = new TF2("cvstl_minus","((0.845757+(0.00186422*x))+(-4.6133e-06*(x*x)))+(3.21723e-09*(x*(x*x)))+y-y",20.,800.0,0,2.4);
        break;
    case e_CSVM:
        _scale = new TF2("csvml","(((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x))))*(y<0.8)+(((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x))))*(0.8<=y)*(y<1.6)+(((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x))))*(1.6<=y)*(y<2.4)",20.,800.0,0,2.4);
        _scale_plus = new TF2("cvsml_plus","(((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x))))*(0.8>y)+(((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x))))*(0.8<=y)*(y<1.6)+(((1.15575+(0.0006933444*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x))))*(1.6<=y)*(y<2.4)",20.,800.0,0,2.4);
        _scale_minus = new TF2("cvsml_minus","(((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x))))*(0.8>y)+(((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x))))*(0.8<=y)*(y<1.6)+(((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x))))*(1.6<=y)*(y<2.4)",20.,800.0,0,2.4);
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

// Toptag Efficiency
//
ToptagEfficiency::ToptagEfficiency() : ToptagFunction() {

}


float ToptagEfficiency::value(const float &jet_pt, const float &jet_eta) const
{
    //From Mistag study - P. Turner
    return 0.199;
}

// TopMistag Efficiency
//
TopMistagEfficiency::TopMistagEfficiency() : ToptagFunction() {

}


float TopMistagEfficiency::value(const float &jet_pt, const float &jet_eta) const
{
    //Measured in mistag study - P. Turner
    return 0.0152;
}


// Btag Efficiency
//
BtagEfficiency::BtagEfficiency(E_BtagType btagtype, E_LeptonSelection leptonsel, bool dosubjets, bool dotopjets) : BtagFunction(btagtype)
{
    const float bins[] = {
        20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 120.0, 160.0, 210.0, 260.0, 320.0, 400.0, 500.0, 600.0, 800.0, 1600.0
    };

    const float CSVTEfficiencies[] = {
      0, 0, 0, 0.516861 , 0.530606 , 0.555624 , 0.541633 , 0.539517 , 0.507292 , 0.424228 , 0.381748 , 0.333086 , 0.26585 , 0.207664 , 0.158442 , 0.135173 , 0.103829
    };
    const float CSVLEfficiencies[] = {
      0, 0, 0, 0.836914 , 0.831639 , 0.84178 , 0.84026 , 0.835832 , 0.8334 , 0.826531 , 0.808195 , 0.777609 , 0.746216 , 0.678185 , 0.666574 , 0.564218 , 0.521006
    };
    const float CSVTEfficiencies_SubJet[] = {
      0.283184 , 0.440065 , 0.528707 , 0.46884 , 0.485505 , 0.563361 , 0.52208 , 0.532883 , 0.509963 , 0.425465 , 0.390519 , 0.350492 , 0.283081 , 0.233869 , 0.104842 , 0.105325 , 0
    };
    const float CSVLEfficiencies_SubJet[] = {
      0.786544 , 0.827889 , 0.834176 , 0.791612 , 0.803744 , 0.833457 , 0.841063 , 0.833662 , 0.844785 , 0.830733 , 0.818157 , 0.785523 , 0.761376 , 0.707712 , 0.55595 , 0.469558 , 0.54417
    };
    const float CSVTEfficiencies_TopJet[] = {
      0,0,0,0,0,0,0,0,0,0, 0.269034 , 0.268951 , 0.245795 , 0.219509 , 0.158286 , 0.158518 , 0.100975
    };
    const float CSVLEfficiencies_TopJet[] = {
      0,0,0,0,0,0,0,0,0,0, 0.7533 , 0.749066 , 0.732398 , 0.710423 , 0.662259 , 0.600003 , 0.542833
    };
    const float CSVTEfficiencies_mu[] = {
      0, 0, 0, 0.527496 , 0.538724 , 0.555376 , 0.540102 , 0.528351 , 0.50122 , 0.41426 , 0.371213 , 0.311879 , 0.264041 , 0.215617 , 0.208605 , 0.150797 , 0.130312
    };
    const float CSVLEfficiencies_mu[] = {
      0, 0, 0, 0.843019 , 0.836389 , 0.847278 , 0.839844 , 0.833074 , 0.835171 , 0.820583 , 0.804997 , 0.762067 , 0.730621 , 0.707983 , 0.655583 , 0.607146 , 0.495463
    };
    const float CSVTEfficiencies_SubJet_mu[] = {
      0.23722 , 0.348922 , 0.433808 , 0.502618 , 0.533373 , 0.517228 , 0.477777 , 0.518199 , 0.505178 , 0.407392 , 0.370565 , 0.327397 , 0.277921 , 0.241946 , 0.166871 , 0.148731 , 0.000610035
    };
    const float CSVLEfficiencies_SubJet_mu[] = {
      0.801842 , 0.762325 , 0.820935 , 0.788438 , 0.813816 , 0.808677 , 0.842798 , 0.821366 , 0.839112 , 0.826879 , 0.810295 , 0.779161 , 0.732937 , 0.693156 , 0.699269 , 0.540767 , 0.194262
    };
    const float CSVTEfficiencies_TopJet_mu[] = {
      0,0,0,0,0,0,0,0,0,0, 0.321782 , 0.290763 , 0.252109 , 0.23444 , 0.192822 , 0.154837 , 0.114785
    };
    const float CSVLEfficiencies_TopJet_mu[] = {
      0,0,0,0,0,0,0,0,0,0, 0.790861 , 0.761854 , 0.73796 , 0.745972 , 0.653736 , 0.617015 , 0.537142
    };
    const float CSVMEfficiencies[] = {
      0, 0, 0, 0.662161, 0.681654, 0.704375, 0.718067, 0.730097, 0.719055, 0.699008, 0.663485, 0.622932, 0.566283, 0.47847, 0.446908, 0.36968, 0.342942
    };
    const float CSVMEfficiencies_SubJet[] = {
      0.456416, 0.557567, 0.687558, 0.621622, 0.729128, 0.743774, 0.720121, 0.742028, 0.726798, 0.709854, 0.671339, 0.609182, 0.543616, 0.449373, 0.210045, 0.217199, 0.457737
    };
    const float CSVMEfficiencies_TopJet[] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.493238, 0.447682, 0.386311, 0.328427
    };
    const float CSVMEfficiencies_mu[] = {
      0, 0, 0, 0.6596, 0.688346, 0.7073, 0.715599, 0.720623, 0.719989, 0.693183, 0.67499, 0.60586, 0.549044, 0.51809, 0.464164, 0.376361, 0.291773
    };
    const float CSVMEfficiencies_SubJet_mu[] = {
      0.355551, 0.481819, 0.617543, 0.723595, 0.755656, 0.728251, 0.671985, 0.754468, 0.740315, 0.707825, 0.666827, 0.623542, 0.585943, 0.519131, 0.464292, 0.394569, 0
    };
    const float CSVMEfficiencies_TopJet_mu[] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.483461, 0.449241, 0.394077, 0.288428
    };
    _bins.assign(bins, bins + 18);
    if (btagtype == e_CSVT && leptonsel == e_Electron && !dosubjets && !dotopjets) {
        _values.assign(CSVTEfficiencies, CSVTEfficiencies + 17);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Muon && !dosubjets && !dotopjets) {
      _values.assign(CSVTEfficiencies_mu, CSVTEfficiencies_mu + 17);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Electron && dosubjets && !dotopjets) {
      _values.assign(CSVTEfficiencies_SubJet, CSVTEfficiencies_SubJet + 17);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Electron && !dosubjets && dotopjets) {
      _values.assign(CSVTEfficiencies_TopJet, CSVTEfficiencies_TopJet + 17);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Muon && dosubjets && !dotopjets) {
      _values.assign(CSVTEfficiencies_SubJet_mu, CSVTEfficiencies_SubJet_mu + 17);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Muon && !dosubjets && dotopjets) {
      _values.assign(CSVTEfficiencies_TopJet_mu, CSVTEfficiencies_TopJet_mu + 17);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Electron && !dosubjets && !dotopjets) {
      _values.assign(CSVLEfficiencies, CSVLEfficiencies + 17);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Muon && !dosubjets && !dotopjets) {
      _values.assign(CSVLEfficiencies_mu, CSVLEfficiencies_mu + 17);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Electron && dosubjets && !dotopjets) {
      _values.assign(CSVLEfficiencies_SubJet, CSVLEfficiencies_SubJet + 17);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Electron && !dosubjets && dotopjets) {
      _values.assign(CSVLEfficiencies_TopJet, CSVLEfficiencies_TopJet + 17);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Muon && dosubjets && !dotopjets) {
      _values.assign(CSVLEfficiencies_SubJet_mu, CSVLEfficiencies_SubJet_mu + 17);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Muon && !dosubjets && dotopjets) {
      _values.assign(CSVLEfficiencies_TopJet_mu, CSVLEfficiencies_TopJet_mu + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Electron && !dosubjets && !dotopjets) {
      _values.assign(CSVMEfficiencies, CSVMEfficiencies + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Muon && !dosubjets && !dotopjets) {
      _values.assign(CSVMEfficiencies_mu, CSVMEfficiencies_mu + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Electron && dosubjets && !dotopjets) {
      _values.assign(CSVMEfficiencies_SubJet, CSVMEfficiencies_SubJet + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Electron && !dosubjets && dotopjets) {
      _values.assign(CSVMEfficiencies_TopJet, CSVMEfficiencies_TopJet + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Muon && dosubjets && !dotopjets) {
      _values.assign(CSVMEfficiencies_SubJet_mu, CSVMEfficiencies_SubJet_mu + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Muon && !dosubjets && dotopjets) {
      _values.assign(CSVMEfficiencies_TopJet_mu, CSVMEfficiencies_TopJet_mu + 17);
    }
    else {
        std::cerr <<  "unsupported b-tagging operating point and lepton selection in BtagEfficiency" <<std::endl;
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


CtagEfficiency::CtagEfficiency(E_BtagType btagtype, E_LeptonSelection leptonsel, bool dosubjets, bool dotopjets) : BtagEfficiency(btagtype, leptonsel, dosubjets, dotopjets)
{
    const float bins[] = {
        20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 120.0, 160.0, 210.0, 260.0, 320.0, 400.0, 1600.0
    };

    const float CSVTEfficiencies[] = {
      0, 0, 0, 0.0652631 , 0.0585445 , 0.0588226 , 0.0536955 , 0.0517468 , 0.0475731 , 0.0290347 , 0.0262906 , 0.0157309 , 0.0166553 , 0.0124923
    };
    const float CSVLEfficiencies[] = {
      0, 0, 0, 0.423634 , 0.394799 , 0.402879 , 0.390526 , 0.374143 , 0.379752 , 0.360757 , 0.347033 , 0.328483 , 0.300221 , 0.279772
    };
    const float CSVTEfficiencies_SubJet[] = {
      0.0631994 , 0.071359 , 0.0799394 , 0.0400092 , 0.0735092 , 0.0679709 , 0.0619224 , 0.0608981 , 0.0440748 , 0.0259931 , 0.0194072 , 0.0260485 , 0.0254116 , 0.0228735
    };
    const float CSVLEfficiencies_SubJet[] = {
      0.489766 , 0.552803 , 0.355127 , 0.449075 , 0.373065 , 0.411385 , 0.380089 , 0.36281 , 0.361701 , 0.361116 , 0.324198 , 0.362666 , 0.299875 , 0.354893
    };
    const float CSVTEfficiencies_TopJet[] = {
      0,0,0,0,0,0,0,0,0,0, 0.0315203 , 0.00788948 , 0.0169692 , 0.00868753
    };
    const float CSVLEfficiencies_TopJet[] = {
      0,0,0,0,0,0,0,0,0,0, 0.217208 , 0.30472 , 0.282458 , 0.269931
    };
    const float CSVTEfficiencies_mu[] = {
      0, 0, 0 , 0.0746531 , 0.0648892 , 0.0661782 , 0.0593803 , 0.0483994 , 0.0424287 , 0.0284464 , 0.0307131 , 0.0275881 , 0.0173931 , 0.0147546
    };
    const float CSVLEfficiencies_mu[] = {
      0, 0, 0 , 0.438554 , 0.409188 , 0.417346 , 0.417047 , 0.385195 , 0.376154 , 0.372287 , 0.342093 , 0.324201 , 0.312492 , 0.252123
    };
    const float CSVTEfficiencies_SubJet_mu[] = {
      0.0471672 , 0.0550835 , 0.0545521 , 0.0750428 , 0.0816817 , 0.0661709 , 0.0699742 , 0.0703206 , 0.0520946 , 0.039097 , 0.0267112 , 0.0331859 , 0.0151289 , 0.0262327
    };
    const float CSVLEfficiencies_SubJet_mu[] = {
      0.404398 , 0.513329 , 0.425533 , 0.403072 , 0.404578 , 0.457312 , 0.409479 , 0.406841 , 0.370388 , 0.361228 , 0.363352 , 0.375948 , 0.268852 , 0.27939
    };
    const float CSVTEfficiencies_TopJet_mu[] = {
      0,0,0,0,0,0,0,0,0,0, 0.0157415 , 0.0172926 , 0.0167701 , 0.0142135
    };
    const float CSVLEfficiencies_TopJet_mu[] = {
      0,0,0,0,0,0,0,0,0,0, 0.282709 , 0.320317 , 0.313382 , 0.241947
    };
    const float CSVMEfficiencies[] = {
      0, 0, 0, 0.186964, 0.194928, 0.193698, 0.199262, 0.207779, 0.206847, 0.178954, 0.178985, 0.153737, 0.136819, 0.11987
    };
    const float CSVMEfficiencies_SubJet[] = {
      0.18236, 0.191935, 0.194408, 0.176161, 0.235158, 0.142888, 0.184881, 0.18078, 0.173559, 0.151429, 0.173081, 0.191432, 0.110484, 0.192341
    };
    const float CSVMEfficiencies_TopJet[] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0915242
    };
    const float CSVMEfficiencies_mu[] = {
      0, 0, 0, 0.196396, 0.190891, 0.209847, 0.202461, 0.199787, 0.19889, 0.187892, 0.162749, 0.126961, 0.138817, 0.105963
    };
    const float CSVMEfficiencies_SubJet_mu[] = {
      0.0894353, 0.191412, 0.163019, 0.099625, 0.215634, 0.187806, 0.193702, 0.248561, 0.208723, 0.227753, 0.204565, 0.153551, 0.1491, 0.113878
    };
    const float CSVMEfficiencies_TopJet_mu[] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0909162
    };

    _bins.assign(bins, bins + 15);
    if (btagtype == e_CSVT && leptonsel == e_Electron && !dosubjets && !dotopjets) {
      _values.assign(CSVTEfficiencies, CSVTEfficiencies + 14);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Muon && !dosubjets && !dotopjets) {
      _values.assign(CSVTEfficiencies_mu, CSVTEfficiencies_mu + 14);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Electron && dosubjets && !dotopjets) {
      _values.assign(CSVTEfficiencies_SubJet, CSVTEfficiencies_SubJet + 14);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Electron && !dosubjets && dotopjets) {
      _values.assign(CSVTEfficiencies_TopJet, CSVTEfficiencies_TopJet + 14);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Muon && dosubjets && !dotopjets) {
      _values.assign(CSVTEfficiencies_SubJet_mu, CSVTEfficiencies_SubJet_mu + 14);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Muon && !dosubjets && dotopjets) {
      _values.assign(CSVTEfficiencies_TopJet_mu, CSVTEfficiencies_TopJet_mu + 14);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Electron && !dosubjets && !dotopjets) {
      _values.assign(CSVLEfficiencies, CSVLEfficiencies + 14);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Muon && !dosubjets && !dotopjets) {
      _values.assign(CSVLEfficiencies_mu, CSVLEfficiencies_mu + 14);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Electron && dosubjets && !dotopjets) {
      _values.assign(CSVLEfficiencies_SubJet, CSVLEfficiencies_SubJet + 14);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Electron && !dosubjets && dotopjets) {
      _values.assign(CSVLEfficiencies_TopJet, CSVLEfficiencies_TopJet + 14);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Muon && dosubjets && !dotopjets) {
      _values.assign(CSVLEfficiencies_SubJet_mu, CSVLEfficiencies_SubJet_mu + 14);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Muon && !dosubjets && dotopjets) {
      _values.assign(CSVLEfficiencies_TopJet_mu, CSVLEfficiencies_TopJet_mu + 14);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Electron && !dosubjets && !dotopjets) {
      _values.assign(CSVMEfficiencies, CSVMEfficiencies + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Muon && !dosubjets && !dotopjets) {
      _values.assign(CSVMEfficiencies_mu, CSVMEfficiencies_mu + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Electron && dosubjets && !dotopjets) {
      _values.assign(CSVMEfficiencies_SubJet, CSVMEfficiencies_SubJet + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Electron && !dosubjets && dotopjets) {
      _values.assign(CSVMEfficiencies_TopJet, CSVMEfficiencies_TopJet + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Muon && dosubjets && !dotopjets) {
      _values.assign(CSVMEfficiencies_SubJet_mu, CSVMEfficiencies_SubJet_mu + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Muon && !dosubjets && dotopjets) {
      _values.assign(CSVMEfficiencies_TopJet_mu, CSVMEfficiencies_TopJet_mu + 17);
    }
    else {
        std::cerr <<  "unsupported b-tagging operating point and lepton selection in CtagEfficiency" <<std::endl;
    }
}


// Light Efficiency
//
LtagEfficiency::LtagEfficiency(E_BtagType btagtype, E_LeptonSelection leptonsel, bool dosubjets, bool dotopjets) : BtagEfficiency(btagtype, leptonsel, dosubjets, dotopjets)
{
    // BtagFunction::BtagFunction(btagtype);

    const float bins[] = {
        20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 120.0, 160.0, 210.0, 260.0, 320.0, 400.0, 500.0, 600.0, 1600.0
    };

    const float CSVTEfficiencies[] = {
      0, 0, 0, 0.00175576 , 0.0013093 , 0.00122353 , 0.00100534 , 0.00119473 , 0.00120241 , 0.00078995 , 0.000833238 , 0.00128178 , 0.000682748 , 0.00160506 , 0.00251959 , 0.00057882
    };
    const float CSVLEfficiencies[] = {
      0, 0, 0, 0.126082 , 0.104035 , 0.105668 , 0.0991957 , 0.0958602 , 0.0892866 , 0.0963189 , 0.0936501 , 0.0909674 , 0.0896909 , 0.0896152 , 0.115579 , 0.102088
    };
    const float CSVTEfficiencies_SubJet[] = {
      0.00259388 , 0.000530493 , 0.00134186 , 0.00305683 , 0.00119375 , 0.00209844 , 0.00132584 , 0.000501977 , 0.0010974 , 0.00112209 , 0.0007239 , 0.001115 , 0.00137592 , 0.00218108 , 0.00162565 , 0
    };
    const float CSVLEfficiencies_SubJet[] = {
      0.215991 , 0.234865 , 0.129691 , 0.105705 , 0.0713213 , 0.0865309 , 0.0650489 , 0.0771265 , 0.0749646 , 0.0839248 , 0.0865248 , 0.0891748 , 0.0995316 , 0.0986173 , 0.137609 , 0.146648
    };
    const float CSVTEfficiencies_TopJet[] = {
      0,0,0,0,0,0,0,0,0,0, 0.000859636 , 0.000757145 , 0.00111742 , 0.000997218 , 0.0040979 , 0.000275211
    };
    const float CSVLEfficiencies_TopJet[] = {
      0,0,0,0,0,0,0,0,0,0, 0.0714881 , 0.0956987 , 0.0887327 , 0.0859661 , 0.0899636 , 0.0965541
    };
    const float CSVTEfficiencies_mu[] = {
      0, 0, 0, 0.00569484 , 0.00548976 , 0.00547509 , 0.00712279 , 0.00675918 , 0.00670482 , 0.00444024 , 0.0045339 , 0.00446893 , 0.00537538 , 0.00461259 , 0.00792868 , 0.00387412
    };
    const float CSVLEfficiencies_mu[] = {
      0, 0, 0, 0.13864 , 0.112884 , 0.110701 , 0.109742 , 0.0995396 , 0.099064 , 0.104746 , 0.101926 , 0.100509 , 0.104804 , 0.108763 , 0.104505 , 0.101038
    };
    const float CSVTEfficiencies_SubJet_mu[] = {
      0.00491591 , 0.00641215 , 0.00440554 , 0.00461913 , 0.00402397 , 0.00373321 , 0.0050391 , 0.00427545 , 0.00440775 , 0.00342722 , 0.00433205 , 0.00420399 , 0.00451348 , 0.00307578 , 0.00939786 , 0.00307568
    };
    const float CSVLEfficiencies_SubJet_mu[] = {
      0.215857 , 0.207552 , 0.128542 , 0.125484 , 0.0958364 , 0.0855152 , 0.0816815 , 0.0764793 , 0.0757557 , 0.0953726 , 0.098633 , 0.0983716 , 0.12283 , 0.114422 , 0.117415 , 0.14174
    };
    const float CSVTEfficiencies_TopJet_mu[] = {
      0,0,0,0,0,0,0,0,0,0, 0.00491419 , 0.00624411 , 0.00641003 , 0.00271659 , 0.00705914 , 0.00264482
    };
    const float CSVLEfficiencies_TopJet_mu[] = {
      0,0,0,0,0,0,0,0,0,0, 0.0924585 , 0.0921337 , 0.103826 , 0.090195 , 0.0914187 , 0.075329
    };
    const float CSVMEfficiencies[] = {
      0, 0, 0, 0.0104451, 0.0109593, 0.0099559, 0.0114463, 0.0128632, 0.0131994, 0.0131001, 0.0132028, 0.016202, 0.0150713, 0.0148206, 0.0266987, 0.0254222
    };
    const float CSVMEfficiencies_SubJet[] = {
      0.00538404, 0.0220808, 0.0101407, 0.0135111, 0.00608758, 0.00634972, 0.00464018, 0.0114969, 0.00766321, 0.0157621, 0.0176389, 0.0130717, 0.0168744, 0.0185147, 0.0159909, 0.0335146
    };
    const float CSVMEfficiencies_TopJet[] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0143229, 0.0223084, 0.0277529
    };
    const float CSVMEfficiencies_mu[] = {
      0, 0, 0, 0.0116975, 0.0104015, 0.0113221, 0.0120317, 0.0102773, 0.0134546, 0.0136533, 0.0131978, 0.0143784, 0.0164759, 0.0203827, 0.022571, 0.0111673
    };
    const float CSVMEfficiencies_SubJet_mu[] = {
      0.00999828, 0.0169453, 0.00321533, 0.00840078, 0.00773902, 0.00201981, 0.011657, 0.00475623, 0.00816656, 0.0138483, 0.0138901, 0.00889166, 0.0193149, 0.0262239, 0.0328485, 0.0207716
    };
    const float CSVMEfficiencies_TopJet_mu[] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.016739, 0.00905475, 0.00753405
    };

    _bins.assign(bins, bins + 17);
    if (btagtype == e_CSVT && leptonsel == e_Electron && !dosubjets && !dotopjets) {
        _values.assign(CSVTEfficiencies, CSVTEfficiencies + 16);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Muon && !dosubjets && !dotopjets) {
      _values.assign(CSVTEfficiencies_mu, CSVTEfficiencies_mu + 16);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Electron && dosubjets && !dotopjets) {
      _values.assign(CSVTEfficiencies_SubJet, CSVTEfficiencies_SubJet + 16);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Electron && !dosubjets && dotopjets) {
      _values.assign(CSVTEfficiencies_TopJet, CSVTEfficiencies_TopJet + 16);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Muon && dosubjets && !dotopjets) {
      _values.assign(CSVTEfficiencies_SubJet_mu, CSVTEfficiencies_SubJet_mu + 16);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Muon && !dosubjets && dotopjets) {
      _values.assign(CSVTEfficiencies_TopJet_mu, CSVTEfficiencies_TopJet_mu + 16);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Electron && !dosubjets && !dotopjets) {
      _values.assign(CSVLEfficiencies, CSVLEfficiencies + 16);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Muon && !dosubjets && !dotopjets) {
      _values.assign(CSVLEfficiencies_mu, CSVLEfficiencies_mu + 16);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Electron && dosubjets && !dotopjets) {
      _values.assign(CSVLEfficiencies_SubJet, CSVLEfficiencies_SubJet + 16);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Electron && !dosubjets && dotopjets) {
      _values.assign(CSVLEfficiencies_TopJet, CSVLEfficiencies_TopJet + 16);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Muon && dosubjets && !dotopjets) {
      _values.assign(CSVLEfficiencies_SubJet_mu, CSVLEfficiencies_SubJet_mu + 16);
    }
    else if (btagtype == e_CSVL && leptonsel == e_Muon && !dosubjets && dotopjets) {
      _values.assign(CSVLEfficiencies_TopJet_mu, CSVLEfficiencies_TopJet_mu + 16);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Electron && !dosubjets && !dotopjets) {
      _values.assign(CSVMEfficiencies, CSVMEfficiencies + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Muon && !dosubjets && !dotopjets) {
      _values.assign(CSVMEfficiencies_mu, CSVMEfficiencies_mu + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Electron && dosubjets && !dotopjets) {
      _values.assign(CSVMEfficiencies_SubJet, CSVMEfficiencies_SubJet + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Electron && !dosubjets && dotopjets) {
      _values.assign(CSVMEfficiencies_TopJet, CSVMEfficiencies_TopJet + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Muon && dosubjets && !dotopjets) {
      _values.assign(CSVMEfficiencies_SubJet_mu, CSVMEfficiencies_SubJet_mu + 17);
    }
    else if (btagtype == e_CSVM && leptonsel == e_Muon && !dosubjets && dotopjets) {
      _values.assign(CSVMEfficiencies_TopJet_mu, CSVMEfficiencies_TopJet_mu + 17);
    }
    else {
        std::cerr <<  "unsupported b-tagging operating point and lepton selection in LtagEfficiency" <<std::endl;
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



