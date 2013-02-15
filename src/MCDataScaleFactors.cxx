#include "include/MCDataScaleFactors.h"
#include <TMath.h>

LeptonScaleFactors::LeptonScaleFactors(std::vector<std::string> correctionlist)
{
    m_syst_shift = e_Default;
    m_muon_unc = false;
    m_ele_unc = false;
    m_tau_unc = false;
 
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
	  if(m_correctionlist[i].first == "MuonRunA") {
            m_mu_trig[etabin].push_back((TGraphAsymmErrors*) file->Get("2012A/SF_" + eta_name + "_TRIG_Mu40"));
            m_mu_id[etabin].push_back((TGraphAsymmErrors*) file->Get("2012A/SF_" + eta_name + "_ID_tight"));
	    isok = true;
	  } else if (m_correctionlist[i].first == "MuonRunB") {
            m_mu_trig[etabin].push_back((TGraphAsymmErrors*) file->Get("2012B/SF_" + eta_name + "_TRIG_Mu40"));
            m_mu_id[etabin].push_back((TGraphAsymmErrors*) file->Get("2012B/SF_" + eta_name + "_ID_tight"));
	    isok = true;
	  } else if (m_correctionlist[i].first == "MuonRunC") {
            m_mu_trig[etabin].push_back((TGraphAsymmErrors*) file->Get("2012C/SF_" + eta_name + "_TRIG_Mu40"));
            m_mu_id[etabin].push_back((TGraphAsymmErrors*) file->Get("2012C/SF_" + eta_name + "_ID_tight"));
	    isok = true;
	  } else if (m_correctionlist[i].first == "MuonRunD") {
            m_mu_trig[etabin].push_back((TGraphAsymmErrors*) file->Get("2012D/SF_" + eta_name + "_TRIG_Mu40"));
            m_mu_id[etabin].push_back((TGraphAsymmErrors*) file->Get("2012D/SF_" + eta_name + "_ID_tight"));	    
	    isok = true;
	  }
	  //isolated muons
	  else if(m_correctionlist[i].first == "IsoMuonRunA") {
            m_mu_trig[etabin].push_back((TGraphAsymmErrors*) file->Get("2012A/SF_" + eta_name + "_TRIG_IsoMu24"));
            m_mu_iso[etabin].push_back((TGraphAsymmErrors*) file->Get("2012A/SF_" + eta_name + "_ISO_tight"));
            m_mu_id[etabin].push_back((TGraphAsymmErrors*) file->Get("2012A/SF_" + eta_name + "_ID_tight"));
	    isok = true;
	  } else if (m_correctionlist[i].first == "IsoMuonRunB") {
            m_mu_trig[etabin].push_back((TGraphAsymmErrors*) file->Get("2012B/SF_" + eta_name + "_TRIG_IsoMu24"));
            m_mu_iso[etabin].push_back((TGraphAsymmErrors*) file->Get("2012B/SF_" + eta_name + "_ISO_tight"));
            m_mu_id[etabin].push_back((TGraphAsymmErrors*) file->Get("2012B/SF_" + eta_name + "_ID_tight"));
	    isok = true;
	  } else if (m_correctionlist[i].first == "IsoMuonRunC") {
            m_mu_trig[etabin].push_back((TGraphAsymmErrors*) file->Get("2012C/SF_" + eta_name + "_TRIG_IsoMu24"));
            m_mu_iso[etabin].push_back((TGraphAsymmErrors*) file->Get("2012C/SF_" + eta_name + "_ISO_tight"));
            m_mu_id[etabin].push_back((TGraphAsymmErrors*) file->Get("2012C/SF_" + eta_name + "_ID_tight"));
	    isok = true;	  
	  } else if (m_correctionlist[i].first == "IsoMuonRunD") {
            m_mu_trig[etabin].push_back((TGraphAsymmErrors*) file->Get("2012D/SF_" + eta_name + "_TRIG_IsoMu24"));
            m_mu_iso[etabin].push_back((TGraphAsymmErrors*) file->Get("2012D/SF_" + eta_name + "_ISO_tight"));
            m_mu_id[etabin].push_back((TGraphAsymmErrors*) file->Get("2012D/SF_" + eta_name + "_ID_tight"));
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

double LeptonScaleFactors::GetElectronWeight()
{
  if (!m_apply) return 1.;
  if (!IsUpToDate()){
    FillWeights();
  }
  double trig = GetElectronTrigWeight();
  return trig;
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

        switch(abs(JetFlavor(&jet))) {
        case 5: // b-quark
            scale_jet = scale(result, jet_pt,
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
            scale_jet = scale(result, jet_pt,
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
            scale_jet = scale(result, jet_pt,
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
                                  const BtagFunction* sf,
                                  const BtagFunction* eff,
                                  const E_SystShift &systematic)
{
    switch(systematic) {
    case e_Default:
        return is_tagged ?
               sf->value(jet_pt) :
               (1 - sf->value(jet_pt) * eff->value(jet_pt)) /
               (1 - eff->value(jet_pt));
        break;

    case e_Up:
        return is_tagged ?
               sf->value_plus(jet_pt) :
               (1 - sf->value_plus(jet_pt) * eff->value_plus(jet_pt)) /
               (1 - eff->value_plus(jet_pt));
        break;

    case e_Down:
        return is_tagged ?
               sf->value_minus(jet_pt) :
               (1 - sf->value_minus(jet_pt) * eff->value_minus(jet_pt)) /
               (1 - eff->value_minus(jet_pt));
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
    // Moriond13 prescription
    const float bins[] = {
        20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800
    };

    const float CSVTErrors[] = {
        0.0567059, 0.0266907, 0.0263491, 0.0342831,
        0.0303327, 0.024608, 0.0333786, 0.0317642,
        0.031102, 0.0295603, 0.0474663, 0.0503182,
        0.0580424, 0.0575776, 0.0769779, 0.0898199 
    };

    switch(btagtype) {
    case e_CSVT: // Moriond13 prescription
        _scale = new TF1("csvtb", "0.869965*((1.+(0.0335062*x))/(1.+(0.0304598*x)))", 20.0, 800.0);
        _bins.assign(bins, bins + 17);
        _errors.assign(CSVTErrors, CSVTErrors + 16);
        break;
    default:
        std::cerr <<  "unsupported b-tagging operating point" <<std::endl;
        break;
    }
}


const unsigned int BtagScale::find_bin(const float &jet_pt) const
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


float BtagScale::value(const float &jet_pt) const
{
    if (_scale->GetXmin() > jet_pt)
        return value(_scale->GetXmin());

    if (_scale->GetXmax() < jet_pt)
        return value(_scale->GetXmax());

    return _scale->Eval(jet_pt);
}


float BtagScale::error(const float &jet_pt) const
{
    if (_scale->GetXmin() > jet_pt)
        return 2.0 * error(_scale->GetXmin());

    if (_scale->GetXmax() < jet_pt)
        return 2.0 * error(_scale->GetXmax());

    return _errors.at(find_bin(jet_pt));
}


// Ctag scale
//
float CtagScale::error(const float &jet_pt) const
{
    return 2 * BtagScale::error(jet_pt);
}

// Light-tag scale
//
LtagScale::LtagScale(E_BtagType btagtype) : BtagFunction(btagtype)
{
    switch(btagtype) {
    case e_CSVT: // Moriond13 prescription
        _scale = new TF1("csvtl","((1.01739+(0.00283619*x))+(-7.93013e-06*(x*x)))+(5.97491e-09*(x*(x*x)))",20.,800.0);
        _scale_plus = new TF1("cvstl_plus","((1.08119+(0.00441909*x))+(-1.18764e-05*(x*x)))+(8.71372e-09*(x*(x*x)))",20.,800.0);
        _scale_minus = new TF1("cvstl_minus","((0.953587+(0.00124872*x))+(-3.97277e-06*(x*x)))+(3.23466e-09*(x*(x*x)))",20.,800.0);
        break;
    default:
        std::cerr <<  "unsupported b-tagging operating point" <<std::endl;
        break;
    }
}


float LtagScale::value(const float &jet_pt) const
{
    if (_scale->GetXmin() > jet_pt)
        return value(_scale->GetXmin());

    if (_scale->GetXmax() < jet_pt)
        return value(_scale->GetXmax());

    return _scale->Eval(jet_pt);
}


float LtagScale::value_plus(const float &jet_pt) const
{
    if (_scale_plus->GetXmin() > jet_pt)
    {
        double error = 2.0*(value_plus(_scale_plus->GetXmin()) - value(_scale_plus->GetXmin()));
        return value(_scale_plus->GetXmin()) + error;
    }

    if (_scale_plus->GetXmax() < jet_pt)
    {
        double error = 2.0*(value_plus(_scale_plus->GetXmax()) - value(_scale_plus->GetXmax()));
        return value(_scale_plus->GetXmax()) + error;
    }

    return _scale_plus->Eval(jet_pt);
}


float LtagScale::value_minus(const float &jet_pt) const
{
    if (_scale_minus->GetXmin() > jet_pt)
    {
        double error = 2.0*(value(_scale->GetXmin()) - value_minus(_scale_minus->GetXmin()));
        double scale = value(_scale_minus->GetXmin()) - error;
        if ( scale >= 0.0 ) return scale;
        return 0.0;
    }

    if (_scale_minus->GetXmax() < jet_pt)
    {
        double error = 2.0*(value(_scale->GetXmax()) - value_minus(_scale_minus->GetXmax()));
        double scale = value(_scale_minus->GetXmax()) - error;
        if ( scale >= 0.0 ) return scale;
        return 0.0;
    }

    return _scale_minus->Eval(jet_pt);
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
      0.48611233734744685, 0.49084495587364141, 0.49557757439983596, 0.50031019292603052, 
      0.50504281145222507, 0.51984588153656797, 0.5111155698055827, 0.50384994428241348,
      0.48720779281653653, 0.40203693185860434, 0.35119940680234535, 0.29170221794956114, 
      0.25324360996755907, 0.1908150859592869, 0.14997899317938629, 0.11390894204727826, 0.053563169122242862
    };

    if (btagtype == e_CSVT && leptonsel == e_Electron) { 
        _bins.assign(bins, bins + 18);
        _values.assign(CSVTEfficiencies, CSVTEfficiencies + 17);
    } 
    else if (btagtype == e_CSVT && leptonsel == e_Muon) { 
      _bins.assign(bins, bins + 18);
      _values.assign(CSVTEfficiencies_mu, CSVTEfficiencies_mu + 17);
    }
    else {
        std::cerr <<  "unsupported b-tagging operating point and lepton selection" <<std::endl;
    }
}


const unsigned int BtagEfficiency::find_bin(const float &jet_pt) const
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


float BtagEfficiency::value(const float &jet_pt) const
{
    return _values.at(find_bin(jet_pt));
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
      0.089953234452144898, 0.083116889562771884, 0.07628054467339887, 0.069444199784025856, 
      0.062607854894652842, 0.066249709912878596, 0.056876624468191909, 0.05002477475011223, 
      0.04508856575948126, 0.030786888568016747, 0.020212080758323085, 0.010082570271345626, 
      0.016672404844844568, 0.007150511494597455
    };

    if (btagtype == e_CSVT && leptonsel == e_Electron) {
        _bins.assign(bins, bins + 15);
        _values.assign(CSVTEfficiencies, CSVTEfficiencies + 14);
    }
    else if (btagtype == e_CSVT && leptonsel == e_Muon) { 
      _bins.assign(bins, bins + 15);
      _values.assign(CSVTEfficiencies_mu, CSVTEfficiencies_mu + 14);
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
      0.012368496488436461, 0.010695314162950109, 0.0090221318374637573, 0.0073489495119774045, 
      0.0056757671864910526, 0.0060836798818149195, 0.0070437434085952808, 0.0060739192548708524, 
      0.0065709163393043229, 0.0056167170501573507, 0.0061268017111299148, 0.0074465430391629462, 
      0.0070243350933687091, 0.0065102712066790842, 0.0042773193409681156, 0.0093149763709977663
    };

    if (btagtype == e_CSVT && leptonsel == e_Electron) {
        _bins.assign(bins, bins + 17);
        _values.assign(CSVTEfficiencies, CSVTEfficiencies + 16);
    } 
    else if (btagtype == e_CSVT && leptonsel == e_Muon) { 
      _bins.assign(bins, bins + 17);
      _values.assign(CSVTEfficiencies_mu, CSVTEfficiencies_mu + 16);
    }
    else {
        std::cerr <<  "unsupported b-tagging operating point and lepton selection" <<std::endl;
    }
}

