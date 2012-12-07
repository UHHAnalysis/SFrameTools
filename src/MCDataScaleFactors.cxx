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

  double IsoMu24triggerRunA[3] = {0.9560, 0.9528, 0.9809};
  double IsoMu24triggerRunB[3] = {0.9798, 0.9618, 0.9814};
  double IsoMu24triggerRunC[3] = {0.9841, 0.9688, 1.0021};

  double TightID_RunAB[3] = {0.9941, 0.9917, 0.9982};
  double TightID_RunC[3] = {0.9934, 0.9903, 0.9979};

  double Isolation_RunAB[3] = {0.9923, 0.9979, 1.0019};
  double Isolation_RunC[3] = {0.9978, 1.0005, 1.0044};

  EventCalc* calc = EventCalc::Instance();

  //Take lepton with highest transverse momentum as reference. Make sure that routine is called after selection of exactly one good lepton.
  Particle* primlep = calc->GetPrimaryLepton();

  if(!primlep){
    std::cout << "WARNING: no primary lepton found in LeptonScaleFactors; return scale factor=1" <<std::endl;
    return 1.;
  }
  
  double eta = primlep->eta();

  int etabin=0;
  if(fabs(eta)<0.9) etabin=0;
  else if(fabs(eta)>=0.9 && fabs(eta)<1.2) etabin=1;
  else if(fabs(eta)>=1.2) etabin=2;


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
      triggerfactor += weight*IsoMu24triggerRunA[etabin];
      isofactor+=weight*Isolation_RunAB[etabin];
      IDfactor+= weight*TightID_RunAB[etabin];   
    }
    else if (m_correctionlist[i].first == "IsoMuonRunB"){
      triggerfactor += weight*IsoMu24triggerRunB[etabin];
      isofactor+=weight*Isolation_RunAB[etabin];
      IDfactor+= weight*TightID_RunAB[etabin];
    }
    else if (m_correctionlist[i].first == "IsoMuonRunC"){
      triggerfactor += weight*IsoMu24triggerRunC[etabin];
      isofactor+=weight*Isolation_RunC[etabin];
      IDfactor+= weight*TightID_RunC[etabin];
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

BTaggingScaleFactors::BTaggingScaleFactors(E_SystShift syst_shift){

  m_syst_shift = syst_shift;
  

  _scale_btag=new BtagScale();
  _eff_btag=new BtagEfficiency();
  
  _scale_ctag=new CtagScale();
  _eff_ctag=new CtagEfficiency();
  
  _scale_light=new LightScale();
  _eff_light=new LightEfficiencyData();
  
}


double BTaggingScaleFactors::GetWeight(){

  EventCalc* calc = EventCalc::Instance();

  std::vector< Jet > *jets =  calc->GetJets();
  if(!jets) return 1;

  double scale_factor = 1.;

  for(unsigned int i=0; i<jets->size(); ++i){

    Jet jet = jets->at(i);
    
    bool result = jet.btag_combinedSecondaryVertex()>0.898;
    double scale_jet = 1;
    
    //if (jet.has_flavor())
      {
	const float jet_pt = jet.pt();
	switch(abs(JetFlavor(&jet)))
	  {
	  case 5: // b-quark
	    scale_jet = scale(result, jet_pt,
			   _scale_btag, _eff_btag,
			   m_syst_shift);
	    break;
	    
	  case 4: // c-quark
	    scale_jet = scale(result, jet_pt,
			       _scale_ctag, _eff_ctag,
			   m_syst_shift);
	    break;
	    
	  case 3: // s-quark
	  case 2: // d-quark
	  case 1: // u-quark
	  case 21: // gluon
		scale_jet = scale_data(result, jet_pt,
				    _scale_light, _eff_light,
				    m_syst_shift);
		break;
		
	  default:
	    break;
	  }
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
    switch(systematic)
    {
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

float BTaggingScaleFactors::scale_data(const bool &is_tagged,
                       const float &jet_pt,
                       const BtagFunction* sf,
                       const BtagFunction* eff,
                       const E_SystShift &systematic)
{
    switch(systematic)
    {
        case e_Default:
            return is_tagged ?
                   sf->value(jet_pt) :
                   (1 - eff->value(jet_pt)) /
                       (1 - eff->value(jet_pt) / sf->value(jet_pt));
            break;

        case e_Up:
            return is_tagged ?
                   sf->value_plus(jet_pt) :
                   (1 - eff->value_plus(jet_pt)) /
                       (1 - eff->value_plus(jet_pt) / sf->value_plus(jet_pt));
            break;

        case e_Down:
            return is_tagged ?
                   sf->value_minus(jet_pt) :
                   (1 - eff->value_minus(jet_pt)) /
                       (1 - eff->value_minus(jet_pt) / sf->value_minus(jet_pt));
            break;

        default:
	    std::cerr <<  "unsupported systematic" <<std::endl;

            break;
    }
    return 1.;
}


// BtagFunction
//
BtagFunction::BtagFunction()
{
    const float bins[] = {
        30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 670
    };

    _bins.assign(bins, bins + 15);
}

const unsigned int BtagFunction::find_bin(const float &jet_pt) const
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


// Btag Scale
//
BtagScale::BtagScale()
{
    const float errors[] = {
        0.0364717, 0.0362281, 0.0232876, 0.0249618, 0.0261482, 0.0290466,
        0.0300033, 0.0453252, 0.0685143, 0.0653621, 0.0712586, 0.094589,
        0.0777011, 0.0866563
    };

    _errors.assign(errors, errors + 14);
}

float BtagScale::value(const float &jet_pt) const
{
    // underflow
    //
    if (30 > jet_pt)
        return value(30);

    // overflow
    //
    if (670 < jet_pt)
        return value(670);

    return 0.901615 *
           (1. + 0.552628 * jet_pt) /
           (1. + 0.547195 * jet_pt);
}

float BtagScale::error(const float &jet_pt) const
{
    if (30 > jet_pt)
        return 0.12;

    if (670 <= jet_pt)
        return 2 * error(669);

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
float LightScale::value(const float &jet_pt) const
{
    if (20 > jet_pt)
        return value(20);

    if (670 < jet_pt)
        return value(670);

    return 0.948463 +
           0.00288102 * jet_pt -
           7.98091e-06 * pow(jet_pt, 2) +
           5.50157e-09 * pow(jet_pt, 3);
}

float LightScale::value_plus(const float &jet_pt) const
{
    if (20 > jet_pt)
        return value_plus(20);

    if (670 < jet_pt)
        return value_plus(670);

    return 0.997077 +
           0.00473953 * jet_pt -
           1.34985e-05 * pow(jet_pt, 2) +
           1.0032e-08 * pow(jet_pt, 3);
}

float LightScale::value_minus(const float &jet_pt) const
{
    if (20 > jet_pt)
        return value_minus(20);

    if (670 < jet_pt)
        return value_minus(670);

    return 0.899715 +
           0.00102278 * jet_pt -
           2.46335e-06 * pow(jet_pt, 2) +
           9.71143e-10 * pow(jet_pt, 3);
}



// Btag Efficiency
//
BtagEfficiency::BtagEfficiency()
{
    const float values[] = {
        0.0, 0.0, 0.45700194476, 0.481780083914, 0.485380425249, 0.515662731626,
        0.498871069179, 0.492350232848, 0.397570152294, 0.32058194111,
        0.271581953854, 0.224112593547, 0.11042330955, 0.123300043702
    };

    _values.assign(values, values + 14);
}

float BtagEfficiency::value(const float &jet_pt) const
{
    return _values.at(find_bin(jet_pt));
}



// Light Efficiency
//
LightEfficiency::LightEfficiency()
{
    const float values[] = {
        0.0, 0.0, 0.00665025786822, 0.00661178343117, 0.00383612052866,
        0.0096833532719, 0.00611911636857, 0.00985837381435, 0.00888843910838,
        0.0152255871887, 0.00625901203913, 0.00862685315311, 0.00729571108855,
        0.01859141652
    };

    _values.assign(values, values + 14);
}

float LightEfficiency::value(const float &jet_pt) const
{
    return _values.at(find_bin(jet_pt));
}



// Light Efficiency from Data
//
float LightEfficiencyData::value(const float &jet_pt) const
{
    if (30 > jet_pt)
        return value(30);

    if (670 < jet_pt)
        return value(670);

    return 0.00315116 * (1 -
                         0.00769281 * jet_pt +
                         2.58066e-05 * pow(jet_pt, 2) -
                         2.02149e-08 * pow(jet_pt, 3));
}



