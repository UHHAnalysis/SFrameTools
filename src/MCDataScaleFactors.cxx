#include "include/MCDataScaleFactors.h"

LeptonScaleFactors::LeptonScaleFactors(std::vector<std::string> correctionlist, E_SystShift syst_shift)
{
    m_syst_shift = syst_shift;

    if(correctionlist.size()%2!=0) {
        std::cerr<< "not a valid list of correction factors given to LeptonScaleFactors" <<std::endl;
    }

    for(unsigned int i=0; i< correctionlist.size()/2; ++i) {
        std::pair<std::string, double> correction (correctionlist[2*i], atof(correctionlist[2*i+1].c_str()));
        //std::cout << "Apply correction " << correction.first << " with factor " << correction.second <<std::endl;
        m_correctionlist.push_back(correction);
    }
}


double LeptonScaleFactors::GetWeight()
{
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

    if(!primlep) {
        std::cout << "WARNING: no primary lepton found in LeptonScaleFactors; return scale factor=1" <<std::endl;
        return 1.;
    }

    double eta = primlep->eta();

    int etabin=0;
    if(fabs(eta)<0.9) etabin=0;
    else if(fabs(eta)>=0.9 && fabs(eta)<1.2) etabin=1;
    else if(fabs(eta)>=1.2) etabin=2;


    for(unsigned int i=0; i<m_correctionlist.size(); ++i) {

        double weight = m_correctionlist[i].second;

        //non isolated muons
        if(m_correctionlist[i].first == "MuonRunA") {
            triggerfactor += weight*Mu40triggerRunA[etabin];
            isofactor+=weight;
            IDfactor+= weight*TightID_RunAB[etabin];
        } else if (m_correctionlist[i].first == "MuonRunB") {
            triggerfactor += weight*Mu40triggerRunB[etabin];
            isofactor+=weight;
            IDfactor+= weight*TightID_RunAB[etabin];
        } else if (m_correctionlist[i].first == "MuonRunC") {
            triggerfactor += weight*Mu40triggerRunC[etabin];
            isofactor+=weight;
            IDfactor+= weight*TightID_RunC[etabin];
        }
        //isolated muons
        else if(m_correctionlist[i].first == "IsoMuonRunA") {
            triggerfactor += weight*IsoMu24triggerRunA[etabin];
            isofactor+=weight*Isolation_RunAB[etabin];
            IDfactor+= weight*TightID_RunAB[etabin];
        } else if (m_correctionlist[i].first == "IsoMuonRunB") {
            triggerfactor += weight*IsoMu24triggerRunB[etabin];
            isofactor+=weight*Isolation_RunAB[etabin];
            IDfactor+= weight*TightID_RunAB[etabin];
        } else if (m_correctionlist[i].first == "IsoMuonRunC") {
            triggerfactor += weight*IsoMu24triggerRunC[etabin];
            isofactor+=weight*Isolation_RunC[etabin];
            IDfactor+= weight*TightID_RunC[etabin];
        }

        else {
            std::cerr<< "No information found for lepton correction named " << m_correctionlist[i].first <<std::endl;
        }

        sumofweights += m_correctionlist[i].second;
    }

    triggerfactor/=sumofweights;
    IDfactor/=sumofweights;
    isofactor/=sumofweights;

    return triggerfactor*IDfactor*isofactor;
}


BTaggingScaleFactors::BTaggingScaleFactors(E_BtagType btagtype, E_LeptonSelection lepton_selection, E_SystShift syst_shift)
{
    m_syst_shift = syst_shift;
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
                              m_syst_shift);
            // std::cout << "b jet pt: " << jet_pt << " is tagged: " << result << " scale: "; 
            // std::cout << _scale_btag->value(jet_pt) << " eff: " << _eff_btag->value(jet_pt); 
            // std::cout << " weight: " << scale_jet << std::endl;
            break;

        case 4: // c-quark
            scale_jet = scale(result, jet_pt,
                              _scale_ctag, _eff_ctag,
                              m_syst_shift);
            // std::cout << "c jet pt: " << jet_pt << " is tagged: " << result << " scale: ";
            // std::cout << _scale_ctag->value(jet_pt) << " eff: " << _eff_ctag->value(jet_pt);
            // std::cout << " weight: " << scale_jet << std::endl;
            break;

        case 3: // s-quark
        case 2: // d-quark
        case 1: // u-quark
        case 21: // gluon
            scale_jet = scale(result, jet_pt,
                              _scale_light, _eff_light,
                              m_syst_shift);
            // std::cout << "l jet pt: " << jet_pt << " is tagged: " << result << " scale: ";
            // std::cout << _scale_light->value(jet_pt) << " eff: " << _eff_light->value(jet_pt);
            // std::cout << " weight: " << scale_jet << std::endl;
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

    if (_scale->GetXmax() <= jet_pt)
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
      0.46025286139927302, 0.47147453428056812, 0.48269620716186323, 0.49391788004315834, 
      0.50513955292445345, 0.5275641768123287, 0.52611555537617938, 0.52000877808816592, 
      0.49670607517271714, 0.41879584426638333, 0.35458829735036534, 0.31980464581467161, 
      0.2633065361800408, 0.15395244857264473, 0.14249076770260585, 0.12832505591928584, 
      0.088225347145104577
    };

    const float CSVTEfficiencies_mu[] = {
      0.4867069914961496, 0.49136872939368514, 0.49603046729122069, 0.50069220518875623, 
      0.50535394308629178, 0.51957801352315347, 0.51118159612471636, 0.50374865320213236, 
      0.48772688493061067, 0.40251157007918764, 0.3518993032067595, 0.29256013153569699, 
      0.25184289481885574, 0.19547194395919484, 0.14962567054357276, 0.11409543692961853, 0.05264742299938039
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
      0.11652306522527095, 0.10250891105181067, 0.0884947568783504, 0.074480602704890125, 
      0.060466448531429851, 0.057013303055346715, 0.053889756298672861, 0.054130869567084745, 
      0.045355750929470529, 0.0295676057229264, 0.029790827791783315, 0.023049530534977217, 
      0.011306632506145586, 0.0030160132589408047
    };

    const float CSVTEfficiencies_mu[] = {
      0.088776771803662724, 0.082308880623651492, 0.075840989443640261, 0.069373098263629029, 
      0.062905207083617798, 0.066150549199851694, 0.056731860275576362, 0.049512588212936318, 
      0.045087707460782625, 0.031411137335074765, 0.020572577595922629, 0.010437956331481727, 
      0.016478496399157043, 0.0070183175039210972
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
      0.00046861311173297324, 0.0017309639668615048, 0.0029933148219900355, 0.0042556656771185666, 
      0.0055180165322470977, 0.0058001789816431223, 0.0043355201565048731, 0.0067166893283312941, 
      0.007193842044811505, 0.0059138053063260108, 0.0051343925902238392, 0.0046618713210304329,
      0.0040317669576414765, 0.0077667262632723812, 0.016098840685614232, 0.0075174218292724223     
    };

    const float CSVTEfficiencies_mu[] = {
      0.013064488336535481, 0.011226378133306248, 0.0093882679300770159, 0.0075501577268477834, 
      0.0057120475236185509, 0.0061491974531369346, 0.0071655856827267305, 0.0063197261605225487, 
      0.006601452080494047, 0.005620369810170083, 0.0063686240585956438, 0.007323803773992962, 
      0.0069554450956424351, 0.0064922032240628327, 0.0042016735455371211, 0.0098530853530043802
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

