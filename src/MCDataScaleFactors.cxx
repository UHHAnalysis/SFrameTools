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

    std::cout << std::endl;

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
        0.44271619165939158, 0.45804384242391016, 0.47337149318842875, 0.48869914395294733, 
        0.50402679471746592, 0.52842961125817112, 0.52339654337656272, 0.51921865895984687, 
        0.49630643341766967, 0.41637903949954108, 0.35762424281015814, 0.31724150562041209,
        0.26751153919597787, 0.16593378465313466, 0.1418335386210324, 0.12439196581392559, 
        0.081187881736808146
    };

    if (btagtype == e_CSVT && leptonsel == e_Electron) { 
        _bins.assign(bins, bins + 18);
        _values.assign(CSVTEfficiencies, CSVTEfficiencies + 17);
    } else {
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
        0.11223218468091323, 0.099345787849841605, 0.086459391018769982, 0.07357299418769836, 
        0.060686597356626737, 0.058987756474951504, 0.05621735538747516, 0.057809589244733527, 
        0.045126899932038556, 0.032870074129705336, 0.028739293547418936, 0.020226473451573173, 
        0.010572298652120057, 0.0028478785823845172
    };

    if (btagtype == e_CSVT && leptonsel == e_Electron) {
        _bins.assign(bins, bins + 15);
        _values.assign(CSVTEfficiencies, CSVTEfficiencies + 14);
    } else {
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
        0.0015019802440593231, 0.0025156090090614237, 0.0035292377740635244, 0.004542866539065625, 
        0.0055564953040677257, 0.0050704875352594421, 0.0044963051359018982, 0.0065689093400694348, 
        0.0065912786708248417, 0.0054336583407874599, 0.0052337972043546791, 0.0048587771285205621,
        0.0040475265394415792, 0.0069754963943397738, 0.017361524530729457, 0.00667538501112912
    };

    if (btagtype == e_CSVT && leptonsel == e_Electron) {
        _bins.assign(bins, bins + 17);
        _values.assign(CSVTEfficiencies, CSVTEfficiencies + 16);
    } else {
        std::cerr <<  "unsupported b-tagging operating point and lepton selection" <<std::endl;
    }
}

