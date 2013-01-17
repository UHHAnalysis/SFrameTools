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

    double Mu40triggerRunA_uncertainty[3] = {0.0013, 0.0042, 0.0034};
    double Mu40triggerRunB_uncertainty[3] = {0.0007, 0.0021, 0.0017};
    double Mu40triggerRunC_uncertainty[3] = {0.0004, 0.0014, 0.0011};

    double IsoMu24triggerRunA[3] = {0.9560, 0.9528, 0.9809};
    double IsoMu24triggerRunB[3] = {0.9798, 0.9618, 0.9814};
    double IsoMu24triggerRunC[3] = {0.9841, 0.9688, 1.0021};

    double IsoMu24triggerRunA_uncertainty[3] = {0.0008, 0.0021, 0.0016};
    double IsoMu24triggerRunB_uncertainty[3] = {0.0004, 0.0010, 0.0008};
    double IsoMu24triggerRunC_uncertainty[3] = {0.0003, 0.0009, 0.0007};

    double TightID_RunAB[3] = {0.9941, 0.9917, 0.9982};
    double TightID_RunC[3] = {0.9934, 0.9903, 0.9979};

    double TightID_RunAB_uncertainty[3] = {0.0003, 0.0005, 0.0004};
    double TightID_RunC_uncertainty[3] = {0.0003, 0.0005, 0.0003};

    double Isolation_RunAB[3] = {0.9923, 0.9979, 1.0019};
    double Isolation_RunC[3] = {0.9978, 1.0005, 1.0044};
    
    double Isolation_RunAB_uncertainty[3] = {0.0002, 0.0003, 0.0002};
    double Isolation_RunC_uncertainty[3] = {0.0002, 0.0003, 0.0002};
    
    if(m_syst_shift==e_Down){
      for(unsigned int i=0; i<3; ++i){
	Mu40triggerRunA[i] -= Mu40triggerRunA_uncertainty[i];
	Mu40triggerRunB[i] -= Mu40triggerRunB_uncertainty[i];
	Mu40triggerRunC[i] -= Mu40triggerRunC_uncertainty[i];
	IsoMu24triggerRunA[i] -= IsoMu24triggerRunA_uncertainty[i];
	IsoMu24triggerRunB[i] -= IsoMu24triggerRunB_uncertainty[i];
	IsoMu24triggerRunC[i] -= IsoMu24triggerRunC_uncertainty[i];
	TightID_RunAB[i] -= TightID_RunAB_uncertainty[i];
	TightID_RunC[i] -= TightID_RunC_uncertainty[i];
	Isolation_RunAB[i] -= Isolation_RunAB_uncertainty[i];
	Isolation_RunC[i] -= Isolation_RunC_uncertainty[i];
      }
    }
    if(m_syst_shift==e_Up){
      for(unsigned int i=0; i<3; ++i){
	Mu40triggerRunA[i] += Mu40triggerRunA_uncertainty[i];
	Mu40triggerRunB[i] += Mu40triggerRunB_uncertainty[i];
	Mu40triggerRunC[i] += Mu40triggerRunC_uncertainty[i];
	IsoMu24triggerRunA[i] += IsoMu24triggerRunA_uncertainty[i];
	IsoMu24triggerRunB[i] += IsoMu24triggerRunB_uncertainty[i];
	IsoMu24triggerRunC[i] += IsoMu24triggerRunC_uncertainty[i];
	TightID_RunAB[i] += TightID_RunAB_uncertainty[i];
	TightID_RunC[i] += TightID_RunC_uncertainty[i];
	Isolation_RunAB[i] += Isolation_RunAB_uncertainty[i];
	Isolation_RunC[i] += Isolation_RunC_uncertainty[i];
      }
    }


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

