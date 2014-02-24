#ifndef MCDataScaleFactors_H
#define MCDataScaleFactors_H

#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include <TF1.h>
#include <TF2.h>

#include "SFrameTools/include/Utils.h"
#include "SFrameTools/include/EventCalc.h"

/**
 *  @short module to apply data-MC lepton scale factors for trigger and ID
 *
 *
 */
class LeptonScaleFactors {
public:
    /**
     * constructor
     *
     * first argument: list of corrections to be applied together with weight factors for luminosity, example: "MuonRunA 1.5 MuonRunB 2.6 MuonRunC 7.8"
     *
     * second argument: systematic shift
     * @see E_SystShift
     */
    LeptonScaleFactors(std::vector<std::string> correctionlist);
    ///Default destructor
    ~LeptonScaleFactors() {};

    ///browse configuration and fill all requested correction factors
    void FillWeights();

    ///return the total weight (muon*electron)
    ///derived from the weighted correction factor for muons (ID, trigger and isolation)
    ///and all electron weights
    double GetWeight();

    ///return the weighted correction factor for muon ID
    double GetMuonIDWeight();

    ///return the weighted correction factor for muon trigger
    double GetMuonTrigWeight();

    ///return the weighted correction factor for muon isolation
    double GetMuonIsoWeight();

    ///return the weighted correction factor for muons (ID, trigger and isolation)
    double GetMuonWeight();

    ///return the weighted correction factor for electron trigger
    double GetElectronTrigWeight();

    ///return the weighted correction factor for electron MVA ID
    double GetElectronMVAIDWeight();

    ///return the weighted correction factor for Ele30_OR_PFJet320 trigger
    double GetElectronORJetTrigWeight(const std::string& sys="none");

    ///return the weighted correction factor for electrons (right now: only trigger)
    double GetElectronWeight();

    ///check if the scale factors are up-to-date, fill them only once per run
    ///implemented for run-dependent scale factors
    bool IsUpToDate();
    
    /// return the bin number of the muon eta bin
    int GetMuonEtaBin(double eta);

    /// return the bin number of a graph corresponding to a certain value of the x coordinate
    int GetBin(double x, TGraphAsymmErrors* graph);
    
    void DoUpVarMuonSF(bool f=true){m_muon_unc=true; m_syst_shift=e_Up;}
    void DoDownVarMuonSF(bool f=true){m_muon_unc=true; m_syst_shift=e_Down;}

    void DoUpVarEleSF(bool f=true){m_ele_unc=true; m_syst_shift=e_Up;}
    void DoDownVarEleSF(bool f=true){m_ele_unc=true; m_syst_shift=e_Down;}

    void DoUpVarTauSF(bool f=true){m_tau_unc=true; m_syst_shift=e_Up;}
    void DoDownVarTauSF(bool f=true){m_tau_unc=true; m_syst_shift=e_Down;}
    
    void DoUpVarTauEleSF(bool f=true){m_tauele_unc=true; m_syst_shift=e_Up;}
    void DoDownVarTauEleSF(bool f=true){m_tauele_unc=true; m_syst_shift=e_Down;}

    void DoUpVarTauEffSF(bool f=true){m_tau_eff_unc=true; m_syst_shift=e_Up;}
    void DoDownVarTauEffSF(bool f=true){m_tau_eff_unc=true; m_syst_shift=e_Down;}
        /// return the scale factor for the fake rate of medium taus
   double GetTauWeight();

   double GetDecayModeFindingWeight();

   
    /// return the scale factor for the tau efficiency
    double GetTauEffUnc();


private:
    E_SystShift m_syst_shift;
    std::vector<std::pair<std::string, double> > m_correctionlist;
    bool m_apply;                   // should any scale factors be applied?
    bool m_muon_unc;                // do shift of muon scale factors 
    bool m_ele_unc;                 // do shift of electron scale factors 
    bool m_tau_unc;                 // do shift of tau scale factors 
    bool m_tauele_unc;               // do shift of e -> tau fake rate
    int m_current_run;              // run for which the scale factors are vali
    bool m_tau_eff_unc;             // do shift of tau efficiency scale factors 
   
    
    std::vector< std::vector<TGraphAsymmErrors*> > m_mu_id;    // two arrays: first index stands for eta bin, second for run period
    std::vector< std::vector<TGraphAsymmErrors*> > m_mu_trig;  // two arrays: first index stands for eta bin, second for run period
    std::vector< std::vector<TGraphAsymmErrors*> > m_mu_iso;   // two arrays: first index stands for eta bin, second for run period
    std::vector<double> m_weights;  // weights for different runs
    std::vector<double> m_ele_trig; // two-parameter function of relative isolation times additional weight
    TH2F* m_ele_mva; //2D histog for Egamma-POG SF for Electron-ID with TrigMVA
};


/**
 *  @short modules to apply data-MC btagging corrections
 *
 *
 */


class BtagFunction {
public:
    BtagFunction(E_BtagType btagtype) {
        m_btagtype = btagtype;
    }

    virtual ~BtagFunction() {
    };

    virtual float value(const float &x, const float &y) const = 0;
    virtual float value_plus(const float &x, const float &y) const = 0;
    virtual float value_minus(const float &x, const float &y) const = 0;

protected:
    E_BtagType m_btagtype;
};


class BtagScale: public BtagFunction {
public:

    BtagScale(E_BtagType);

    virtual float value(const float &jet_pt, const float &jet_eta) const;
    virtual float value_plus(const float &jet_pt, const float &jet_eta) const {
      return value(jet_pt, jet_eta) + error(jet_pt, jet_eta);
    }

    virtual float value_minus(const float &jet_pt, const float &jet_eta) const {
      const float value_ = value(jet_pt,jet_eta) - error(jet_pt,jet_eta);
        return value_ > 0 ? value_ : 0;
    }

protected:

    virtual float error(const float &jet_pt, const float &jet_eta) const;

private:

    TF1 * _scale;
    std::vector<float> _errors;
    std::vector<float> _bins;
    const unsigned int find_bin(const float &jet_pt, const float &jet_eta) const;
};


class CtagScale: public BtagScale {
public:

    CtagScale(E_BtagType btagtype) : BtagScale(btagtype) {}

protected:

    virtual float error(const float &jet_pt, const float &jet_eta) const;

};


class LtagScale: public BtagFunction {
public:

    LtagScale(E_BtagType btagtype);

    virtual float value(const float &jet_pt, const float &jet_eta ) const;
    virtual float value_plus(const float &jet_pt, const float &jet_eta) const;
    virtual float value_minus(const float &jet_pt, const float &jet_eta) const;

private:

    TF2 * _scale;
    TF2 * _scale_plus;
    TF2 * _scale_minus;

};


class BtagEfficiency: public BtagFunction {
public:

    BtagEfficiency(E_BtagType, E_LeptonSelection);

    virtual float value(const float &jet_p, const float &jet_etat) const;
    virtual float value_plus(const float &jet_pt, const float &jet_eta) const {
      return value(jet_pt,jet_eta);
    }

    virtual float value_minus(const float &jet_pt, const float &jet_eta) const {
      return value(jet_pt,jet_eta);
    }

protected:

    const unsigned int find_bin(const float &jet_pt, const float &jet_eta) const;
    std::vector<float> _values;
    std::vector<float> _bins;

};


class CtagEfficiency: public BtagEfficiency {
public:

    CtagEfficiency(E_BtagType, E_LeptonSelection);

};


class LtagEfficiency: public BtagEfficiency {
public:

    LtagEfficiency(E_BtagType, E_LeptonSelection);

};


/**
 *  @short module to apply data-MC scale factors for b tagging
 *
 *
 */
class BTaggingScaleFactors {
public:
    /**
     * constructor
     *
     * second argument: systematic shift
     * @see E_SystShift
     */
    BTaggingScaleFactors(
        E_BtagType, E_LeptonSelection, E_SystShift sys_bjets=e_Default, E_SystShift sys_ljets=e_Default
    );
    ///Default destructor
    ~BTaggingScaleFactors() {};

    ///return the weighted correction factor
    double GetWeight();

private:

    E_SystShift m_sys_bjets;
    E_SystShift m_sys_ljets;
    E_BtagType m_btagtype;
    E_LeptonSelection m_lepton_selection;

    float scale(const bool &is_tagged,
                const float &jet_pt,
                const float &jet_eta,
                const BtagFunction* sf,
                const BtagFunction* eff,
                const E_SystShift &sytematic);

    float scale_data(const bool &is_tagged,
                     const float &jet_pt,
		     const float &jet_eta,
                     const BtagFunction* sf,
                     const BtagFunction* eff,
                     const E_SystShift &sytematic);

    BtagFunction* _scale_btag;
    BtagFunction* _eff_btag;

    BtagFunction* _scale_ctag;
    BtagFunction* _eff_ctag;

    BtagFunction* _scale_light;
    BtagFunction* _eff_light;
};

class JetpTReweightingInWJets {
public:
   
   JetpTReweightingInWJets();
   ///Default destructor
   ~JetpTReweightingInWJets() {};
   
   ///return the weighted correction factor
   double GetWeight();
   
   void DoUpVarJetSF(bool f=true){m_jetpTreweigting_unc=true; m_syst_shift=e_Up;}
   void DoDownVarJetSF(bool f=true){m_jetpTreweigting_unc=true; m_syst_shift=e_Down;}
   
   
private:
   E_SystShift m_syst_shift;
   bool m_jetpTreweigting_unc; // do shift of jet pT reweighting in W+jets
};




#endif
