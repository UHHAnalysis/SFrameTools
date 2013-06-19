#ifndef Utils_H
#define Utils_H

#include "Objects.h"
#include "BaseCycleContainer.h"
#include "TVector3.h"
#include <limits>
#include <algorithm>
#include "EventCalc.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVectorXYZE;

struct HigherPt {
    bool operator() (const Particle& j1, const Particle& j2) const {
        return j1.pt() > j2.pt();
    };
};

/**
 * list of lepton selection
 * @see  Cleaner::JetEnergyResolutionShifter
 */
enum E_LeptonSelection {
    e_Electron, /**< Electron selection */
    e_Muon, /**< Muon selection */
};

/**
 * list of systematic shift options
 * @see Cleaner::JetEnergyResolutionShifter 
 * @see Cleaner::JetLeptonSubtractor
*/
enum E_SystShift {
    e_Default, /**< default, no systematic shift is applied */
    e_Up, /**< up variation */
    e_Down /**< down variation */
};

/**
 * list of possible systematic uncertainties
*/
enum E_SysUnc {
  e_None=0,  /**< no systematic shift is applied */
  e_JEC,     /**< jet energy scale uncertainty */
  e_JER,     /**< jet resolution uncertainty */
  e_MuonSF,  /**< muon scale factor uncertainty */
  e_EleSF,   /**< electron scale factor uncertainty */
  e_TauSF,   /**< tau scale factor uncertainty */
  e_PDF,     /**< PDF uncertainty */
};

/**
 * list of supported b-tagging working points
 * @see  NBTagSelection
*/
enum E_BtagType {
    e_CSVT, /**< Combined Secondary Vertex tagger, tight working point */
    e_CSVM, /**< Combined Secondary Vertex tagger, medium working point */
    e_CSVL, /**< Combined Secondary Vertex tagger, loose working point */
    e_JPT,  /**< Jet Probability tagger, tight working point */
    e_JPM, /**< Jet Probability tagger, medium working point */
    e_JPL /**< Jet Probability tagger, loose working point */
};

enum E_EventFlavor {
    e_BFlavor, 
    e_CFlavor, 
    e_LFlavor
};


int subJetBTag(TopJet topjet, E_BtagType type);

bool HiggsTag(TopJet topjet, E_BtagType type1, E_BtagType type2);

bool HepTopTagFull(TopJet topjet);

/// return if a jet is tagged given a btag operating point
bool IsTagged(Jet &, E_BtagType);


/**
 * basic top tagging routine
 *
 * Test, if the jet topjet fulfills the top-tagging criteria. The tagging variables mjet, nsubjets and mmin are determined and returned.
 * @see NTopTagSelection
 */

bool variableHepTopTag(TopJet topjet, double ptJetMin = 200., double massWindowLower = 0.85, double massWindowUpper = 1.15, double cutCondition2 = 0.35, double cutCondition3 = 0.35);

bool HepTopTag(TopJet topjet);
bool TopTag(TopJet topjet, double& mjet, int& nsubjets, double& mmin);
bool variableTopTag(TopJet topjet, double &mjet, int &nsubjets, double &mmin, double mminLower = 50., double mjetLower = 140., double mjetUpper = 250.);

/**
 * W tagging routine
 *
 * @see NWTagSelection
 */
bool WTag(TopJet prunedjet, double& mjet, int &nsubjets, double& massdrop);

//double HTlep(const BaseCycleContainer *bcc); ->moved to EventCalc

/// returns the next jet in the list jets to the particle p in the eta-phi plane
Jet* nextJet(const Particle *p, std::vector<Jet> *jets);
/// relative momentum of the particle p with respect to the next jet in the list of jets
double pTrel(const Particle *p, std::vector<Jet> *jets);
/// minimal distance in the eta-phi plane of particle p with respect to the closest jet from the list of jets
double deltaRmin(const Particle *p, std::vector<Jet> *jets);

/// transfroms LorentzVector to TVector3
TVector3 toVector(LorentzVector v4);

/// transfroms LorentzVector to TVector3
TVector3 toVector(LorentzVectorXYZE v4);

/// converts LorentzVector into LorentzVector
LorentzVectorXYZE toXYZ(LorentzVector v4);

/// converts LorentzVector into LorentzVector
LorentzVector toPtEtaPhi(LorentzVectorXYZE v4);

double deltaR(LorentzVector v1, LorentzVector v2);

double double_infinity();
int int_infinity();

/// x^p
int myPow(int x, unsigned int p) ;

#endif
