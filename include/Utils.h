#ifndef Utils_H
#define Utils_H

#include "SFrameTools/include/Objects.h"
#include "SFrameTools/include/fwd.h"
#include "SFrameTools/include/boost_includes.h" // for shared_array

#include "TVector3.h"
#include <limits>
#include <algorithm>
#include <memory>
#include <TF1.h>


#define DEPRECATED __attribute__ ((deprecated))
#define DEPRECATED_MSG(msg) __attribute__ ((deprecated(msg)))

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
  e_subJEC,     /**< jet energy scale uncertainty subjets */
  e_JER,     /**< jet resolution uncertainty */
  e_subJER,     /**< jet resolution uncertainty */
  e_MuonSF,  /**< muon scale factor uncertainty */
  e_EleSF,   /**< electron scale factor uncertainty */
  e_TauSF,   /**< tau scale factor uncertainty */
  e_TauEleSF, /**< e-> tau fake rate scale factor uncertainty */
  e_TauEffSF, /**< tau identification scale factor uncertainty */
  e_TauEnergy,/**< tau energy scale uncertainty */
  e_PDF,     /**< PDF uncertainty */
  e_TER     /**< tau energy resolution */

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

/**
 * list of possible jet types, needed to determine the jet radius in jet shape function
 * @see EventCalc::IntegratedJetShape
*/
enum E_JetType {
    e_CA8, 
    e_CA15, 
    e_AK5
};

/** \brief Format information in a table.
 * 
 */
class TableOutput{
public:
    
    /// construct, setting the column headings. Note that this defines now and for always the number of columns.
    explicit TableOutput(const std::vector<std::string> & header);
    
    /// Append a row to the end of the table; the number of columns must match the one defined at the time of construction
    void add_row(const std::vector<std::string> & row);
    
    /// print the table to the supplied stream.
    void print(std::ostream & out);
    
private:
    
    size_t ncols;
    std::vector<std::string> header;
    std::vector<std::vector<std::string> > rows;
};


// make a logarithmic binning. The resulting array conatins n_bins+1 entries with the lower and upper
// end of n_bins bins, suitable for passing it to the TH1 constructor.
boost::shared_array<double> log_binning(size_t n_bins, double xmin, double xmax);

//tmva tagger variables
Double_t HelicityAngle(TLorentzVector j1, TLorentzVector j2);

float HiggsBRweight();

int subJetBTag(TopJet topjet, E_BtagType type, TString mode="default",TString filename="");

int subJetBTagOne(TopJet topjet, E_BtagType type, TString mode="default", TString filename="", int whichsub=-1);

int subJetBTagTop(TopJet topjet, E_BtagType type, TString mode="default",TString filename="");

bool HiggsTag(TopJet topjet, E_BtagType type1, E_BtagType type2, TString mode="default",TString filename="");

bool HepTopTagMatch(TopJet topjet);

float HepTopTagMatchMass(TopJet topjet);

float HepTopTagMatchPt(TopJet topjet);

/* int subJetBTag(TopJet topjet, E_BtagType type); */

/* bool HiggsTag(TopJet topjet, E_BtagType type1, E_BtagType type2); */

bool HepTopTagFull(TopJet topjet, std::vector<PFParticle>* allparts);

/// return if a jet is tagged given a btag operating point
bool IsTagged(Jet &, E_BtagType);


/**
 * basic top tagging routine
 *
 * Test, if the jet topjet fulfills the top-tagging criteria. The tagging variables mjet, nsubjets and mmin are determined and returned.
 * @see NTopTagSelection
 */

bool variableHepTopTag(TopJet topjet, double ptJetMin = 200., double massWindowLower = 0.85, double massWindowUpper = 1.15, double cutCondition2 = 0.35, double cutCondition3 = 0.35);
bool variableHepTopTagWithMatch(TopJet topjet, double ptJetMin = 200., double massWindowLower = 0.85, double massWindowUpper = 1.15, double cutCondition2 = 0.35, double cutCondition3 = 0.35);

double HiggsMassFromSubjets(TopJet topjet);
double HiggsMassFromBTaggedSubjets(TopJet topjet, E_BtagType type, TString mode="default", TString filename="");

double WMassWithMatch(TopJet topjet);

double HepTopTagPairwiseMassWithMatch1(TopJet topjet, double ptJetMin = 200., double massWindowLower = 0.85, double massWindowUpper = 1.15, double cutCondition2 = 0.35, double cutCondition3 = 0.35);

double HepTopTagPairwiseMassWithMatch2(TopJet topjet, double ptJetMin = 200., double massWindowLower = 0.85, double massWindowUpper = 1.15, double cutCondition2 = 0.35, double cutCondition3 = 0.35);


bool HepTopTag(TopJet topjet);
bool HepTopTagWithMatch(TopJet topjet);
bool HepTopTagInverted(TopJet topjet);
bool TopTag(TopJet topjet, double& mjet, int& nsubjets, double& mmin);
bool variableTopTag(TopJet topjet, double &mjet, int &nsubjets, double &mmin, double mminLower = 50., double mjetLower = 140., double mjetUpper = 250.);

/**
 * W tagging routine
 *
 * @see NWTagSelection
 */
bool WTag(TopJet prunedjet, double& mjet, int &nsubjets, double& massdrop);

// return flavor of subjet, mimic algorithmic definition meant only for some tests
int SubJetFlavor(TopJet topjet, int whichsub, float matchcone=0.3);

float relIsoMuon(const Muon & mu, float deltaR = 0.4); //DEPRECATED_MSG("use the version with EventCalc argument");

/// re-calculate the PF relative isolation of a muon in a cone radius deltaR, using the PFCandidates in event
float relIsoMuon(EventCalc & event, const Muon & mu, float deltaR = 0.4);

//float relIso(const Particle & particle, float deltaR = 0.4); //DEPRECATED_MSG("use the version with EventCalc argument");

/// re-calculate the PF relative isolation of a muon in a cone radius deltaR, using the PFCandidates in event
float relIso(EventCalc* event, const Particle & particle, float deltaR = 0.4);

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

/// distance in eta-phi space
double deltaR(LorentzVector v1, LorentzVector v2);

/// absolute value of Delta_phi, xi should be in [-pi,pi]
double deltaPhiAbs(double x1, double x2);

double double_infinity();
int int_infinity();

/// x^p
int myPow(int x, unsigned int p) ;

//overload
double pTrel(const Particle & p1, const Particle & p2);
double pTrel(const LorentzVector & p1,const LorentzVector & p2);

#endif
