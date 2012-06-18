#ifndef Utils_H
#define Utils_H

#include "Objects.h"
#include "BaseCycleContainer.h"
#include "TVector3.h"
#include <limits>

struct HigherPt {
  bool operator() (const Particle& j1, const Particle& j2) const {
    return j1.pt() > j2.pt();
  };
};

/**
 * list of supported b-tagging working points
 * @see  NBTagSelection
*/
enum E_BtagType{
  e_CSVT, /**< Combined Secondary Vertex tagger, tight working point */
  e_CSVM, /**< Combined Secondary Vertex tagger, medium working point */
  e_CSVL, /**< Combined Secondary Vertex tagger, loose working point */
  e_JPT,  /**< Jet Probability tagger, tight working point */
  e_JPM, /**< Jet Probability tagger, medium working point */
  e_JPL /**< Jet Probability tagger, loose working point */
};

/**
 * basic top tagging routine
 * 
 * Test, if the jet topjet fulfills the top-tagging criteria. The tagging variables mjet, nsubjets and mmin are determined and returned.
 * @see NTopTagSelection
 */
bool TopTag(TopJet topjet, double& mjet, int& nsubjets, double& mmin);

//double HTlep(const BaseCycleContainer *bcc); ->moved to EventCalc

/// returns the next jet in the list jets to the particle p in the eta-phi plane
Jet* nextJet(const Particle *p, std::vector<Jet> *jets);
/// relative momentum of the particle p with respect to the next jet in the list of jets 
double pTrel(const Particle *p, std::vector<Jet> *jets);
/// minimal distance in the eta-phi plane of particle p with respect to the closest jet from the list of jets
double deltaRmin(const Particle *p, std::vector<Jet> *jets);

double double_infinity();
int int_infinity();

#endif
