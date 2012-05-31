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

enum E_BtagType{e_CSVT,e_CSVM,e_CSVL,e_JPT,e_JPM,e_JPL};

bool TopTag(TopJet topjet, double& mjet, int& nsubjets, double& mmin);

double HTlep(const BaseCycleContainer *bcc);
Jet* nextJet(const Particle *p, std::vector<Jet> *jets);
double pTrel(const Particle *p, std::vector<Jet> *jets);
double deltaRmin(const Particle *p, std::vector<Jet> *jets);

double double_infinity();
int int_infinity();

#endif
