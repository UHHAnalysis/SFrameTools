#ifndef _SubJetTagger_H
#define _SubJetTagger_H

#include "SFrameTools/include/Objects.h"
#include <map> 

using namespace std;

class SubJetTagger
{
 public:
  SubJetTagger(){};
  virtual ~SubJetTagger(){};
  virtual bool Tag(const TopJet& topjet) = 0;
  virtual map<string, double> TagVar()=0;

  //discriminator value is supposed to be treated as part of the TagVariables

};



class CMSTopTagger : public SubJetTagger
{
 public:
  CMSTopTagger(const double& mminLower=50., const double& mjetLower=140., const double& mjetUpper=250.);
  ~CMSTopTagger(){};

  bool Tag (const TopJet& topjet);
  map<string,double> TagVar();


 private:

  double m_mjet;
  double m_nsubjets;
  double m_mmin;
  double m_mminLower;
  double m_mjetLower; 
  double m_mjetUpper;

};

#endif // _SubJetTagger_H




