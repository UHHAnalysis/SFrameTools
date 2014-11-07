#ifndef SubJetTagger_H
#define SubJetTagger_H

#include "SFrameTools/include/Objects.h"
#include "SFrameTools/include/Utils.h"
#include <map>
#include <string>

class SubJetTagger {
 public:
  virtual ~SubJetTagger();
  virtual bool Tag(const TopJet& topjet) = 0;
  virtual std::map<std::string, double> TagVar()=0;

  //discriminator value is supposed to be treated as part of the TagVariables

};


/** \brief CMS top tagger, with optional tau3/tau2 cut
 * 
 * The default values of the constructor correspond to the default values used in CMS, without tau3/tau2 cut.
 */
class CMSTopTagger : public SubJetTagger {
public:
  
  CMSTopTagger(double mminLower=50., double mjetLower=140., double mjetUpper=250., double tau32Upper = double_infinity());
  ~CMSTopTagger();

  bool Tag (const TopJet& topjet);
  void SetTau32Cut(double cut=0.7);
  // available variables: "nsubjets", "mmin" and "mjet"
  std::map<std::string,double> TagVar();


 private:

  double m_mjet;
  double m_nsubjets;
  double m_mmin;
  double m_mminLower;
  double m_mjetLower;
  double m_mjetUpper;
  double m_tau32Upper;

};

#endif // _SubJetTagger_H




