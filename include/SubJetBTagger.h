#ifndef _SubJetBTagger_H
#define _SubJetBTagger_H

#include "include/SubJetTagger.h"
#include "include/Utils.h"
#include "TH1F.h"
#include "TFile.h"

using namespace std;

class SubJetBTagger : public SubJetTagger
{
 public:
  SubJetBTagger(const E_BtagType& i_type, const TString& i_mode= "default", const TString& i_filename="", const int& i_whichsub = -1);
  ~SubJetBTagger(){};
  bool Tag(const TopJet& topjet);
  map<string, double> TagVar();

 private:


  int cvsCalculator(const int& flav,  const E_BtagType&  type, const float& subeta);
  void fillHistos(const int& flav,  const E_BtagType&  type, const float& subeta);

  double discriminator_cut;

  int nBTagsSub;

  E_BtagType  type;
  TString mode; 
  TString filename;
  bool dosf;

  int whichsub;

  TF1 *csv;
  TF1 *csvu;
  TF1 *csvd;

  TH1F *numpt;
  TH1F *denpt;
  TH1F *effipt;
  TH1F *errbc;


  TFile *file_mc;

};



#endif // _SubJetBTagger_H
