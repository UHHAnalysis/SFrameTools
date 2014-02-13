// Dear emacs, this is -*- c++ -*-
// $Id: SFrameTools_LinkDef.h,v 1.7 2013/06/19 16:15:46 peiffer Exp $
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;

// Add the declarations of your cycles, and any other classes for which you
// want to generate a dictionary, here. The usual format is:
//
// #pragma link C++ class MySuperClass+;

#pragma link C++ class EventCalc+;
#pragma link C++ class FJet+;
#pragma link C++ class BaseHists+;
#pragma link C++ class LuminosityHandler+;
#pragma link C++ struct BaseCycleContainer+;
//#pragma link C++ struct Selection+;
#pragma link C++ class AnalysisModuleRunner+;
#pragma link C++ struct PUWeightProducer+;
#pragma link C++ class BaseHists+;
#pragma link C++ class LeptonScaleFactors+;
#pragma link C++ class TriggerWeight+;
#pragma link C++ class ReconstructionHypothesis+;
#pragma link C++ class TopFitCalc+;
#pragma link C++ class std::vector<ReconstructionHypothesis>+;
#pragma link C++ struct qualityflag+;

#pragma link C++ class LuminosityUtils::RunNr_LbNr+;
#pragma link C++ class LuminosityUtils::LumiBinInfo+;
#pragma link C++ class pair<LuminosityUtils::RunNr_LbNr,double>+;
#pragma link C++ class pair<LuminosityUtils::RunNr_LbNr,int>+;
#pragma link C++ class pair<int,LuminosityUtils::LumiBinInfo*>+;
#pragma link C++ class pair<double,unsigned int>+;
#pragma link C++ class ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<Double32_t> >+;
#pragma link C++ class ROOT::Math::PtEtaPhiE4D<Double32_t>+;

// needed by PROOF (don't really know why)
#pragma link C++ class pair<const string,string*>;
#pragma link C++ class pair<const string,int*>;
#pragma link C++ class pair<const string,double*>;
#pragma link C++ class pair<const string,bool*>;
#pragma link C++ class pair<const string,vector<string>*>;
#pragma link C++ class pair<const string,vector<int>*> ;
#pragma link C++ class pair<const string,vector<double>*>;
#pragma link C++ class pair<const string,vector<bool>*>;
#pragma link C++ class pair<pair<string,string>,TH1*>;
#pragma link C++ class TObjLink;



#endif // __CINT__
