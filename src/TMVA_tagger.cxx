#include "include/TMVA_tagger.h"
#include "NtupleWriter/include/JetProps.h"

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include "SFrameTools/include/EventCalc.h"

using namespace std;

TMVA_tagger* TMVA_tagger::m_instance = NULL;

TMVA_tagger* TMVA_tagger::Instance()
{
  // Get a pointer to the object handler.
  // This is the only way to access this class, 
  // since it's a singleton. This method is accessible
  // from everywhere.

  if (m_instance == NULL){
    m_instance = new TMVA_tagger();
  }
  return m_instance; 
   
}


TMVA_tagger::TMVA_tagger() : m_logger("TMVA_tagger")
{
  // constructor: initialise all variables

  
  m_logger << DEBUG << "Constructor called." << SLogger::endmsg;
  
     //Reset();

  
}



TMVA_tagger::~TMVA_tagger()
{
  // default destructor
  
}

bool TMVA_tagger::IsTobiasTagged(TopJet topjet){
  TMVA_tagger::push_variables(topjet);
  //cuts
  // if(Nsubjets>1 && Nsubjets<4)
    //if(subjet1pt<300)
  // if(tau1>0.36 && tau1<0.64)
  //	if(t3t2>0 && t3t2<0.72)
	//  if(tau1pruned>0.32 && tau1pruned<0.64)
	 if(Subjet12mass>80 && Subjet12mass<110)
	{
	    TMVA_tagger::reset_variables();
	    return true;
	  }
      //   else return false;
      // else return false;
      //	else return false;
      //  else return false;
  //else return false;
  else return false;
  
}

bool TMVA_tagger::IsTagged(TString taggername, TopJet topjet, double cutvalue, double &mva_value){
  // double mva_value;
  if(reader==NULL) TMVA_tagger::Set_Reader(taggername);
  TMVA_tagger::push_variables(topjet);
  mva_value=TMVA_tagger::GetMVA_value();
  TMVA_tagger::reset_variables();
  return mva_value>cutvalue;
}

void TMVA_tagger::reset_variables(){
 Nsubjets=0; //done
 Topjetmass=0; //done
 Subjet12mass=0; //done
 Subjet13mass=0; //done
 Subjet23mass=0; //done
 t3t2=0; //done
 t2t1=0; //done
 Q_valatility=0; //done
 npv=0; 
 Subjet12mass_pruned=0; //done
 Subjet13mass_pruned=0; //done
 Subjet23mass_pruned=0; //done
 t3t2p=0;
 t2t1p=0;
 subjet1pt=0; //sjpt1 done
 subjet2pt=0;//done	 
 subjet3pt=0; //done
 tau1=0; //done
 tau2=0; //done
 tau3=0; //done
 tau4=0;
 Topjetmass_pruned=0; //done
 tau1pruned=0; //done t1p
 tau2pruned=0; //done t2p
 tau3pruned=0; //done t3p
 HelAng12=0; //done HA12
 HelAng13=0; //done HA13
 HelAng23=0; //done HA23
 P02=0;
 P04=0;
 P06=0;
 P08=0;
 P10=0;
 w=0;
 st1=0;
 st2=0;
 st3=0;

}


void TMVA_tagger::Set_Reader(TString option) {
  reader=new TMVA::Reader();

if(option=="TLflat"){
    reader->AddVariable( "TopJet_mass", &Topjetmass );
    reader->AddVariable( "Subjets12_mass", &Subjet12mass);
    reader->AddVariable( "Subjets23_mass", &Subjet23mass);
    reader->AddVariable( "TopJet_Nsubjets", &Nsubjets);
    reader->AddVariable( "TopJet_Qjets_volatility", &Q_valatility);
    reader->AddVariable( "Subjet1pT := sqrt(pow(Subjet1_px,2)+pow(Subjet1_py,2))", &subjet1pt);
    reader->AddVariable( "TopJet_pruned_tau3", &tau3pruned);
    reader->AddVariable( "t3/t2_pruned := TopJet_pruned_tau3/TopJet_pruned_tau2", &t3t2p);
    reader->AddVariable( "t2/t1_pruned := TopJet_pruned_tau2/TopJet_pruned_tau1", &t2t1p );
    reader->AddVariable( "HelAng13", &HelAng13);
      reader->BookMVA("BDTG", "/nfs/dust/cms/user/tlapsien/TMVA/TopTag/weights/flatTL_BDTG.weights.xml");
    // reader->BookMVA("BDTG", "/nfs/dust/cms/user/tlapsien/TMVA/TopTag/weights/highptTL_BDTG.weights.xml");

  }



  if(option=="NPVweight"){
    reader->AddVariable( "TopJet_mass", &Topjetmass );
    reader->AddVariable( "Subjets12_mass", &Subjet12mass);
    reader->AddVariable( "Subjets23_mass", &Subjet23mass);
    reader->AddVariable( "TopJet_Nsubjets", &Nsubjets);
    reader->AddVariable( "TopJet_Qjets_volatility", &Q_valatility);
    reader->AddVariable( "Subjet1pT := sqrt(pow(Subjet1_px,2)+pow(Subjet1_py,2))", &subjet1pt);
    reader->AddVariable( "TopJet_pruned_tau3", &tau3pruned);
    reader->AddVariable( "t3/t2_pruned := TopJet_pruned_tau3/TopJet_pruned_tau2", &t3t2p);
    reader->AddVariable( "t2/t1_pruned := TopJet_pruned_tau2/TopJet_pruned_tau1", &t2t1p );
    reader->AddVariable( "HelAng13", &HelAng13);

    //  reader->BookMVA("BDTG", "/nfs/dust/cms/user/tlapsien/TMVA/weights/NPVweight_BDTG.weights.xml");
    reader->BookMVA("BDTG", "/nfs/dust/cms/user/tlapsien/TMVA/TopTag/weights/secondTL_BDTG.weights.xml");
 

  }
  
 if(option=="Uncorr+3"){
    reader->AddVariable( "TopJet_mass", &Topjetmass );
    reader->AddVariable( "Subjets12_mass", &Subjet12mass);
    reader->AddVariable( "Subjets23_mass", &Subjet23mass);
    reader->AddVariable( "TopJet_Nsubjets", &Nsubjets);
    reader->AddVariable( "TopJet_Qjets_volatility", &Q_valatility);
    reader->AddVariable( "Event_NPV", &npv);
    reader->AddVariable( "Subjet1pT := sqrt(pow(Subjet1_px,2)+pow(Subjet1_py,2))", &subjet1pt);
    reader->AddVariable( "TopJet_pruned_tau3", &tau3pruned);
    reader->AddVariable( "t3/t2_pruned := TopJet_pruned_tau3/TopJet_pruned_tau2", &t3t2p);
    reader->AddVariable( "t2/t1_pruned := TopJet_pruned_tau2/TopJet_pruned_tau1", &t2t1p );
    reader->AddVariable( "HelAng13", &HelAng13);

    reader->BookMVA("BDTG", "/nfs/dust/cms/user/tlapsien/TMVA/weights/Uncorr+3_BDTG.weights.xml");
 }

  if(option=="pT"){
    reader->AddVariable( "TopJet_mass", &Topjetmass );
    reader->AddVariable( "Subjets12_mass", &Subjet12mass);
    reader->AddVariable( "Subjets23_mass", &Subjet23mass);
    reader->AddVariable( "TopJet_Nsubjets", &Nsubjets);
    reader->AddVariable( "TopJet_Qjets_volatility", &Q_valatility);
    reader->AddVariable( "Event_NPV", &npv);
    reader->AddVariable( "TopJet_pruned_tau3", &tau3pruned);
    reader->AddVariable( "t3/t2_pruned := TopJet_pruned_tau3/TopJet_pruned_tau2", &t3t2p);
    reader->AddVariable( "t2/t1_pruned := TopJet_pruned_tau2/TopJet_pruned_tau1", &t2t1p );
    reader->AddVariable( "HelAng13", &HelAng13);

    reader->BookMVA("BDTG", "/nfs/dust/cms/user/tlapsien/TMVA/weights/pT_BDTG.weights.xml");
 

  }

 
  if(option=="bestown70"){
     reader->AddVariable( "Subjets12_mass", &Subjet12mass);
    reader->AddVariable( "t2/t1 := TopJet_tau2/TopJet_tau1", &t2t1);
    reader->AddVariable( "Subjets23_pruned_mass", &Subjet23mass_pruned);
    reader->AddVariable( "TopJet_pruned_tau1", &tau1pruned);
    reader->AddVariable( "HelAng12", &HelAng12);


    reader->BookMVA("BDTG", "/nfs/dust/cms/user/tlapsien/TMVA/weights/bestown70_BDTG.weights.xml");
 
    
  }


  if(option=="Q"){
    reader->AddVariable( "TopJet_mass", &Topjetmass );
    reader->AddVariable( "TopJet_Nsubjets", &Nsubjets);
    reader->AddVariable( "Subjets12_mass", &Subjet12mass);
    reader->AddVariable( "Subjets13_mass", &Subjet13mass);
    reader->AddVariable( "Subjets23_mass", &Subjet23mass);
    reader->AddVariable("t3/t2 := TopJet_tau3/TopJet_tau2", &t3t2);
    reader->AddVariable( "TopJet_Qjets_volatility", &Q_valatility);
    
    reader->BookMVA("BDTG", "/nfs/dust/cms/user/tlapsien/TMVA/weights/Qjets_BDTG.weights.xml");
  }
if(option=="Qref_weight"){
     reader->AddVariable( "TopJet_mass", &Topjetmass );
    reader->AddVariable( "Subjets12_mass", &Subjet12mass);
    reader->AddVariable( "Subjets13_mass", &Subjet23mass);
    reader->AddVariable( "Subjets23_mass", &Subjet23mass);
    reader->AddVariable( "TopJet_Nsubjets", &Nsubjets);
    reader->AddVariable( "TopJet_Qjets_volatility", &Q_valatility);
    reader->AddVariable("t3/t2 := TopJet_tau3/TopJet_tau2", &t3t2);
    reader->BookMVA("BDTG", "/nfs/dust/cms/user/tlapsien/TMVA/weights/Qref_BDTG.weights.xml");
   }

if(option=="CMS"){
    reader->AddVariable( "TopJet_mass", &Topjetmass );
    reader->AddVariable( "Subjets12_mass", &Subjet12mass);
    reader->AddVariable( "Subjets13_mass", &Subjet13mass);
    reader->AddVariable( "Subjets23_mass", &Subjet23mass);
    reader->AddVariable( "TopJet_Nsubjets", &Nsubjets);
    
    reader->BookMVA("BDTG", "/nfs/dust/cms/user/tlapsien/TMVA/weights/cmsref_weight_BDTG.weights.xml");
 }

// reader2=reader;
 }

Double_t TMVA_tagger::GetMVA_value(){
  // std::cout<<Tmass<<std::endl;
  
  Double_t mvaValue = reader->EvaluateMVA("BDTG");
  return mvaValue;

}

void TMVA_tagger::push_variables(TopJet topjet){
  // Tmass=(topjet.v4()).M();
  EventCalc* calc = EventCalc::Instance();
  npv = calc->GetPrimaryVertices()->size();
  Nsubjets=topjet.numberOfDaughters();
  LorentzVector allsubjets(0,0,0,0);
  
  for(int j=0; j<topjet.numberOfDaughters(); ++j){
    allsubjets += topjet.subjets()[j].v4();
  }

  if(!allsubjets.isTimelike()){
    Topjetmass=0;
  } else {
    Topjetmass = allsubjets.M();
  }
  
  //subjet mass
  std::vector<Particle> subjet = topjet.subjets();
  sort(subjet.begin(),subjet.end(),HigherPt());
  if(Nsubjets>=3){
    Subjet12mass=(subjet[0].v4()+subjet[1].v4()).mass();
    Subjet13mass=(subjet[0].v4()+subjet[2].v4()).mass();
    Subjet23mass=(subjet[1].v4()+subjet[2].v4()).mass();
  }
   if(Nsubjets>=2){
    Subjet12mass=(subjet[0].v4()+subjet[1].v4()).mass();
     }


   //nsubjetiness
   JetProps jp(&topjet, calc->GetPFParticles() );
   tau1 = jp.GetNsubjettiness(1, Njettiness::onepass_kt_axes, 1., 0.8);
  tau2 = jp.GetNsubjettiness(2, Njettiness::onepass_kt_axes, 1., 0.8);
  tau3 = jp.GetNsubjettiness(3, Njettiness::onepass_kt_axes, 1., 0.8);
  
  //pruned subjetiness
  tau1pruned = jp.GetPrunedNsubjettiness(1, Njettiness::onepass_kt_axes, 1., 0.8);
  tau2pruned = jp.GetPrunedNsubjettiness(2, Njettiness::onepass_kt_axes, 1., 0.8);
  tau3pruned = jp.GetPrunedNsubjettiness(3, Njettiness::onepass_kt_axes, 1., 0.8);
  //ratio subjetiness
  t3t2=tau3/tau2;
  t2t1=tau2/tau1;
  t3t2p=tau3pruned/tau2pruned;
  t2t1p=tau2pruned/tau1pruned;
  Q_valatility = topjet.qjets_volatility();
  

 // calculate pruned masses
  std::vector<fastjet::PseudoJet> jets = jp.GetFastJet(2.0);   // something large to make sure jet is inside radius
  if(jets.empty()){
      m_logger << WARNING << "TMVATreeFiller::FillTopJetProperties: no jet found!" << SLogger::endmsg; 
  }
  else{
    fastjet::PseudoJet pjet = jp.GetPrunedJet(jets[0]);
    std::vector<fastjet::PseudoJet> prunedsubjets;
    if (pjet.constituents().size()>=2){
        prunedsubjets = pjet.exclusive_subjets(2);
    }
    if (pjet.constituents().size()>=3){
        prunedsubjets = pjet.exclusive_subjets(3);
    }

    unsigned int pnsubs = prunedsubjets.size();

    fastjet::PseudoJet psubjets(0,0,0,0);
    
    for(unsigned int j=0; j<pnsubs; ++j){
        psubjets += pjet.pieces()[j];
    }

    Topjetmass_pruned = psubjets.m();

    if (pnsubs>=2) {
        Subjet12mass_pruned = (prunedsubjets[0]+prunedsubjets[1]).m();
    }

    if (pnsubs>=3) {
        Subjet13mass_pruned = (prunedsubjets[0]+prunedsubjets[2]).m();
        Subjet23mass_pruned = (prunedsubjets[1]+prunedsubjets[2]).m();
    }
  }
  
  //HelicityAngles && subjet pt
  if(Nsubjets>0) subjet1pt=sqrt(pow(subjet[0].v4().Px(),2)+pow(subjet[0].v4().Py(),2));
  if(Nsubjets>1){
    TLorentzVector subjet1(subjet[0].v4().Px(),subjet[0].v4().Py(),subjet[0].v4().Pz(),subjet[0].v4().E());
    TLorentzVector subjet2(subjet[1].v4().Px(),subjet[1].v4().Py(),subjet[1].v4().Pz(),subjet[1].v4().E());
    subjet1pt=sqrt(pow(subjet[0].v4().Px(),2)+pow(subjet[0].v4().Py(),2));
    subjet2pt=sqrt(pow(subjet[1].v4().Px(),2)+pow(subjet[1].v4().Py(),2));
    HelAng12= HelicityAngle(subjet1,subjet2);
  }
  if(Nsubjets>2){
    TLorentzVector subjet1(subjet[0].v4().Px(),subjet[0].v4().Py(),subjet[0].v4().Pz(),subjet[0].v4().E());
    TLorentzVector subjet2(subjet[1].v4().Px(),subjet[1].v4().Py(),subjet[1].v4().Pz(),subjet[1].v4().E());
    TLorentzVector subjet3(subjet[2].v4().Px(),subjet[2].v4().Py(),subjet[2].v4().Pz(),subjet[2].v4().E());
    subjet1pt=sqrt(pow(subjet[0].v4().Px(),2)+pow(subjet[0].v4().Py(),2));
    subjet2pt=sqrt(pow(subjet[1].v4().Px(),2)+pow(subjet[1].v4().Py(),2));
    subjet3pt=sqrt(pow(subjet[2].v4().Px(),2)+pow(subjet[2].v4().Py(),2));
    HelAng12= HelicityAngle(subjet1,subjet2);
    HelAng13= HelicityAngle(subjet1,subjet3);
    HelAng23= HelicityAngle(subjet2,subjet3);
  }



}

