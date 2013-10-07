#include "include/Utils.h"
#include "NtupleWriter/include/JetProps.h"

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
namespace external {
#include "include/HEPTopTagger.h"
}
#include "SFrameTools/include/EventCalc.h"

using namespace std;


bool HepTopTagMatch(TopJet topjet){

   EventCalc* calc = EventCalc::Instance();

   BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

   double deltarmin = double_infinity();

   TopJet nextjet;

   for(unsigned int it=0; it<bcc->toptagjets->size();++it){

     TopJet top4jet=bcc->toptagjets->at(it);

     if(top4jet.deltaR(topjet) < deltarmin){
       deltarmin = top4jet.deltaR(topjet);
       nextjet = top4jet;
     }

  }

   if(deltarmin<0.3){

     return HepTopTag(nextjet);

   }
   else return 0;

}

float HepTopTagMatchMass(TopJet topjet){

   EventCalc* calc = EventCalc::Instance();

   BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

   double deltarmin = double_infinity();

   TopJet nextjet;

   for(unsigned int it=0; it<bcc->toptagjets->size();++it){

     TopJet top4jet=bcc->toptagjets->at(it);

     if(top4jet.deltaR(topjet) < deltarmin){
       deltarmin = top4jet.deltaR(topjet);
       nextjet = top4jet;
     }

  }

   if(deltarmin>=0.3){

     return -99999.;

   }
   else return nextjet.v4().M();

}

std::unique_ptr<double[]> log_binning(size_t n_bins, double xmin, double xmax){
    assert(xmin > 0 && xmin < xmax);
    std::unique_ptr<double[]> result(new double[n_bins + 1]);
    
    // equidistant binning og log-scale means always use the same ratio between bin borders:
    double ratio = pow(xmax / xmin, 1.0 / n_bins);
    result[0] = xmin;
    for(size_t i=1;i <= n_bins; ++i){
        result[i] = ratio * result[i-1];
    }
    return result;
}


float HiggsBRweight(){

  float addweight=1.;

  EventCalc* calc = EventCalc::Instance();
 
  // ObjectHandler* objs = ObjectHandler::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  TFile *fileweight = new
TFile("/scratch/hh/dust/naf/cms/user/imarches/newSFrame/SFrame/SFrameAnalysis/config/correctHiggsBR.root",
"READ");

  TH1F *Higgs_BR=(TH1F*)fileweight->Get("Higgs_BR");

  for(unsigned int ig=0; ig<bcc->genparticles->size(); ++ig){
   
    GenParticle genp = bcc->genparticles->at(ig);
   
    if (abs(genp.pdgId()) == 25){
     
      const GenParticle *dau1;
      const GenParticle *dau2;
     
      dau1=genp.daughter(bcc->genparticles,1);
      dau2=genp.daughter(bcc->genparticles,2);
     
      int flav1=abs(dau1->pdgId());
      int flav2=abs(dau2->pdgId());
     
      if (dau1&&dau2){
   
    if(flav1==3&&flav2==3){
      addweight=addweight*Higgs_BR->GetBinContent(1);
    }
    else if(flav1==4&&flav2==4){
      addweight=addweight*Higgs_BR->GetBinContent(2);
    }
    else if(flav1==5&&flav2==5){
      addweight=addweight*Higgs_BR->GetBinContent(3);
    }
    else if(flav1==13&&flav2==13){
      addweight=addweight*Higgs_BR->GetBinContent(4);
    }
    else if(flav1==15&&flav2==15){
      addweight=addweight*Higgs_BR->GetBinContent(5);
    }
    else if(flav1==21&&flav2==21){
      addweight=addweight*Higgs_BR->GetBinContent(6);
    }
    else if(flav1==22&&flav2==22){
      addweight=addweight*Higgs_BR->GetBinContent(7);
    }
    else if(flav1==23&&flav2==23){
      addweight=addweight*Higgs_BR->GetBinContent(8);
    }
    else if(flav1==24&&flav2==24){
      addweight=addweight*Higgs_BR->GetBinContent(9);
    }
    else if((flav1==22&&flav2==23)||(flav1==23&&flav2==22)){
      addweight=addweight*Higgs_BR->GetBinContent(10);
    }
    else{
      cout << "Forgotten Higgs decays into " << flav1 << " " << flav2 <<
endl;
    }
   
      }
     
    }//is an Higgs
   
  }//gen particles loop

  fileweight->Close();

  delete fileweight;

  return addweight;

}

#include "TRandom3.h"

//subjet b-tagger, returns number of b-tagged subjets

int subJetBTag(TopJet topjet, E_BtagType type, TString mode, TString filename){

  //Modes:
  //default --> no SF
  //mean --> SF
  //lightup,lightdown,bcup,bcdown --> SF systematics evaluation

  int nBTagsSub=0;

  TString syst=mode;

  float discriminator_cut;

  bool dosf=1;

  if(type==e_CSVL) discriminator_cut = 0.244;
  if(type==e_CSVM) discriminator_cut = 0.679;
  if(type==e_CSVT) discriminator_cut = 0.898;
 
  std::vector<Particle> subjets_top;
  std::vector<float> btagsub_combinedSecondaryVertex_top;
  std::vector<int> flavorsub_top;

  subjets_top=topjet.subjets();
  btagsub_combinedSecondaryVertex_top=topjet.btagsub_combinedSecondaryVertex();
  flavorsub_top=topjet.flavorsub();
 
  if(filename==""&&mode!="default"){
    std::cout << "ATTENTION!!! Asked b-tagging SF, but no efficiencies provided! SF will NOT be applied!!!" << std::endl;
  }

  if(filename==""||mode=="default"){
    dosf=0;
  }
  
  if(mode!="default"&&mode!="mean"&&mode!="lightup"&&mode!="lightdown"&&mode!="bcup"&&mode!="bcdown"){
    std::cout << "ATTENTION!!! B-tagging SF mode "<< mode <<" not known! Will NOT perform any SF re-weighting!" << std::endl;
    dosf=0;
  }

  TFile *file_mc;

  if(dosf){
    file_mc = new TFile(filename);
    
    file_mc->cd();
  }
  
  for(unsigned int i=0; i < btagsub_combinedSecondaryVertex_top.size(); ++i){

    Particle subjet=subjets_top[i];
    int flav = flavorsub_top[i];
    float subpt=subjet.pt();
    float subeta=subjet.eta();
    float test=btagsub_combinedSecondaryVertex_top[i];

    if(!dosf){
      if(test>discriminator_cut){
	nBTagsSub += 1;
      }     
      continue;
    }

    TF1 *csv;
    TF1 *csvu;
    TF1 *csvd;
    
    TH1F *numpt;
    TH1F *denpt;
    TH1F *effipt;
    TH1F *errbc;

    TRandom3* rand;

    double bc_bins[] = {
      20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800
    };
    errbc=new TH1F("errbc","shift bc", 16, bc_bins);

    int possible=1;
    
    if(type==e_CSVL){
      if(abs(flav)==4||abs(flav)==5){
	csv = new TF1("csv", "0.981149*((1.+(-0.000713295*x))/(1.+(-0.000703264*x)))", 20.0, 800.0);
	double SFbc_error[] = {
	  0.0484285,
	  0.0126178,
	  0.0120027,
	  0.0141137,
	  0.0145441,
	  0.0131145,
	  0.0168479,
	  0.0160836,
	  0.0126209,
	  0.0136017,
	  0.019182,
	  0.0198805,
	  0.0386531,
	  0.0392831,
	  0.0481008,
	  0.0474291 };
	for(Int_t ibc=1;ibc<17;ibc++){
	  errbc->SetBinContent(ibc,SFbc_error[ibc-1]);
	}
      }
      else if(abs(flav)==1||abs(flav)==2||abs(flav)==3||abs(flav)==21){
	if(fabs(subeta)>=0&&fabs(subeta)<=0.5)
	  {
	    csv = new TF1("csv", "((1.04901+(0.00152181*x))+(-3.43568e-06*(x*x)))+(2.17219e-09*(x*(x*x)))", 20.0, 800.0);
	    csvu = new TF1("csvu", "((1.12424+(0.00201136*x))+(-4.64021e-06*(x*x)))+(2.97219e-09*(x*(x*x)))", 20.0, 800.0);
	    csvd = new TF1("csvd", "((0.973773+(0.00103049*x))+(-2.2277e-06*(x*x)))+(1.37208e-09*(x*(x*x)))", 20.0, 800.0);
	  }
	else if(fabs(subeta)>0.5&&fabs(subeta)<=1.0)
	  {
	    csv = new TF1("csv", "((0.991915+(0.00172552*x))+(-3.92652e-06*(x*x)))+(2.56816e-09*(x*(x*x)))", 20.0, 800.0);
	    csvu = new TF1("csvu", "((1.06231+(0.00215815*x))+(-4.9844e-06*(x*x)))+(3.27623e-09*(x*(x*x)))", 20.0, 800.0);
	    csvd = new TF1("csvd", "((0.921518+(0.00129098*x))+(-2.86488e-06*(x*x)))+(1.86022e-09*(x*(x*x)))", 20.0, 800.0);
	  }
	else if(fabs(subeta)>1.0&&fabs(subeta)<=1.5)
	  {
	    csv = new TF1("csv", "((0.962127+(0.00192796*x))+(-4.53385e-06*(x*x)))+(3.0605e-09*(x*(x*x)))", 20.0, 800.0);
	    csvu = new TF1("csvu", "((1.02883+(0.00231985*x))+(-5.57924e-06*(x*x)))+(3.81235e-09*(x*(x*x)))", 20.0, 800.0);
	    csvd = new TF1("csvd", "((0.895419+(0.00153387*x))+(-3.48409e-06*(x*x)))+(2.30899e-09*(x*(x*x)))", 20.0, 800.0);
	  }
	else if(fabs(subeta)>1.5&&fabs(subeta)<=2.4)
	  {
	    csv = new TF1("csv", "((1.06121+(0.000332747*x))+(-8.81201e-07*(x*x)))+(7.43896e-10*(x*(x*x)))", 20.0, 700.0);
	    csvu = new TF1("csvu", "((1.1388+(0.000468418*x))+(-1.36341e-06*(x*x)))+(1.19256e-09*(x*(x*x)))", 20.0, 700.0);
	    csvd = new TF1("csvd", "((0.983607+(0.000196747*x))+(-3.98327e-07*(x*x)))+(2.95764e-10*(x*(x*x)))", 20.0, 700.0);
	  }
	else{
	  possible=0;
	}
      }
      else{
	possible=0;
      }
    }
    if(type==e_CSVM){
      if(abs(flav)==4||abs(flav)==5){
	csv = new TF1("csv", "0.726981*((1.+(0.253238*x))/(1.+(0.188389*x)))", 20.0, 800.0);
	double SFbc_error[] = {
	  0.0554504,
	  0.0209663,
	  0.0207019,
	  0.0230073,
	  0.0208719,
	  0.0200453,
	  0.0264232,
	  0.0240102,
	  0.0229375,
	  0.0184615,
	  0.0216242,
	  0.0248119,
	  0.0465748,
	  0.0474666,
	  0.0718173,
	  0.0717567 };
	for(Int_t ibc=1;ibc<17;ibc++){
	  errbc->SetBinContent(ibc,SFbc_error[ibc-1]);
	}
      }
      else if(abs(flav)==1||abs(flav)==2||abs(flav)==3||abs(flav)==21){
	if(fabs(subeta)>=0&&fabs(subeta)<=0.8)
	  {
	    csv = new TF1("csv", "((1.06238+(0.00198635*x))+(-4.89082e-06*(x*x)))+(3.29312e-09*(x*(x*x)))", 20.0, 800.0);
	    csvu = new TF1("csvu", "((1.15201+(0.00292575*x))+(-7.41497e-06*(x*x)))+(5.0512e-09*(x*(x*x)))", 20.0, 800.0);
	    csvd = new TF1("csvd", "((0.972746+(0.00104424*x))+(-2.36081e-06*(x*x)))+(1.53438e-09*(x*(x*x)))", 20.0, 800.0);
	  }
	else if(fabs(subeta)>0.8&&fabs(subeta)<=1.6)
	  {
	    csv = new TF1("csv", "((1.08048+(0.00110831*x))+(-2.96189e-06*(x*x)))+(2.16266e-09*(x*(x*x)))", 20.0, 800.0);
	    csvu = new TF1("csvu", "((1.17735+(0.00156533*x))+(-4.32257e-06*(x*x)))+(3.18197e-09*(x*(x*x)))", 20.0, 800.0);
	    csvd = new TF1("csvd", "((0.9836+(0.000649761*x))+(-1.59773e-06*(x*x)))+(1.14324e-09*(x*(x*x)))", 20.0, 800.0);
	  }
	else if(fabs(subeta)>1.6&&fabs(subeta)<=2.4)
	  {
	    csv = new TF1("csv", "((1.09145+(0.000687171*x))+(-2.45054e-06*(x*x)))+(1.7844e-09*(x*(x*x)))", 20.0, 700.0);
	    csvu = new TF1("csvu", "((1.17671+(0.0010147*x))+(-3.66269e-06*(x*x)))+(2.88425e-09*(x*(x*x)))", 20.0, 700.0);
	    csvd = new TF1("csvd", "((1.00616+(0.000358884*x))+(-1.23768e-06*(x*x)))+(6.86678e-10*(x*(x*x)))", 20.0, 700.0);
	  }
	else{
	  possible=0;
	}
      }
      else{
	possible=0;
      }
    }
    if(type==e_CSVT){
      if(abs(flav)==4||abs(flav)==5){
	csv = new TF1("csv", "0.869965*((1.+(0.0335062*x))/(1.+(0.0304598*x)))", 20.0, 800.0);
	double SFbc_error[] = {
	  0.0567059,
	  0.0266907,
	  0.0263491,
	  0.0342831,
	  0.0303327,
	  0.024608,
	  0.0333786,
	  0.0317642,
	  0.031102,
	  0.0295603,
	  0.0474663,
	  0.0503182,
	  0.0580424,
	  0.0575776,
	  0.0769779,
	  0.0898199 };
	for(Int_t ibc=1;ibc<17;ibc++){
	  errbc->SetBinContent(ibc,SFbc_error[ibc-1]);
	}
      }
      else if(abs(flav)==1||abs(flav)==2||abs(flav)==3||abs(flav)==21){
	if(fabs(subeta)>=0&&fabs(subeta)<=2.4)
	  {
	    csv = new TF1("csv", "((1.01739+(0.00283619*x))+(-7.93013e-06*(x*x)))+(5.97491e-09*(x*(x*x)))", 20.0, 700.0);
	    csvu = new TF1("csvu", "((1.08119+(0.00441909*x))+(-1.18764e-05*(x*x)))+(8.71372e-09*(x*(x*x)))", 20.0, 700.0);
	    csvd = new TF1("csvd", "((0.953587+(0.00124872*x))+(-3.97277e-06*(x*x)))+(3.23466e-09*(x*(x*x)))", 20.0, 700.0);
	  }
	else{
	  possible=0;
	}
      }
      else{
	possible=0;
      }
    }
    
    if(!possible){
      if(test>discriminator_cut){
	nBTagsSub += 1;
      }

      delete errbc;
      continue;
    }
    
    double SF=0;
    
    if (csv->GetXmin() > subpt){
      SF=csv->Eval(csv->GetXmin());
    }
    else if (csv->GetXmax() < subpt){ 
      SF=csv->Eval(csv->GetXmax());
    }
    else{
      SF=csv->Eval(subpt);
    } 
    
    double addSF=0;
    
    int bin;
    
    int doubleunc=0;
    
    if((syst=="bcup"||syst=="bcdown")&&(abs(flav)==5||abs(flav)==4)){
      if(subpt>=800){
	bin=errbc->GetXaxis()->GetNbins();
	doubleunc=1;
      }
      else if(subpt<20){
	bin=1;
	doubleunc=1;
      }
      else{
	bin=errbc->GetXaxis()->FindBin(subpt);
      }
      addSF=errbc->GetBinContent(bin);
      if(doubleunc){
	addSF=2*addSF;
      }
      if(abs(flav)==4){
	addSF=2*addSF;
      }
      if(syst=="bcdown"){
	addSF=-addSF;
      }
    }
    
    SF=SF+addSF;
    
    if(syst=="lightup"&&(abs(flav)==1||abs(flav)==2||abs(flav)==3||abs(flav)==21)){
      if (csvu->GetXmin() > subpt){
	SF=csvu->Eval(csvu->GetXmin());
      }
      else if (csvu->GetXmax() < subpt){ 
	SF=csvu->Eval(csvu->GetXmax());
      }
      else{
	SF=csvu->Eval(subpt);
      } 
    }
    
    if(syst=="lightdown"&&(abs(flav)==1||abs(flav)==2||abs(flav)==3||abs(flav)==21)){
      if (csvd->GetXmin() > subpt){
	SF=csvd->Eval(csvd->GetXmin());
      }
      else if (csvd->GetXmax() < subpt){ 
	SF=csvd->Eval(csvd->GetXmax());
      }
      else{
	SF=csvd->Eval(subpt);
      } 
    }
    
    if(type==e_CSVL){
      if(abs(flav)==5){
	double jetpt_bins[] = {
	  20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800, 1600
	};
	denpt=(TH1F*) file_mc->Get("BTagEff/pt_bJet");
	numpt=(TH1F*) file_mc->Get("BTagEff/pt_bJet_bTagL");
	effipt=new TH1F("effipt","effi pt", 17, jetpt_bins);
      }
      if(abs(flav)==4){
	double jetpt_bins[] = {
	  20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 1600
	};
	denpt=(TH1F*) file_mc->Get("BTagEff/pt_cJet");
	numpt=(TH1F*) file_mc->Get("BTagEff/pt_cJet_bTagL");
	effipt=new TH1F("effipt","effi pt", 14, jetpt_bins);
      }
      if(abs(flav)==1||abs(flav)==2||abs(flav)==3||abs(flav)==21){
	double jetpt_bins[] = {
	  20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 1600
	};
	if(fabs(subeta)>=0&&fabs(subeta)<=0.5)
	  {
	    denpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e1L");
	    numpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e1_bTagL");
	  }
	if(fabs(subeta)>0.5&&fabs(subeta)<=1.0)
	  {
	    denpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e2L");
	    numpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e2_bTagL");
	  }
	if(fabs(subeta)>1.0&&fabs(subeta)<=1.5)
	  {
	    denpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e3L");
	    numpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e3_bTagL");
	  }
	if(fabs(subeta)>1.5&&fabs(subeta)<=2.4)
	  {
	    denpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e3L");
	    numpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e3_bTagL");
	  }
	
	effipt=new TH1F("effipt","effi pt", 16, jetpt_bins);
      } 
    }

    if(type==e_CSVM){
      if(abs(flav)==5){
	double jetpt_bins[] = {
	  20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800, 1600
	};
	denpt=(TH1F*) file_mc->Get("BTagEff/pt_bJet");
	numpt=(TH1F*) file_mc->Get("BTagEff/pt_bJet_bTagM");
	effipt=new TH1F("effipt","effi pt", 17, jetpt_bins);
      }
      if(abs(flav)==4){
	double jetpt_bins[] = {
	  20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 1600
	};
	denpt=(TH1F*) file_mc->Get("BTagEff/pt_cJet");
	numpt=(TH1F*) file_mc->Get("BTagEff/pt_cJet_bTagM");
	effipt=new TH1F("effipt","effi pt", 14, jetpt_bins);
      }
      if(abs(flav)==1||abs(flav)==2||abs(flav)==3||abs(flav)==21){
	double jetpt_bins[] = {
	  20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 1600
	};
	if(fabs(subeta)>=0&&fabs(subeta)<=0.8)
	  {
	    denpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e1M");
	    numpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e1_bTagM");
	  }
	if(fabs(subeta)>0.8&&fabs(subeta)<=1.6)
	  {
	    denpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e2M");
	    numpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e2_bTagM");
	  }
	if(fabs(subeta)>1.6&&fabs(subeta)<=2.4)
	  {
	    denpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e3M");
	    numpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e3_bTagM");
	  }
	effipt=new TH1F("effipt","effi pt", 16, jetpt_bins);
      } 
    }

    if(type==e_CSVT){
      if(abs(flav)==5){
	double jetpt_bins[] = {
	  20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800, 1600
	};
	denpt=(TH1F*) file_mc->Get("BTagEff/pt_bJet");
	numpt=(TH1F*) file_mc->Get("BTagEff/pt_bJet_bTagT");
	effipt=new TH1F("effipt","effi pt", 17, jetpt_bins);
      }
      if(abs(flav)==4){
	double jetpt_bins[] = {
	  20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 1600
	};
	denpt=(TH1F*) file_mc->Get("BTagEff/pt_cJet");
	numpt=(TH1F*) file_mc->Get("BTagEff/pt_cJet_bTagT");
	effipt=new TH1F("effipt","effi pt", 14, jetpt_bins);
      }
      if(abs(flav)==1||abs(flav)==2||abs(flav)==3||abs(flav)==21){
	double jetpt_bins[] = {
	  20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 1600
	};
	if(fabs(subeta)>=0&&fabs(subeta)<=2.4)
	  {
	    denpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e1T");
	    numpt=(TH1F*) file_mc->Get("BTagEff/pt_lJet_e1_bTagT");
	  }
	effipt=new TH1F("effipt","effi pt", 16, jetpt_bins);
      } 
    }

    effipt->Divide(numpt,denpt,1,1,"B");
    
    float effi;
    
    if(subpt>=1600){
      bin=effipt->GetXaxis()->GetNbins();
    }
    else if(subpt<20){
      bin=1;
    }
    else{
      bin=effipt->GetXaxis()->FindBin(subpt);
    }
    
    effi=effipt->GetBinContent(bin);
    
    double etaseed = subeta;
    double sin_eta = sin(etaseed*1000000);
    int seed = abs(static_cast<int>(sin_eta*100000));
    
    rand = new TRandom3(seed);
    
    float coin = rand->Uniform(1.);
    
    bool istagged=false;
    
    if(test>discriminator_cut) istagged=true;
    
    bool newtag = istagged;
    
    if(SF > 1){
      
      if( !istagged ) {
	
	//fraction of jets that needs to be upgraded
	float mistagPercent = (1.0 - SF) / (1.0 - (SF/effi) );
	
	//upgrade to tagged
	if( coin < mistagPercent ) {newtag = true;}
      }
      
    }else{
      
      //downgrade tagged to untagged
      if( istagged && coin > SF ) {newtag = false;}
      
    }
    
    //if(istagged!=newtag)
      // cout << "Flavor " << flav << " was " << istagged << " is " << newtag << endl;
    

    if(newtag){
      nBTagsSub += 1;
    }

    if(abs(flav)==1||abs(flav)==2||abs(flav)==3||abs(flav)==21){
      if(fabs(subeta)>=0&&fabs(subeta)<=2.4){

	delete csvu;
	delete csvd;

      }
    }

    delete csv;
    delete numpt;
    delete denpt;
    delete effipt;
    delete errbc;
    delete rand;
    
  }
  
  if(dosf){
    file_mc->Close();   
    delete file_mc;
  }

  return nBTagsSub;

}



// int subJetBTag(TopJet topjet, E_BtagType type){
//   int nBTagsSub = 0;
//   float discriminator_cut;
//   if(type==e_CSVL) discriminator_cut = 0.244;
//   if(type==e_CSVM) discriminator_cut = 0.679;
//   if(type==e_CSVT) discriminator_cut = 0.898;

//   //Create a vector of subjets
//   std::vector<Particle> subjets_top;
//   //Create a float vector of the subjets discriminators
//   std::vector<float> btagsub_combinedSecondaryVertex_top;

//   //Fill the vector of subjets with the subjets of a topjet
//   subjets_top=topjet.subjets();
//   //Fill the vector of discriminators with the discriminators of the subjets of a certain topjet
//   btagsub_combinedSecondaryVertex_top=topjet.btagsub_combinedSecondaryVertex();

//   //Looping over subjets and checking if they are b-tagged
//   for(unsigned int i=0; i < btagsub_combinedSecondaryVertex_top.size(); ++i){
//     float test=btagsub_combinedSecondaryVertex_top[i];
//     if(test>discriminator_cut){
//       nBTagsSub += 1;
//       //This means it is b-tagged
//     }
//   }
//   return nBTagsSub;
// }


bool HiggsTag(TopJet topjet, E_BtagType type1, E_BtagType type2, TString mode, TString filename){
  int nBTagsSub1 = 0;
  int nBTagsSub2 = 0;
  
  nBTagsSub1 = subJetBTag(topjet, type1, mode, filename);
  nBTagsSub2 = subJetBTag(topjet, type2, mode, filename);

 
  if (type1 == type2 &&  nBTagsSub1>= 2) return true;
  if (type1 > type2 && nBTagsSub1!=0 && nBTagsSub2 >=2) return true;
  if (type1 < type2 && nBTagsSub1 >= 2 && nBTagsSub2 != 0)return true;
  else return false;
}







// bool HiggsTag(TopJet topjet, E_BtagType type1, E_BtagType type2){
//   int nBTagsSub1 = 0;
//   int nBTagsSub2 = 0;
//   float discriminator_cut1;
//   float discriminator_cut2;
//   if(type1==e_CSVL) discriminator_cut1 = 0.244;
//   if(type1==e_CSVM) discriminator_cut1 = 0.679;
//   if(type1==e_CSVT) discriminator_cut1 = 0.898;
//   if(type2==e_CSVL) discriminator_cut2 = 0.244;
//   if(type2==e_CSVM) discriminator_cut2 = 0.679;
//   if(type2==e_CSVT) discriminator_cut2 = 0.898;
//   // std::cout << "discriminator_cut1: " <<  discriminator_cut1 << " discriminator_cut2: "<< discriminator_cut1 << std::endl;
//   //Create a vector of subjets
//   std::vector<Particle> subjets_top;
//   //Create a float vector of the subjets discriminators
//   std::vector<float> btagsub_combinedSecondaryVertex_top;

//   //Fill the vector of subjets with the subjets of a topjet
//   subjets_top=topjet.subjets();
//   //Fill the vector of discriminators with the discriminators of the subjets of a certain topjet
//   btagsub_combinedSecondaryVertex_top=topjet.btagsub_combinedSecondaryVertex();

//   //Looping over subjets and checking if they are b-tagged
//   for(unsigned int i=0; i < btagsub_combinedSecondaryVertex_top.size(); ++i){
//     float test=btagsub_combinedSecondaryVertex_top[i];
//     if (nBTagsSub1 != 0 && test>discriminator_cut2) nBTagsSub2 =+ 1;
//     if(test>discriminator_cut1){
//       if(test>discriminator_cut2 && nBTagsSub2==0) nBTagsSub2+=1;
//       else nBTagsSub1 += 1;      //This means it is b-tagged  
 
//     }
//   }
//   if (nBTagsSub1!=0 && nBTagsSub2!=0) return true;
//   else return false;
// }


bool HepTopTagFull(TopJet topjet, std::vector<PFParticle>* allparts){

  //Transform the SFrame TopJet object in a fastjet::PseudoJet

  if (!allparts) return false;

  fastjet::ClusterSequence* JetFinder;
  fastjet::JetDefinition* JetDef ;

  JetProps jp(&topjet, allparts);
  std::vector<fastjet::PseudoJet> jetpart = jp.GetJetConstituents();
 
  //Clustering definition
  double conesize=3;
  JetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm,conesize); 

  JetFinder = new fastjet::ClusterSequence(jetpart, *JetDef);

  std::vector<fastjet::PseudoJet> tops = JetFinder->inclusive_jets(10.);

  if (tops.size() != 1){
    std::cout << "Problem! Doesn't give exactly one jet!!!!!!!!!!!!!!Gives " << tops.size() << " jets" << std::endl;
    delete JetFinder;
    delete JetDef;
    return false;
  }

  std::vector<fastjet::PseudoJet> SortedJets = sorted_by_pt(tops);


  //Run the HEPTopTagger
  external::HEPTopTagger tagger(*JetFinder, SortedJets[0]);

  //Mass window to be applied in a second step
  tagger.set_top_range(0.0, 10000.0);
  tagger.set_mass_drop_threshold(0.8);
  tagger.set_max_subjet_mass(30);

  tagger.run_tagger();

  delete JetFinder;
  delete JetDef;
 
  if (tagger.is_masscut_passed()) return true;
  else return false;

  return true;
 
}



// global function to define a tagged jet

bool IsTagged(Jet & jet, E_BtagType type)
{
    if(type==e_CSVL && jet.btag_combinedSecondaryVertex()>0.244) return true;
    if(type==e_CSVM && jet.btag_combinedSecondaryVertex()>0.679) return true;
    if(type==e_CSVT && jet.btag_combinedSecondaryVertex()>0.898) return true;
    if(type==e_JPL && jet.btag_jetProbability()>0.275) return true;
    if(type==e_JPM && jet.btag_jetProbability()>0.545) return true;
    if(type==e_JPT && jet.btag_jetProbability()>0.790) return true;       
    
    return false;
}

//variable HEP Tagger from Rebekka

bool variableHepTopTagWithMatch(TopJet topjet, double ptJetMin, double massWindowLower, double massWindowUpper, double cutCondition2, double cutCondition3)
{

  //Taking the top tag from the proper jet collection

  EventCalc* calc = EventCalc::Instance();
  
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  
  double deltarmin = double_infinity();
  
  TopJet nextjet;
  
  for(unsigned int it=0; it<bcc->toptagjets->size();++it){

    TopJet top4jet=bcc->toptagjets->at(it);
    
    if(top4jet.deltaR(topjet) < deltarmin){
      deltarmin = top4jet.deltaR(topjet);
      nextjet = top4jet;
    }
    
  }
  
  if(deltarmin>=0.3) return 0;

  double mjet;
  double ptjet;
  int nsubjets;
  
  double topmass=172.3;
  double wmass=80.4;
  
  nsubjets=nextjet.numberOfDaughters();
  
  LorentzVector allsubjets(0,0,0,0);
  
  for(int j=0; j<nextjet.numberOfDaughters(); ++j) {
    allsubjets += nextjet.subjets()[j].v4();
  }
  if(!allsubjets.isTimelike()) {
    mjet=0;
    return false;
  }
  
  mjet = allsubjets.M();
  ptjet= nextjet.pt();
    
  double m12, m13, m23;
  
  //The subjetcs have to be three
  if(nsubjets==3) {
    
    std::vector<Particle> subjets = nextjet.subjets();
    sort(subjets.begin(), subjets.end(), HigherPt());
    
    m12 = 0;
    if( (subjets[0].v4()+subjets[1].v4()).isTimelike())
      m12=(subjets[0].v4()+subjets[1].v4()).M();
    m13 = 0;
    if( (subjets[0].v4()+subjets[2].v4()).isTimelike() )
      m13=(subjets[0].v4()+subjets[2].v4()).M();
    m23 = 0;
    if( (subjets[1].v4()+subjets[2].v4()).isTimelike()  )
      m23 = (subjets[1].v4()+subjets[2].v4()).M();
    
  } else {
    return false;
  }
  
  double rmin=massWindowLower*wmass/topmass;
  double rmax=massWindowUpper*wmass/topmass;
  
  int keep=0;
  
  //Conditions on the subjects: at least one has to be true
  //1 condition
  if(atan(m13/m12)>0.2 && atan(m13/m12)<1.3 && m23/mjet>rmin && m23/mjet<rmax) keep=1;
  
  //2 condition
  double cond2left=pow(rmin,2)*(1+pow((m13/m12),2));
  double cond2cent=1-pow(m23/mjet,2);
  double cond2right=pow(rmax,2)*(1+pow(m13/m12,2));
  
  if(cond2left<cond2cent && cond2cent<cond2right && m23/mjet>cutCondition2) keep=1;
  
  //3 condition
  double cond3left=pow(rmin,2)*(1+pow((m12/m13),2));
  double cond3cent=1-pow(m23/mjet,2);
  double cond3right=pow(rmax,2)*(1+pow(m12/m13,2));
  
  if(cond3left<cond3cent && cond3cent<cond3right && m23/mjet>cutCondition3) keep=1;
  
  if( nextjet.v4().M() < 140 || nextjet.v4().M() > 250) keep=0;
  
  //Final requirement: at least one of the three subjets conditions and total pt
  if(keep==1 && ptjet>ptJetMin) {
    return true;
  } else {
    return false;
  }
  
}



//variable HEP Tagger from Rebekka

bool variableHepTopTag(TopJet topjet, double ptJetMin, double massWindowLower, double massWindowUpper, double cutCondition2, double cutCondition3)
{

  TopJet nextjet=topjet;
  
  double mjet;
  double ptjet;
  int nsubjets;
  
  double topmass=172.3;
  double wmass=80.4;
  
  nsubjets=nextjet.numberOfDaughters();
  
  LorentzVector allsubjets(0,0,0,0);
  
  for(int j=0; j<nextjet.numberOfDaughters(); ++j) {
    allsubjets += nextjet.subjets()[j].v4();
  }
  if(!allsubjets.isTimelike()) {
    mjet=0;
    return false;
  }
  
  mjet = allsubjets.M();
  ptjet= nextjet.pt();
    
  double m12, m13, m23;
  
  //The subjetcs have to be three
  if(nsubjets==3) {
    
    std::vector<Particle> subjets = nextjet.subjets();
    sort(subjets.begin(), subjets.end(), HigherPt());
    
    m12 = 0;
    if( (subjets[0].v4()+subjets[1].v4()).isTimelike())
      m12=(subjets[0].v4()+subjets[1].v4()).M();
    m13 = 0;
    if( (subjets[0].v4()+subjets[2].v4()).isTimelike() )
      m13=(subjets[0].v4()+subjets[2].v4()).M();
    m23 = 0;
    if( (subjets[1].v4()+subjets[2].v4()).isTimelike()  )
      m23 = (subjets[1].v4()+subjets[2].v4()).M();
    
  } else {
    return false;
  }
  
  double rmin=massWindowLower*wmass/topmass;
  double rmax=massWindowUpper*wmass/topmass;
  
  int keep=0;
  
  //Conditions on the subjects: at least one has to be true
  //1 condition
  if(atan(m13/m12)>0.2 && atan(m13/m12)<1.3 && m23/mjet>rmin && m23/mjet<rmax) keep=1;
  
  //2 condition
  double cond2left=pow(rmin,2)*(1+pow((m13/m12),2));
  double cond2cent=1-pow(m23/mjet,2);
  double cond2right=pow(rmax,2)*(1+pow(m13/m12,2));
  
  if(cond2left<cond2cent && cond2cent<cond2right && m23/mjet>cutCondition2) keep=1;
  
  //3 condition
  double cond3left=pow(rmin,2)*(1+pow((m12/m13),2));
  double cond3cent=1-pow(m23/mjet,2);
  double cond3right=pow(rmax,2)*(1+pow(m12/m13,2));
  
  if(cond3left<cond3cent && cond3cent<cond3right && m23/mjet>cutCondition3) keep=1;
  
  if( mjet < 140 || mjet > 250) keep=0;
  
  //Final requirement: at least one of the three subjets conditions and total pt
  if(keep==1 && ptjet>ptJetMin) {
    return true;
  } else {
    return false;
  }
  
}






bool HepTopTagInverted(TopJet topjet)
{
  //Taking the top tag from the proper jet collection

  EventCalc* calc = EventCalc::Instance();
  
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  
  double deltarmin = double_infinity();
  
  TopJet nextjet;
  
  for(unsigned int it=0; it<bcc->toptagjets->size();++it){

    TopJet top4jet=bcc->toptagjets->at(it);
    
    if(top4jet.deltaR(topjet) < deltarmin){
      deltarmin = top4jet.deltaR(topjet);
      nextjet = top4jet;
    }
    
  }
  
  if(deltarmin>=0.3) return 0;


    double mjet;
    double ptjet;
    int nsubjets;

    double topmass=172.3;
    double wmass=80.4;

    nsubjets=nextjet.numberOfDaughters();

    LorentzVector allsubjets(0,0,0,0);

    for(int j=0; j<nextjet.numberOfDaughters(); ++j) {
        allsubjets += nextjet.subjets()[j].v4();
    }
    if(!allsubjets.isTimelike()) {
        mjet=0;
        return false;
    }

    mjet = allsubjets.M();
    ptjet= nextjet.pt();

    double m12, m13, m23;

    //The subjetcs have to be three
    if(nsubjets==3) {

        std::vector<Particle> subjets = nextjet.subjets();
        sort(subjets.begin(), subjets.end(), HigherPt());

        m12 = 0;
        if( (subjets[0].v4()+subjets[1].v4()).isTimelike())
            m12=(subjets[0].v4()+subjets[1].v4()).M();
        m13 = 0;
        if( (subjets[0].v4()+subjets[2].v4()).isTimelike() )
            m13=(subjets[0].v4()+subjets[2].v4()).M();
        m23 = 0;
        if( (subjets[1].v4()+subjets[2].v4()).isTimelike()  )
            m23 = (subjets[1].v4()+subjets[2].v4()).M();

    } else {
        return false;
    }

    double rmin=0.85*wmass/topmass;
    double rmax=1.15*wmass/topmass;

    int keep=0;


    if (m23/mjet <= 0.35) keep = 1; 

    //invert top mass window
    if( mjet > 140 && mjet < 250) keep=0;


    //Final requirement: at least one of the three subjets conditions and total pt
    if(keep==1 && ptjet>200.) {
        return true;
    } else {
        return false;
    }

}




//HEP Tagger from Ivan

bool HepTopTag(TopJet topjet)
{

  //call variable tagger with default parameters
  return variableHepTopTag(topjet);

}

//HEP Tagger from Ivan

bool HepTopTagWithMatch(TopJet topjet)
{

  //call variable tagger with default parameters
  return variableHepTopTagWithMatch(topjet);

}

//default values (mminLower=50., mjetLower=140, mjetUpper=250.) defined in Utils.h
bool variableTopTag(TopJet topjet, double &mjet, int &nsubjets, double &mmin, double mminLower, double mjetLower, double mjetUpper)
{

    nsubjets=topjet.numberOfDaughters();

    LorentzVector allsubjets(0,0,0,0);

    for(int j=0; j<topjet.numberOfDaughters(); ++j) {
        allsubjets += topjet.subjets()[j].v4();
    }
    if(!allsubjets.isTimelike()) {
        mjet=0;
        mmin=0;
//  mminLower=50;
//   mjetLower=140;
//   mjetUpper=250;
        return false;
    }

    mjet = allsubjets.M();

    if(nsubjets>=3) {

        std::vector<Particle> subjets = topjet.subjets();
        sort(subjets.begin(), subjets.end(), HigherPt());

        double m01 = 0;
        if( (subjets[0].v4()+subjets[1].v4()).isTimelike())
            m01=(subjets[0].v4()+subjets[1].v4()).M();
        double m02 = 0;
        if( (subjets[0].v4()+subjets[2].v4()).isTimelike() )
            m02=(subjets[0].v4()+subjets[2].v4()).M();
        double m12 = 0;
        if( (subjets[1].v4()+subjets[2].v4()).isTimelike() )
            m12 = (subjets[1].v4()+subjets[2].v4()).M();

        //minimum pairwise mass
        mmin = std::min(m01,std::min(m02,m12));
    }

    //at least 3 sub-jets
    if(nsubjets<3) return false;
    //minimum pairwise mass > 50 GeV/c^2
    if(mmin<mminLower) return false;
    //jet mass between 140 and 250 GeV/c^2
    if(mjet<mjetLower || mjet>mjetUpper) return false;

    return true;
}


bool TopTag(TopJet topjet,  double &mjet, int &nsubjets, double &mmin)
{
  //call variable tagger with default parameters
  return variableTopTag(topjet, mjet, nsubjets, mmin);

}
 
Jet* nextJet(const Particle *p, std::vector<Jet> *jets)
{

    double deltarmin = double_infinity();
    Jet* nextjet=0;
    for(unsigned int i=0; i<jets->size(); ++i) {
        Jet ji = jets->at(i);
	if (fabs(p->pt() - ji.pt())<1e-8) continue; // skip identical particle
        if(jets->at(i).deltaR(*p) < deltarmin) {
            deltarmin = jets->at(i).deltaR(*p);
            nextjet = &jets->at(i);
        }
    }

    return nextjet;
}

bool WTag(TopJet prunedjet,  double& mjet, int &nsubjets, double& massdrop)
{

    nsubjets=prunedjet.numberOfDaughters();

    mjet = 0;
    if(prunedjet.v4().isTimelike())
        mjet = prunedjet.v4().M();

    //calculate mass drop for first sub-jet ordered in pt
    massdrop = 0;
    if(nsubjets>=1 && mjet>0) {

        std::vector< Particle > subjets = prunedjet.subjets();
        sort(subjets.begin(), subjets.end(), HigherPt());

        double m1 = 0;
        if(subjets[0].v4().isTimelike())
            m1 = subjets[0].v4().M();

        massdrop = m1/mjet;
    }

    //at least 2 sub-jets
    if(nsubjets<2) return false;
    //60 GeV < pruned jet mass < 100 GeV
    if(mjet <= 60 || mjet >= 100) return false;
    //mass drop < 0.4
    if(massdrop>=0.4) return false;

    return true;

}

float relIsoMuon(EventCalc & event, const Muon & mu, float deltaR){
  float chargedHadronIso=0;
  float neutralHadronIso=0;
  float photonIso=0;
  float puiso=0;

  for(PFParticle & pfp : *event.GetIsoPFParticles()){
      auto dr = pfp.deltaR(mu);
      if(dr < deltaR){
         if(pfp.particleID() == PFParticle::eH && pfp.pt()>0.0 && dr > 0.0001 ) chargedHadronIso += pfp.pt();
         if(pfp.particleID() == PFParticle::eH0 && pfp.pt()>0.5 && dr > 0.01) neutralHadronIso += pfp.pt();
         if(pfp.particleID() == PFParticle::eGamma && pfp.pt()>0.5 && dr > 0.01) photonIso += pfp.pt();
      }
  }
  
  for(PFParticle & pfp : *event.GetPUIsoPFParticles()){
    auto dr = pfp.deltaR(mu);
    if(dr<deltaR ){
      if(pfp.particleID() == PFParticle::eH && pfp.pt()>0.5 && dr>0.01 ) puiso += pfp.pt();
    }
  }
  
  return (chargedHadronIso + std::max( 0.0, neutralHadronIso + photonIso - 0.5*puiso))/mu.pt();
}

float relIsoMuon(const Muon & mu, float deltaR ){
     return relIsoMuon(*EventCalc::Instance(), mu, deltaR);
}

double pTrel(const Particle *p, std::vector<Jet> *jets)
{

    double ptrel=0;
    Jet* nextjet =  nextJet(p,jets);
    if (!nextjet) return ptrel;

    TVector3 p3(p->v4().Px(),p->v4().Py(),p->v4().Pz());
    TVector3 jet3(nextjet->v4().Px(),nextjet->v4().Py(),nextjet->v4().Pz());

    if(p3.Mag()!=0 && jet3.Mag()!=0) {
        double sin_alpha = (p3.Cross(jet3)).Mag()/p3.Mag()/jet3.Mag();
        ptrel = p3.Mag()*sin_alpha;
    } else {
        std::cout << "something strange happend in the ptrel calculation: either lepton or jet momentum is 0" <<std::endl;
    }

    return ptrel;
}

double deltaRmin(const Particle *p, std::vector<Jet> *jets)
{
    Jet* j = nextJet(p,jets);
    double dr = 999.;
    if (j) dr = j->deltaR(*p);
    return dr;
}

TVector3 toVector(LorentzVector v4)
{

    TVector3 v3(0,0,0);
    v3.SetX(v4.X());
    v3.SetY(v4.Y());
    v3.SetZ(v4.Z());
    return v3;
}

TVector3 toVector(LorentzVectorXYZE v4)
{

    TVector3 v3(0,0,0);
    v3.SetX(v4.X());
    v3.SetY(v4.Y());
    v3.SetZ(v4.Z());
    return v3;
}

LorentzVectorXYZE toXYZ(LorentzVector v4)
{

    LorentzVectorXYZE v4_new(0,0,0,0);
    v4_new.SetPx(v4.Px());
    v4_new.SetPy(v4.Py());
    v4_new.SetPz(v4.Pz());
    v4_new.SetE(v4.E());
    return v4_new;
}

LorentzVector toPtEtaPhi(LorentzVectorXYZE v4)
{

    LorentzVector v4_new(0,0,0,0);
    v4_new.SetPt(v4.Pt());
    v4_new.SetEta(v4.Eta());
    v4_new.SetPhi(v4.Phi());
    v4_new.SetE(v4.E());
    return v4_new;
}

double deltaR(LorentzVector v1, LorentzVector v2)
{

    Particle p1;
    p1.set_v4(v1);
    Particle p2;
    p2.set_v4(v2);
    return p1.deltaR(p2);
}

double double_infinity()
{
    return std::numeric_limits<double>::infinity() ;
}

int int_infinity()
{
    return std::numeric_limits<int>::max() ;
}

int myPow(int x, unsigned int p)
{
    int i = 1;
    for (unsigned int j = 1; j <= p; j++)  i *= x;
    return i;
}


double pTrel(const Particle & p1, const Particle & p2){
  TVector3 p3(p1.v4().Px(),p1.v4().Py(),p1.v4().Pz());
  TVector3 p4(p2.v4().Px(),p2.v4().Py(),p2.v4().Pz());
  
  if(p3.Mag()!=0 && p4.Mag()!=0) {
    double sin_alpha = (p3.Cross(p4)).Mag()/p3.Mag()/p4.Mag();
    return p3.Mag()*sin_alpha;
  } else {
    throw runtime_error("pTrel: something strange happend in the ptrel calculation: either lepton or jet momentum is 0");
  }
}
 
double pTrel(const LorentzVector & p1, const LorentzVector & p2){
  TVector3 p3(p1.Px(),p1.Py(),p1.Pz());
  TVector3 p4(p2.Px(),p2.Py(),p2.Pz());
  
  if(p3.Mag()!=0 && p4.Mag()!=0) {
    double sin_alpha = (p3.Cross(p4)).Mag()/p3.Mag()/p4.Mag();
    return p3.Mag()*sin_alpha;
  } else {
    throw runtime_error("something strange happend in the ptrel calculation: either lepton or jet momentum is 0");
  }
}

TableOutput::TableOutput(vector<std::string> header_): ncols(header_.size()), header(std::move(header_)){
}

void TableOutput::print(ostream & out){
    const string sep = " | ";
    // calculate column sizes, given by maximum entry:
    vector<size_t> colsize(ncols, 0);
    for(size_t i=0; i<ncols; ++i){
        colsize[i] = header[i].size();
    }
    for(auto & row: rows){
        for(size_t i=0; i<ncols; ++i){
            colsize[i] = max(colsize[i], row[i].size());
        }
    }
    size_t total_width = 0;
    for(size_t i=0; i<ncols; ++i){
        total_width += sep.size() + colsize[i];
    }
    total_width += sep.size();
    
    // output helpers for a horizontal line and for a single row:
    auto hline = [&](){
            out << " ";
            for(size_t i=1; i+1<total_width; ++i){
                out << "-";
            }
            out << endl;};
    auto out_row = [&](const vector<string> & row){
        for(size_t i=0; i<ncols; ++i){
            out << sep << row[i];
            size_t nfill = colsize[i] - row[i].size();
            for(size_t j=0; j<nfill; ++j){
                out << " ";
            }
        }
        out << sep << endl;
    };

    // output table:
    hline();
    out_row(header);
    hline();
    for(auto & row : rows){
        out_row(row);
    }
    hline();
}

// entries.size() == ncols must hold
void TableOutput::add_row(vector<string> row){
    assert(ncols == row.size());
    rows.push_back(std::move(row));
}


