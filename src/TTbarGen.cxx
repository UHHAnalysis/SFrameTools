#include "include/TTbarGen.h"

using namespace std;

TTbarGen::TTbarGen()
{
  ObjectHandler* objs = ObjectHandler::Instance();
  BaseCycleContainer* bcc = objs->GetBaseCycleContainer();
  m_pdgId1 = 0;
  m_pdgId2 = 0;
  for(unsigned int i=0; i<bcc->genparticles->size(); ++i)
   {
     GenParticle genp = bcc->genparticles->at(i);

     if ( genp.pdgId() == 6 ) 
       {
	 m_Top = genp;
	 if (genp.daughter(bcc->genparticles,1) && genp.daughter(bcc->genparticles,2))
	   {
	     if (genp.daughter(bcc->genparticles,1)->pdgId() == 24) 
	       {
		 m_indexW = genp.daughter(bcc->genparticles,1)->index();
		 m_WTop = bcc->genparticles->at(m_indexW);
		 m_indexb = genp.daughter(bcc->genparticles,2)->index();
		 m_bTop = bcc->genparticles->at(m_indexb);
	       }
	     if (genp.daughter(bcc->genparticles,2)->pdgId() == 24)
	       {
		 m_indexW = genp.daughter(bcc->genparticles,2)->index();
		 m_WTop = bcc->genparticles->at(m_indexW);
		 m_indexb = genp.daughter(bcc->genparticles,1)->index();
		 m_bTop = bcc->genparticles->at(m_indexb);
	       }
	   }
       }
     if ( genp.pdgId() == -6 ) 
       {
	 m_Antitop = genp;
	 if (genp.daughter(bcc->genparticles,1) && (genp.daughter(bcc->genparticles,2)))
	   {
	     if (genp.daughter(bcc->genparticles,1)->pdgId() == -24) 
	       {
		 m_indexW = genp.daughter(bcc->genparticles,1)->index();
		 m_WAntitop = bcc->genparticles->at(m_indexW);
		 m_indexb = genp.daughter(bcc->genparticles,2)->index();
		 m_bAntitop = bcc->genparticles->at(m_indexb);
	       }
	     if (genp.daughter(bcc->genparticles,2)->pdgId() == -24)
	       {
		 m_indexW = genp.daughter(bcc->genparticles,2)->index();
		 m_WAntitop = bcc->genparticles->at(m_indexW);
		 m_indexb = genp.daughter(bcc->genparticles,1)->index();
		 m_bAntitop = bcc->genparticles->at(m_indexb);
	       }
	   }
       }
     if ( genp.pdgId() == 24 && genp.mother(bcc->genparticles,1)->pdgId() == 6)
       {
	 if (genp.daughter(bcc->genparticles,1) && genp.daughter(bcc->genparticles,2))
	   {
	     m_index1 = genp.daughter(bcc->genparticles,1)->index();
	     m_Wdecay1 = bcc->genparticles->at(m_index1);
	     m_pdgId1= m_Wdecay1.pdgId();
	     m_index2 = genp.daughter(bcc->genparticles,2)->index();
	     m_Wdecay2 = bcc->genparticles->at(m_index2);
	     
	   }
       }
     if ( genp.pdgId() == -24 && genp.mother(bcc->genparticles,1)->pdgId() == -6)
       {
	 if (genp.daughter(bcc->genparticles,1) && (genp.daughter(bcc->genparticles,2)))
	   {
	     m_index1 = genp.daughter(bcc->genparticles,1)->index();
	     m_WMinusdecay1 = bcc->genparticles->at(m_index1);
	     m_pdgId2= m_WMinusdecay1.pdgId(); 
	     m_index2 = genp.daughter(bcc->genparticles,2)->index();
	     m_WMinusdecay2 = bcc->genparticles->at(m_index2);
	   }
       }
   }

  // W not linked correctly -> W decay products have top as mother
  if(m_pdgId1==0 || m_pdgId2==0 ){
    //search all particles with top as mother; store that ones which are not W or b (or equivalent light quark)
    bool notfilled1=true;
    bool notfilled2=true;

    for(unsigned int i=0; i<bcc->genparticles->size(); ++i)
      {
	GenParticle genp = bcc->genparticles->at(i);
	

	if(genp.mother(bcc->genparticles,1)){
	  if(genp.mother(bcc->genparticles,1)->pdgId()==6 && abs(genp.pdgId())!=24 &&  genp.pdgId()!=1 && genp.pdgId()!=3 && genp.pdgId()!=5 ){
	    
	    if(notfilled1){
	      m_Wdecay1 = genp;
	      m_pdgId1 = genp.pdgId();
	      notfilled1=false;
	    }
	    else{
	      m_Wdecay2 = genp;
	    }
	  }
	  

	  if(genp.mother(bcc->genparticles,1)->pdgId()==-6 && abs(genp.pdgId())!=24 && genp.pdgId()!=-1 && genp.pdgId()!=-3 && genp.pdgId()!=-5 ){
	    if(notfilled2){
	      m_WMinusdecay1 = genp;
	      m_pdgId2 = genp.pdgId();
	      notfilled2=false;
	    }
	    else{
	      m_WMinusdecay2 = genp;
	    }
	  }
	  if(genp.mother(bcc->genparticles,1)->pdgId()==6 && abs(genp.pdgId())!=24 &&  (genp.pdgId()==1 || genp.pdgId()==3 || genp.pdgId()==5) ){
	    m_bTop = genp;
	  }
	  if(genp.mother(bcc->genparticles,1)->pdgId()==-6 && abs(genp.pdgId())!=24 &&  (genp.pdgId()==-1 || genp.pdgId()==-3 || genp.pdgId()==-5) ){
	    m_bAntitop = genp;
	  }
	}

      }
    //fill the missing W bosons
    GenParticle Wplus, Wminus;

    LorentzVector Wp, Wm;
    Wp.SetPxPyPzE( Wdecay1().v4().Px()+ Wdecay2().v4().Px(),  Wdecay1().v4().Py()+ Wdecay2().v4().Py(),  Wdecay1().v4().Pz()+ Wdecay2().v4().Pz(),  Wdecay1().v4().E()+ Wdecay2().v4().E() );
    Wm.SetPxPyPzE( WMinusdecay1().v4().Px()+ WMinusdecay2().v4().Px(),  WMinusdecay1().v4().Py()+ WMinusdecay2().v4().Py(),  WMinusdecay1().v4().Pz()+ WMinusdecay2().v4().Pz(),  WMinusdecay1().v4().E()+ WMinusdecay2().v4().E() );

    Wplus.set_v4( Wp );
    Wplus.set_charge(1.);
    Wplus.set_pdgId(24);
    Wplus.set_status(3);
    Wminus.set_v4( Wm );
    Wminus.set_charge(-1.);
    Wminus.set_pdgId(-24);
    Wminus.set_status(3);
    m_WTop = Wplus;
    m_WAntitop = Wminus;
  }


}   
 

TTbarGen::~TTbarGen()
{
  // default destructor, does nothing
}

GenParticle TTbarGen::Top()
{
  return m_Top;
}

GenParticle TTbarGen::Antitop()
{
  return m_Antitop;
} 

GenParticle TTbarGen::WTop()
{
  return m_WTop;
}

GenParticle TTbarGen::WAntitop()
{
  return m_WAntitop;
}

GenParticle TTbarGen::bTop()
{
  return m_bTop;
}

GenParticle TTbarGen::bAntitop()
{
  return m_bAntitop;
} 

GenParticle TTbarGen::Wdecay1()
{
  return m_Wdecay1;
} 

GenParticle TTbarGen::Wdecay2()
{
  return m_Wdecay2;
} 

GenParticle TTbarGen::WMinusdecay1()
{
  return m_WMinusdecay1;
} 

GenParticle TTbarGen::WMinusdecay2()
{
  return m_WMinusdecay2;
} 

TTbarGen::E_DecayChannel TTbarGen::DecayChannel()
{
  if ((abs(m_pdgId1)==11 || abs(m_pdgId1)==12) && (abs(m_pdgId2)==11 || abs(m_pdgId2)==12)) m_type = e_ee;
  if ((abs(m_pdgId1)==13 || abs(m_pdgId1)==14) && (abs(m_pdgId2)==13 || abs(m_pdgId2)==14)) m_type = e_mumu;
  if ((abs(m_pdgId1)==15 || abs(m_pdgId1)==16) && (abs(m_pdgId2)==15 || abs(m_pdgId2)==16)) m_type = e_tautau;
  if ((abs(m_pdgId1)<=5) && (abs(m_pdgId2)<=5)) m_type = e_had;
  if ((abs(m_pdgId1)==11 || abs(m_pdgId1)==12) && (abs(m_pdgId2)<=5)) m_type = e_ehad;
  if ((abs(m_pdgId2)==11 || abs(m_pdgId2)==12) && (abs(m_pdgId1)<=5)) m_type = e_ehad;
  if ((abs(m_pdgId1)==13 || abs(m_pdgId1)==14) && (abs(m_pdgId2)<=5)) m_type = e_muhad;
  if ((abs(m_pdgId2)==13 || abs(m_pdgId2)==14) && (abs(m_pdgId1)<=5)) m_type = e_muhad;
  if ((abs(m_pdgId1)==15 || abs(m_pdgId1)==16) && (abs(m_pdgId2)<=5)) m_type = e_tauhad;
  if ((abs(m_pdgId2)==15 || abs(m_pdgId2)==16) && (abs(m_pdgId1)<=5)) m_type = e_tauhad;
  if ((abs(m_pdgId1)==11 || abs(m_pdgId1)==12) && (abs(m_pdgId2)==13 || abs(m_pdgId2)==14)) m_type = e_emu;
  if ((abs(m_pdgId2)==11 || abs(m_pdgId2)==12) && (abs(m_pdgId1)==13 || abs(m_pdgId1)==14)) m_type = e_emu;
  if ((abs(m_pdgId1)==11 || abs(m_pdgId1)==12) && (abs(m_pdgId2)==15 || abs(m_pdgId2)==16)) m_type = e_etau;
  if ((abs(m_pdgId2)==11 || abs(m_pdgId2)==12) && (abs(m_pdgId1)==15 || abs(m_pdgId1)==16)) m_type = e_etau;
  if ((abs(m_pdgId1)==13 || abs(m_pdgId1)==14) && (abs(m_pdgId2)==15 || abs(m_pdgId2)==16)) m_type = e_mutau;
  if ((abs(m_pdgId2)==13 || abs(m_pdgId2)==14) && (abs(m_pdgId1)==15 || abs(m_pdgId1)==16)) m_type = e_mutau;
  return m_type;
}

bool TTbarGen::IsTopHadronicDecay()
{
  if (abs(m_pdgId1)<=5) return true;
  else return false;
}

bool TTbarGen::IsAntiTopHadronicDecay()
{
  if (abs(m_pdgId2)<=5) return true;
  else return false;
}

GenParticle TTbarGen::ChargedLepton(){
  GenParticle lepton;  
  if (m_type != e_ehad &&  m_type != e_muhad  && m_type!= e_tauhad){
    std::cerr << "This is no l+jets ttbar event, no charged lepton found" <<std::endl;
    return lepton;
  }

  int nlepton=0;
  if(abs(Wdecay1().pdgId())==11 || abs(Wdecay1().pdgId())==13 || abs(Wdecay1().pdgId())==15){
    lepton = Wdecay1();
    nlepton++;
  }
  if(abs(Wdecay2().pdgId())==11 || abs(Wdecay2().pdgId())==13 || abs(Wdecay2().pdgId())==15){
    lepton = Wdecay2();
    nlepton++;
  }
  if(abs(WMinusdecay1().pdgId())==11 || abs(WMinusdecay1().pdgId())==13 || abs(WMinusdecay1().pdgId())==15){
    lepton = WMinusdecay1();
    nlepton++;
  }
  if(abs(WMinusdecay2().pdgId())==11 || abs(WMinusdecay2().pdgId())==13 || abs(WMinusdecay2().pdgId())==15){
    lepton = WMinusdecay2();
    nlepton++;
  }
  if(nlepton!=1) std::cerr << "Not exactly one lepton found " << nlepton<< std::endl;

  return lepton;

}

GenParticle TTbarGen::Neutrino(){
  GenParticle neutrino;
  
  if (m_type != e_ehad &&  m_type != e_muhad  && m_type!= e_tauhad){
    std::cerr << "This is no l+jets ttbar event, no neutrino found" <<std::endl;
    return neutrino;
  }

  int nneutrino=0;
  if(abs(Wdecay1().pdgId())==12 || abs(Wdecay1().pdgId())==14 || abs(Wdecay1().pdgId())==16){
    neutrino = Wdecay1();
    nneutrino++;
  }
  if(abs(Wdecay2().pdgId())==12 || abs(Wdecay2().pdgId())==14 || abs(Wdecay2().pdgId())==16){
    neutrino = Wdecay2();
    nneutrino++;
  }
  if(abs(WMinusdecay1().pdgId())==12 || abs(WMinusdecay1().pdgId())==14 || abs(WMinusdecay1().pdgId())==16){
    neutrino = WMinusdecay1();
    nneutrino++;
  }
  if(abs(WMinusdecay2().pdgId())==12 || abs(WMinusdecay2().pdgId())==14 || abs(WMinusdecay2().pdgId())==16){
    neutrino = WMinusdecay2();
    nneutrino++;
  }
  if(nneutrino!=1) std::cerr << "Not exactly one neutrino found " << nneutrino<< std::endl;

  return neutrino;

}

GenParticle TTbarGen::TopLep(){
 
  if(ChargedLepton().charge()>0) return Top();
  else return Antitop();

}

GenParticle TTbarGen::TopHad(){

  if(ChargedLepton().charge()<0) return Top();
  else return Antitop();

}


GenParticle TTbarGen::BLep(){
 
  if(ChargedLepton().charge()>0) return bTop();
  else return bAntitop();

}

GenParticle TTbarGen::BHad(){
 
  if(ChargedLepton().charge()<0) return bTop();
  else return bAntitop();

}

GenParticle TTbarGen::WLep(){
 
  if(ChargedLepton().charge()>0) return WTop();
  else return WAntitop();

}

GenParticle TTbarGen::WHad(){
 
  if(ChargedLepton().charge()<0) return WTop();
  else return WAntitop();

}

GenParticle TTbarGen::Q1(){
 
  if(ChargedLepton().charge()>0) return WMinusdecay1();
  else return Wdecay1();

}

GenParticle TTbarGen::Q2(){
 
  if(ChargedLepton().charge()>0) return WMinusdecay2();
  else return Wdecay2();

}

