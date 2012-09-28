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
	     m_index2 = genp.daughter(bcc->genparticles,1)->index();
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
	     m_index2 = genp.daughter(bcc->genparticles,1)->index();
	     m_WMinusdecay2 = bcc->genparticles->at(m_index2);
	   }
       }
   }
  
  // W not linked correctly -> W decay products have top as mother
  if(m_pdgId1==0 || m_pdgId2==0 ){
    
    //search all particles with top as mother; store that ones which are not W or b (or equivalent light quark)
    for(unsigned int i=0; i<bcc->genparticles->size(); ++i)
      {
	GenParticle genp = bcc->genparticles->at(i);
	
	bool notfilled1=true;
	bool notfilled2=true;

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
	}

      }
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
