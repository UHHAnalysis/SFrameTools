#include "include/TTbarGen.h"

using namespace std;

TTbarGen::TTbarGen()
{
  ObjectHandler* objs = ObjectHandler::Instance();
  BaseCycleContainer* bcc = objs->GetBaseCycleContainer();
  
  for(unsigned int i=0; i<bcc->genparticles->size(); ++i)
   {
     GenParticle genp = bcc->genparticles->at(i);
     if ( genp.pdgId() == 6 ) 
       {
	 m_Top = genp;
	 if (genp.daughter(bcc->genparticles,1) && genp.daughter(bcc->genparticles,2))
	   {
	     int indexW;
	     int indexb;
	     if (genp.daughter(bcc->genparticles,1)->pdgId() == 24) 
	       {
		 indexW = genp.daughter(bcc->genparticles,1)->index();
		 m_WTop = bcc->genparticles->at(indexW);
	       }
	     if (genp.daughter(bcc->genparticles,1)->pdgId() == 5)  
	       {
		 indexb = genp.daughter(bcc->genparticles,1)->index();
		 m_bTop = bcc->genparticles->at(indexb);
	       }
	     if (genp.daughter(bcc->genparticles,2)->pdgId() == 24)
	       {
		 indexW = genp.daughter(bcc->genparticles,2)->index();
		 m_WTop = bcc->genparticles->at(indexW);
	       }
	     if (genp.daughter(bcc->genparticles,2)->pdgId() == 5)  
	       {
		 indexb = genp.daughter(bcc->genparticles,2)->index();
		 m_bTop = bcc->genparticles->at(indexb);
	       }
	   }
       }
     if ( genp.pdgId() == -6 ) 
       {
	 m_Antitop = genp;
	 if (genp.daughter(bcc->genparticles,1) && (genp.daughter(bcc->genparticles,2)))
	   {
	     int indexW;
	     int indexb;
	     if (genp.daughter(bcc->genparticles,1)->pdgId() == -24) 
	       {
		 indexW = genp.daughter(bcc->genparticles,1)->index();
		 m_WAntitop = bcc->genparticles->at(indexW);
	       }
	     if (genp.daughter(bcc->genparticles,1)->pdgId() == -5)  
	       {
		 indexb = genp.daughter(bcc->genparticles,1)->index();
		 m_bAntitop = bcc->genparticles->at(indexb);
	       }
	     if (genp.daughter(bcc->genparticles,2)->pdgId() == -24)
	       {
		 indexW = genp.daughter(bcc->genparticles,2)->index();
		 m_WAntitop = bcc->genparticles->at(indexW);
	       }
	     if (genp.daughter(bcc->genparticles,2)->pdgId() == -5)  
	       {
		 indexb = genp.daughter(bcc->genparticles,2)->index();
		 m_bAntitop = bcc->genparticles->at(indexb);
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



