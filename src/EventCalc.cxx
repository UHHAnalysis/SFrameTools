// Dear emacs, this is -*- c++ -*-

#include "include/EventCalc.h"

#include <iostream>

EventCalc* EventCalc::m_instance = NULL;

EventCalc* EventCalc::Instance()
{
  // Get a pointer to the object handler.
  // This is the only way to access this class, 
  // since it's a singleton. This method is accessible
  // from everywhere.

  if (m_instance == NULL){
    m_instance = new EventCalc();
  }
  return m_instance;    
}

EventCalc::EventCalc() : m_logger( "EventCalc" )
{
  // constructor: initialise all variables
  m_logger << DEBUG << "Constructor called." << SLogger::endmsg;
  m_bcc = NULL;
  m_lumi = NULL;
}

EventCalc::~EventCalc()
{
  // default destructor
}

void EventCalc::Reset()
{
  // reset: set all booleans to false
  // this has to be done at the beginning of each event
  // after a call to reset all quantities will be re-calculated 
  // when they are accessed
  
  // (re-)set the pointers using the ObjectHandler
  ObjectHandler* objs = ObjectHandler::Instance();
  m_bcc = objs->GetBaseCycleContainer();
  m_lumi = objs->GetLumiHandler();

  // reset booleans
  b_HT = false;
  b_HTlep = false;

}

BaseCycleContainer* EventCalc::GetBaseCycleContainer()
{
  // return the pointer to the container with all objects
  if (!m_bcc){
    m_logger << WARNING << "Pointer to BaseCycleContainer is NULL." << SLogger::endmsg;
  }
  return m_bcc;
}

LuminosityHandler* EventCalc::GetLumiHandler()
{
  // return the pointer to the container with all objects
  if (!m_lumi){
    m_logger << WARNING << "Pointer to LumiHandler is NULL." << SLogger::endmsg;
  }
  return m_lumi;
}

double EventCalc::GetHT()
{
  // calculate HT, which is defined as the scalar sum of all
  // jets, leptons and missing transverse momentum in the event
  if (!b_HT){

    b_HT = true;
    m_HT = 0;

    // add lepton pt and MET 
    m_HT += GetHTlep();

    // sum over pt of all jets
    if(m_bcc->jets){
      for(unsigned int i=0; i<m_bcc->jets->size(); ++i){
	m_HT += m_bcc->jets->at(i).pt();
      }
    }

  }
  return m_HT;
}

double EventCalc::GetHTlep()
{
  // calculate HT_lep, which is defined as the scalar sum of all
  // leptons and missing transverse momentum in the event
  if (!b_HTlep){

    b_HTlep = true;
    m_HTlep=0;

    // sum over pt of all electrons
    if(m_bcc->electrons){
      for(unsigned int i=0; i<m_bcc->electrons->size(); ++i){
	m_HTlep += m_bcc->electrons->at(i).pt();
      }
    }

    // sum over pt of all muons
    if(m_bcc->muons){
      for(unsigned int i=0; i<m_bcc->muons->size(); ++i){
	m_HTlep += m_bcc->muons->at(i).pt();
      }
    }

    // sum over pt of all taus
    if(m_bcc->taus){
      for(unsigned int i=0; i<m_bcc->taus->size(); ++i){
	m_HTlep += m_bcc->taus->at(i).pt();
      }
    }
    
    // add MET
    if(m_bcc->met) m_HTlep += m_bcc->met->pt();   

  }
  return m_HTlep;
}

void EventCalc::PrintEventContent(){

  m_logger << INFO << "----------------- event content -----------------" << SLogger::endmsg;
  m_logger << INFO << "run: " << m_bcc->run <<  "   lumi block:" << m_bcc->luminosityBlock << "   event: " << m_bcc->event << SLogger::endmsg;
  m_logger << INFO << "MET = " << m_bcc->met->pt() << "   METphi = " << m_bcc->met->phi() <<  "   HTlep = " << GetHTlep() << SLogger::endmsg;
  if(m_bcc->electrons){m_logger << INFO << "Electrons:" << SLogger::endmsg;
    for(unsigned int i=0; i<m_bcc->electrons->size(); ++i){
      m_logger << INFO << "     " << i+1 << " pt = " << m_bcc->electrons->at(i).pt() <<"   eta = " << m_bcc->electrons->at(i).eta() << SLogger::endmsg;
    }
  }
  if(m_bcc->muons){m_logger << INFO << "Muons:" << SLogger::endmsg;
    for(unsigned int i=0; i<m_bcc->muons->size(); ++i){
      m_logger << INFO << "     " << i+1 << " pt = " << m_bcc->muons->at(i).pt() <<"   eta = " << m_bcc->muons->at(i).eta() << SLogger::endmsg;
    }
  }
  if(m_bcc->taus){m_logger << INFO << "Taus:" << SLogger::endmsg;
    for(unsigned int i=0; i<m_bcc->taus->size(); ++i){
      m_logger << INFO << "     " << i+1 << " pt = " << m_bcc->taus->at(i).pt() <<"   eta = " << m_bcc->taus->at(i).eta() << SLogger::endmsg;
    }
  }
  if(m_bcc->jets){m_logger << INFO << "Jets:" << SLogger::endmsg;
    for(unsigned int i=0; i<m_bcc->jets->size(); ++i){
      m_logger << INFO << "     " << i+1 << " pt = " << m_bcc->jets->at(i).pt() <<"   eta = " << m_bcc->jets->at(i).eta() << "   genjet_pt = " << m_bcc->jets->at(i).genjet_pt() << "   genjet_eta = " <<   m_bcc->jets->at(i).genjet_eta() <<SLogger::endmsg;
    }
  }
  if(m_bcc->topjets){m_logger << INFO << "TopJets:" << SLogger::endmsg;
    for(unsigned int i=0; i<m_bcc->topjets->size(); ++i){
      m_logger << INFO << "     " << i+1 << " pt = " << m_bcc->topjets->at(i).pt() <<"   eta = " << m_bcc->topjets->at(i).eta() << SLogger::endmsg;
    }
  }
  if(m_bcc->photons){m_logger << INFO << "Photons:" << SLogger::endmsg;
    for(unsigned int i=0; i<m_bcc->photons->size(); ++i){
      m_logger << INFO << "     " << i+1 << " pt = " << m_bcc->photons->at(i).pt() <<"   eta = " << m_bcc->photons->at(i).eta() << SLogger::endmsg;
    }
  }
}
