#include <iomanip>

#include "../include/Selection.h"

Selection::Selection(std::string name):
  m_logger ( name.c_str() ){

  m_name = name.c_str();
  clearSelectionModulesList();
  Ntotal=0;
  m_isactive = true;
}

void Selection::resetCutFlow(){
  Ntotal=0;
  for(unsigned int i=0; i<m_cuts.size(); ++i){
    m_cutflow[i]=0;
  }
}

void Selection::addSelectionModule(SelectionModule* sel){
  m_logger << DEBUG << "Adding selection module: " << sel->description() << SLogger::endmsg;
  m_cuts.push_back(sel);
  m_cutflow.push_back(0);
}

void Selection::clearSelectionModulesList(){
  m_cuts.clear();
  m_cutflow.clear();
  Ntotal=0;
}

bool Selection::passSelection(BaseCycleContainer *bcc){

  if (!m_isactive) return true; // always true if the selection is not active

  Ntotal++;
  if(m_cuts.size()!=m_cutflow.size()){
    m_logger << WARNING << "size of cut list != number of entries in cut flow table "<< SLogger::endmsg;
  }
  for(unsigned int i=0; i<m_cuts.size(); ++i){
    if(!m_cuts[i]->pass(bcc)) return false;
    m_cutflow[i]++;
  }
  return true;
}

bool Selection::passInvertedSelection(BaseCycleContainer *bcc){

  if (!m_isactive) return true; // always true if the selection is not active
  
  for(unsigned int i=0; i<m_cuts.size(); ++i){
    if(!m_cuts[i]->pass(bcc)) return true;
  }
  return false;
}

bool Selection::passSelection(){
  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  return passSelection(bcc);
}

bool Selection::passInvertedSelection(){
  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  return passInvertedSelection(bcc);
}

void Selection::printCutFlow(){

  using namespace std;

  m_logger << INFO << "--------------- Cut Flow Table of Selection " << m_name << " ---------------"<< SLogger::endmsg;
  if (!m_isactive){
    m_logger << INFO <<  "Selection was not active." << SLogger::endmsg;
    m_logger << INFO << "-----------------------------------------------------------------------------"<< SLogger::endmsg;
    return;
  }

  if(m_cuts.size()!=m_cutflow.size()){
    m_logger << WARNING << "size of cut list != number of entries in cut flow table "<< SLogger::endmsg;
  }
  else{
      m_logger << INFO << setw(12) << "Events" <<  " |  Description" << SLogger::endmsg;
      m_logger << INFO << "-------------+-----------------------------------------------------"<< SLogger::endmsg;
      m_logger << INFO << setw(12) << Ntotal << " | Events entered the selection. " << SLogger::endmsg;
    for(unsigned int i=0; i<m_cuts.size(); ++i){
      m_logger << INFO << setw(12) << m_cutflow[i] << " | left after: " << m_cuts[i]->description() << SLogger::endmsg;
    }
  }
  m_logger << INFO << "-------------+---------------------------------------------------------------"<< SLogger::endmsg;

}
