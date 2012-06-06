#include <iomanip>

#include "../include/Selection.h"

Selection::Selection(std::string name):
  m_logger ( name.c_str() ){

  m_name = name.c_str();
  clearSelectionModulesList();
  Ntotal=0;
}

void Selection::addSelectionModule(SelectionModule* sel){
  m_cuts.push_back(sel);
  m_cutflow.push_back(0);
}

void Selection::clearSelectionModulesList(){
  m_cuts.clear();
  m_cutflow.clear();
  Ntotal=0;
}

bool Selection::passSelection(BaseCycleContainer *bcc){
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
  
  for(unsigned int i=0; i<m_cuts.size(); ++i){
    if(!m_cuts[i]->pass(bcc)) return true;
  }
  return false;
}

bool Selection::passSelection(){
  ObjectHandler* objs = ObjectHandler::Instance();
  BaseCycleContainer* bcc = objs->GetBaseCycleContainer();
  return passSelection(bcc);
}

bool Selection::passInvertedSelection(){
  ObjectHandler* objs = ObjectHandler::Instance();
  BaseCycleContainer* bcc = objs->GetBaseCycleContainer();
  return passInvertedSelection(bcc);
}

void Selection::printCutFlow(){

  using namespace std;

  m_logger << INFO << "-------------------------- Cut Flow Table -------------------------"<< SLogger::endmsg;
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
  m_logger << INFO << "-------------+-----------------------------------------------------"<< SLogger::endmsg;

}
