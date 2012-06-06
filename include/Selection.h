// Dear emacs, this is -*- c++ -*-
#ifndef Selection_H
#define Selection_H

#include <TString.h>
#include "Objects.h"
#include "BaseCycleContainer.h"
#include "ObjectHandler.h"
#include "../../core/include/SLogger.h"

class SelectionModule{
 public:
  SelectionModule(){};
  virtual ~SelectionModule(){};
   
  virtual bool pass(BaseCycleContainer*)=0;
  virtual std::string description()=0;
};


class Selection{
 public:
  
  //Selection(){};
  Selection(std::string name);
  ~Selection(){};
  
  void addSelectionModule(SelectionModule*);
  void clearSelectionModulesList();
  bool passSelection(BaseCycleContainer *bcc);
  bool passInvertedSelection(BaseCycleContainer *bcc);
  bool passSelection();
  bool passInvertedSelection();

  void printCutFlow();

  TString GetName() {return m_name;}
  void SetName(TString name) { m_name = name;}

 private:
  TString m_name;
  mutable SLogger m_logger;
  std::vector<SelectionModule*> m_cuts;
  std::vector<int> m_cutflow;
  int Ntotal;

};



#endif // Selection_H
