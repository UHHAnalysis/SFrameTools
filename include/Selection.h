// Dear emacs, this is -*- c++ -*-
#ifndef Selection_H
#define Selection_H

#include <TString.h>
#include "Objects.h"
#include "BaseCycleContainer.h"
#include "ObjectHandler.h"
#include "../../core/include/SLogger.h"

/**
 *  @short basic class for selection modules
 *  @author Thomas Peiffer
 */

class SelectionModule{
 public:
  //Default constructor
  SelectionModule(){};
  //Default destructor
  virtual ~SelectionModule(){};
   
  virtual bool pass(BaseCycleContainer*)=0;
  virtual std::string description()=0;
};

/**
 *  @short basic class for a selection
 *
 *         contains a list of SelectionModules
 *         and produces cut flow tables
 *
 *  @author Thomas Peiffer
 */

class Selection{
 public:
  
  //Selection(){};
  /// Default constructor
  Selection(std::string name);
  /// Default destructor
  ~Selection(){};
  
  /// add a new SelectionModule to the event selection
  void addSelectionModule(SelectionModule*);
  /// remove all booked SelectionModules
  void clearSelectionModulesList();
  /// set all entries in the cut flow table to 0
  void resetCutFlow();

  /**
   * loop over all booked SelectionModules and call the pass() routine of each module for a given BaseCycleContainer
   * @see passSelection()
   * @see SelectionModule::pass(BaseCycleContainer*)
   */
  bool passSelection(BaseCycleContainer *bcc);
  bool passInvertedSelection(BaseCycleContainer *bcc);
  /**
   * same as passSelection(BaseCycleContainer *bcc) but takes BaseCycleContainer from ObjectHandler
   * @see passSelection(BaseCycleContainer *bcc)
   */
  bool passSelection();
  bool passInvertedSelection();

  /**
   * enable or disable the selection, the selection is enabled by default
   */
  void DisableSelection(){m_isactive = false;}
  void EnableSelection(){m_isactive = true;}
  bool IsActive(){return m_isactive;}

  /// print the cut-flow table for the booked SelectionModules, to be called in the EndInputData routine of a cycle
  void printCutFlow();

  TString GetName() {return m_name;}
  void SetName(TString name) { m_name = name;}

 private:
  TString m_name;
  mutable SLogger m_logger;
  std::vector<SelectionModule*> m_cuts;
  std::vector<int> m_cutflow;
  int Ntotal;
  bool m_isactive;

};



#endif // Selection_H
