// Dear emacs, this is -*- c++ -*-
#ifndef Selection_H
#define Selection_H

#include <TString.h>

#include "SFrameTools/include/Objects.h"
#include "SFrameTools/include/BaseCycleContainer.h"
#include "SFrameTools/include/EventCalc.h"
#include "core/include/SLogger.h"

/**
 *  @short basic class for selection modules
 *  @author Thomas Peiffer
 * 
 * New classes chould implement pass(EventCalc &), the pass(BaseCycleContainer*) is only kept
 * for backward compatibility.
 */
class SelectionModule{
public:
   virtual ~SelectionModule(){}
  
   // return whether or not the event passes the selection
   virtual bool pass(BaseCycleContainer*) {
       assert(false);
   }
   
   // please override this for new classes. The default implementation
   // provided here is for backwards-compatibility and calls the BaseCycleContainer* variant
   virtual bool pass(EventCalc & event);
   
   virtual std::string description() = 0;
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
  /// Default constructor
  Selection(const std::string & name);
  
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
   * same as passSelection(BaseCycleContainer *bcc) but takes BaseCycleContainer from EventCalc
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
