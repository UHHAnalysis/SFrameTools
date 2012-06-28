#ifndef TTbarGen_H
#define TTbarGen_H

// ROOT include(s):
#include <TObject.h>
#include <TString.h>
#include <iostream>

// Local include(s):
#include "include/ObjectHandler.h"

/**
 *   Example class for booking and filling histograms
 *
 *   This class books and fills a collection of histograms.
 *   It should have a unique name, such that the histograms
 *   of multiple instances of this class are ordered in the
 *   output file. 
 *   Always sort your histograms and used methods topically.
 *   This example collection can be used for data and reconstructed
 *   MC events.
 *   
 *   @version $Revision: 1.1 $
 */

class TTbarGen
{
  
 public:
  /// Named constructor
  TTbarGen();
  
  /// Default destructor
  ~TTbarGen();

  GenParticle Top(); 
  GenParticle Antitop();
  GenParticle WTop(); 
  GenParticle WAntitop(); 
  GenParticle bTop(); 
  GenParticle bAntitop(); 
  
 private:

  GenParticle m_Top; 
  GenParticle m_Antitop;
  GenParticle m_WTop; 
  GenParticle m_WAntitop;
  GenParticle m_bTop; 
  GenParticle m_bAntitop;
  
}; // class TTbarGen


#endif // TTbarGen_H
