#ifndef TTbarGen_H
#define TTbarGen_H

// Local include(s):
#include "include/ObjectHandler.h"

/**
 *   Class for ttbar generator truth information
 */

class TTbarGen
{
  
 public:
  /// Default constructor, loops over generator particle list and fills the relevant particles
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
