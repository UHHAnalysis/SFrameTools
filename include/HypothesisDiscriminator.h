#ifndef HypothesisDiscriminator_H
#define HypothesisDiscriminator_H

#include "BaseCycleContainer.h"
#include "ObjectHandler.h"
#include "Utils.h"
#include "ReconstructionHypothesis.h"
//#include "../../core/include/SLogger.h"
#include "EventCalc.h"

/**
 *  @short basic class to select a reconstruction hypothesis with a certain discriminator value
 *
 * Each class that inherits from the HypothesisDiscriminator class provides fucntions to fill disriminator values for every ReconstructionHypothesis in the actual BaseCycleContainer
 * and to select the hypothesis with the best discriminator value.
 *
 */

class HypothesisDiscriminator{
 public:
  ///Default constructor
  HypothesisDiscriminator(std::string label_name);
  ///Default destructor
  virtual ~HypothesisDiscriminator(){};

  /// return the ReconstructionHypothesis with the best discriminator value (discriminator values have to be filled with FillDiscriminatorValues before calling this routine).
  virtual ReconstructionHypothesis* GetBestHypothesis()=0;

  /// determine the discriminator values and store them in the list of reconstruction hypotheses in the BaseCycleContainer
  virtual void FillDiscriminatorValues();

  std::string GetLabel(){return m_label;}

 protected:

  ReconstructionHypothesis* GetHypWithSmallestDiscriminator();
  ReconstructionHypothesis* GetHypWithLargestDiscriminator();

  mutable SLogger m_logger;

  string m_label;
  bool m_filled;

};

/**
 * BestPossible hypothesis discriminator, only works on MC
 */
class BestPossibleDiscriminator: public HypothesisDiscriminator{
 public:
  BestPossibleDiscriminator():HypothesisDiscriminator("BestPossible"){};
  ~BestPossibleDiscriminator(){};

};

/**
 * hypothesis discriminator Chi^2 
 */
class Chi2Discriminator: public HypothesisDiscriminator{
 public:
  Chi2Discriminator():HypothesisDiscriminator("Chi2"){};
  ~Chi2Discriminator(){};

  virtual void FillDiscriminatorValues();
  
  virtual ReconstructionHypothesis* GetBestHypothesis();

  
};

#endif
