// Dear emacs, this is -*- c++ -*-
#ifndef JetProps_H
#define JetProps_H

#include "Objects.h"
#include "SLogger.h"
#include "Utils.h"

#include "Njettiness.h"


/**
 *  @class for the calculation of jet properties
 *   most of the functions here use FastJet to calculate
 *   observables like N-subjettiness, mass-drop,
 *   sub-jets or jet grooming variables
 *   
 *   If you use this functionality, make sure that the 
 *   jet constituents are provided
 *
 *
 *  @author Roman Kogler
 */

class JetProps
{
 public:
  JetProps();
  JetProps(TopJet* jet);
  ~JetProps();

  double GetNsubjettiness(int N, Njettiness::AxesMode mode, double beta, double R0, double Rcutoff=std::numeric_limits<double>::max());

  double GetPrunedNsubjettiness(int N, Njettiness::AxesMode mode, double beta, double R0, double Rcutoff=std::numeric_limits<double>::max());

  std::vector<fastjet::PseudoJet> GetFastJet(double R0, fastjet::JetAlgorithm jet_algo=fastjet::cambridge_algorithm);
  fastjet::PseudoJet GetPrunedJet(fastjet::PseudoJet injet);

  double GetQjetVolatility(int seed, double R0);

  double FindMean( std::vector< double > qjetmasses );
  double FindRMS( std::vector< double > qjetmasses );

  std::vector<fastjet::PseudoJet> GetJetConstituents();

 private:

  TopJet* m_jet;
  SLogger m_logger;
  fastjet::ClusterSequence* m_JetFinder;
  fastjet::JetDefinition* m_JetDef ;

};


#endif // JetProps_H
