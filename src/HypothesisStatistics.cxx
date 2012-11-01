#include "include/HypothesisStatistics.h"

HypothesisStatistics::HypothesisStatistics(std::string name):
  m_logger ( name.c_str() ){

  Reset();

}

void HypothesisStatistics::Reset(){
  m_ntotal=0;
  m_identical=0;
  m_toplep_identical=0;
  m_tophad_identical=0;
}


void HypothesisStatistics::FillHyps(ReconstructionHypothesis *hyp1, ReconstructionHypothesis *hyp2){

  m_ntotal++;
  if(hyp1==hyp2) m_identical++;
  if(hyp1->toplep_v4() == hyp2->toplep_v4()) m_toplep_identical++;
  if(hyp1->tophad_v4() == hyp2->tophad_v4()) m_tophad_identical++; 
  
}

void HypothesisStatistics::PrintStatistics(){
  m_logger << INFO << "-------------------------- Reconstruction Hypothesis Statistics -------------------------"<< SLogger::endmsg;
  m_logger << INFO << "number of analysed events:   " << m_ntotal << SLogger::endmsg;
  m_logger << INFO << "hypotheses are identical:    " << 100.*(double)m_identical/((double)m_ntotal)<< "%" << SLogger::endmsg;
  m_logger << INFO << "leptonic tops are identical: " << 100.*(double)m_toplep_identical/((double)m_ntotal)<< "%" << SLogger::endmsg;
  m_logger << INFO << "hadronic tops are identical: " << 100.*(double)m_tophad_identical/((double)m_ntotal)<< "%" << SLogger::endmsg;

}
