#include "HypothesisDiscriminator.h"

HypothesisDiscriminator::HypothesisDiscriminator(std::string label_name):  m_logger ( "HypothesisDiscriminator" ){

  m_label = label_name;
  m_filled = false;
}

void HypothesisDiscriminator::FillDiscriminatorValues(){
  ObjectHandler* objs = ObjectHandler::Instance();
  BaseCycleContainer* bcc = objs->GetBaseCycleContainer();
  if(!bcc->recoHyps || bcc->recoHyps->size()==0){
    m_logger << ERROR << "No reconstruction hypotheses in this event."<< SLogger::endmsg;
  }
  m_filled = bcc->recoHyps->at(0).has_discriminator(m_label);
}


ReconstructionHypothesis* HypothesisDiscriminator::GetHypWithSmallestDiscriminator(){

  ObjectHandler* objs = ObjectHandler::Instance();
  BaseCycleContainer* bcc = objs->GetBaseCycleContainer();

  float min=double_infinity();

  ReconstructionHypothesis* hyp = NULL;

  for(unsigned int i=0; i<bcc->recoHyps->size(); ++i){
    float discr=bcc->recoHyps->at(i).discriminator(m_label);
    if(discr<min){
      min = discr;
      hyp = &bcc->recoHyps->at(i);
    }
  }

  return hyp;

}

ReconstructionHypothesis* HypothesisDiscriminator::GetHypWithLargestDiscriminator(){

  ObjectHandler* objs = ObjectHandler::Instance();
  BaseCycleContainer* bcc = objs->GetBaseCycleContainer();

  float max=-1*double_infinity();

  ReconstructionHypothesis* hyp = NULL;

  for(unsigned int i=0; i<bcc->recoHyps->size(); ++i){
    float discr=bcc->recoHyps->at(i).discriminator(m_label);
    if(discr>max){
      max = discr;
      hyp = &bcc->recoHyps->at(i);
    }
  }

  return hyp;

}


void Chi2Discriminator::FillDiscriminatorValues(){
  HypothesisDiscriminator::FillDiscriminatorValues();
  if(m_filled) return;

  ObjectHandler* objs = ObjectHandler::Instance();
  BaseCycleContainer* bcc = objs->GetBaseCycleContainer();

  for(unsigned int i=0; i<bcc->recoHyps->size(); ++i){

    ReconstructionHypothesis* hyp =  &bcc->recoHyps->at(i);

    const double mass_thad = 179;
    const double mass_thad_sigma = 14;
    
    const double mass_tlep = 172;
    const double mass_tlep_sigma = 18;
    
//     const double drsum_thad = 0.7;
//     const double drsum_thad_sigma = 0.3;
    
//     const double drsum_tlep = 1.3;
//     const double drsum_tlep_sigma = 0.5;

//     double tlep_deltar = deltaR(hyp->toplep_v4(), bcc->jets->at(hyp->blep_index()).v4() ) + deltaR(hyp->toplep_v4(), hyp->neutrino_v4()) + deltaR(hyp->toplep_v4(), hyp->lepton().v4());
//     double thad_deltar = 0.0;

//     for(unsigned int j=0; j< hyp->tophad_jets_indices().size();++j){
//       thad_deltar += deltaR(hyp->tophad_v4(), bcc->jets->at(hyp->tophad_jets_indices().at(j)).v4());
//     }

    double mass_thad_rec=0;
    double mass_tlep_rec=0;

    if(hyp->toplep_v4().isTimelike()){
      mass_tlep_rec=hyp->toplep_v4().mass();
    }
    else{
      mass_tlep_rec= -sqrt(-hyp->toplep_v4().mass2());
    }
    if(hyp->tophad_v4().isTimelike()){
      mass_thad_rec=hyp->tophad_v4().mass();
    }
    else{
      mass_thad_rec= -sqrt(-hyp->tophad_v4().mass2());
    }


    double chi2_thad = pow((mass_thad_rec - mass_thad) / mass_thad_sigma, 2);
//     if( hyp->tophad_jets_indices().size()>1) 
//       chi2_thad += pow((thad_deltar - drsum_thad) / drsum_thad_sigma, 2);
    double chi2_tlep = pow((mass_tlep_rec - mass_tlep) / mass_tlep_sigma, 2); //+ pow((tlep_deltar - drsum_tlep) / drsum_tlep_sigma, 2);
    
    // make the chi2 on the leptonic side very high if there is more than 1 jet assigned:
//     if(hyp->toplep_jets_indices().size() > 1){
//       chi2_tlep += 1000.0;
//     }

    hyp->add_qualityflag(m_label, chi2_thad + chi2_tlep); //m_label="Chi2"
    hyp->add_qualityflag("Chi2_tlep", chi2_tlep);
    hyp->add_qualityflag("Chi2_thad", chi2_thad);
  }

}

ReconstructionHypothesis* Chi2Discriminator::GetBestHypothesis(){

  return GetHypWithSmallestDiscriminator();

}


void BestPossibleDiscriminator::FillDiscriminatorValues(){
  HypothesisDiscriminator::FillDiscriminatorValues();
  if(m_filled) return;

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  for(unsigned int i=0; i<bcc->recoHyps->size(); ++i){

    ReconstructionHypothesis* hyp =  &bcc->recoHyps->at(i);

    float bp_deltar=0;

    if(calc->GetGenParticles ())
      bp_deltar = deltaR(calc->GetTTbarGen()->Top().v4(),hyp->top_v4()) + deltaR(calc->GetTTbarGen()->Antitop().v4(),hyp->antitop_v4());
    
    hyp->add_qualityflag(m_label, bp_deltar);

  }

}

ReconstructionHypothesis* BestPossibleDiscriminator::GetBestHypothesis(){

  return GetHypWithSmallestDiscriminator();

}
