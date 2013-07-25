#include "include/HypothesisDiscriminator.h"

HypothesisDiscriminator::HypothesisDiscriminator(std::string label_name):  m_logger ( "HypothesisDiscriminator" ){

  m_label = label_name;
  m_filled = false;
}

void HypothesisDiscriminator::FillDiscriminatorValues(){
  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();
  if(!bcc->recoHyps || bcc->recoHyps->size()==0){
    m_logger << ERROR << "No reconstruction hypotheses in this event."<< SLogger::endmsg;
  }
  m_filled = bcc->recoHyps->at(0).has_discriminator(m_label);
}


ReconstructionHypothesis* HypothesisDiscriminator::GetHypWithSmallestDiscriminator(){

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

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

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

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

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  for(unsigned int i=0; i<bcc->recoHyps->size(); ++i){

    ReconstructionHypothesis* hyp =  &bcc->recoHyps->at(i);

    const double mass_thad = 181;
    const double mass_thad_sigma = 15;
    
    const double mass_tlep = 174;
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


void SumDeltaRDiscriminator::FillDiscriminatorValues(){
  HypothesisDiscriminator::FillDiscriminatorValues();
  if(m_filled) return;

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  for(unsigned int i=0; i<bcc->recoHyps->size(); ++i){

    ReconstructionHypothesis* hyp =  &bcc->recoHyps->at(i);

    
    double tlep_deltar = deltaR(hyp->toplep_v4(), bcc->jets->at(hyp->blep_index()).v4() ) + deltaR(hyp->toplep_v4(), hyp->neutrino_v4()) + deltaR(hyp->toplep_v4(), hyp->lepton().v4());
//     double thad_deltar = 0.0;

//      for(unsigned int j=0; j< hyp->tophad_jets_indices().size();++j){
//        thad_deltar += deltaR(hyp->tophad_v4(), bcc->jets->at(hyp->tophad_jets_indices().at(j)).v4());
//      }

    double tlep_thad_deltar = deltaR(hyp->toplep_v4(),hyp->toplep_v4());

    hyp->add_qualityflag(m_label, tlep_deltar - 0.0001*tlep_thad_deltar); 

  }

}

ReconstructionHypothesis* SumDeltaRDiscriminator::GetBestHypothesis(){

  return GetHypWithSmallestDiscriminator();

}

void CorrectMatchDiscriminator::FillDiscriminatorValues(){
  HypothesisDiscriminator::FillDiscriminatorValues();
  if(m_filled) return;

  EventCalc* calc = EventCalc::Instance();
  BaseCycleContainer* bcc = calc->GetBaseCycleContainer();

  double bestdr = 99999.;

  for(unsigned int i=0; i<bcc->recoHyps->size(); ++i){

    ReconstructionHypothesis* hyp =  &bcc->recoHyps->at(i);

    float correct_dr=9999;

    if(calc->GetGenParticles ()){
      //consider l+jets events only
      if( calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_ehad || calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_muhad){

	//list of jets at hadronic top
	std::vector<Jet> hadr_jets;
	for(unsigned int i=0; i<hyp->tophad_jets_indices().size(); ++i){
	  unsigned int index = hyp->tophad_jets_indices().at(i);
	  hadr_jets.push_back( bcc->jets->at( index ) );
	}
	//list of jets at leptonic top
	std::vector<Jet> lept_jets;
	for(unsigned int i=0; i<hyp->toplep_jets_indices().size(); ++i){
	  unsigned int index = hyp->toplep_jets_indices().at(i);
	  lept_jets.push_back( bcc->jets->at( index ) );
	}

	//index lists of jets that can be matched to partons
	std::vector<int> matched_hadr_jets;
	std::vector<int> matched_lept_jets;

	//add b jets
	//top is decaying leptonically and anti-top hadronically
	int index_l, index_h;
	if( abs(calc->GetTTbarGen()->Wdecay1().pdgId()) >=11 && abs(calc->GetTTbarGen()->Wdecay1().pdgId())<=16){
	  correct_dr = isMatched(calc->GetTTbarGen()->bTop(),lept_jets, index_l) + isMatched(calc->GetTTbarGen()->bAntitop(),hadr_jets, index_h);
	  matched_lept_jets.push_back(index_l);
	  matched_hadr_jets.push_back(index_h);  
	}
	//vice versa
	else if (  abs(calc->GetTTbarGen()->WMinusdecay1().pdgId()) >=11 && abs(calc->GetTTbarGen()->WMinusdecay1().pdgId())<=16 ) {
	  correct_dr = isMatched(calc->GetTTbarGen()->bTop(),hadr_jets, index_h) + isMatched(calc->GetTTbarGen()->bAntitop(),lept_jets, index_l); 
	  matched_lept_jets.push_back(index_l);
	  matched_hadr_jets.push_back(index_h);  
	}
	else{
	  std::cerr << "This should not happen" <<std::endl;
	}


	//add quarks from W decays
	int n_hadr_jets=0;
	if(abs(calc->GetTTbarGen()->Wdecay1().pdgId())<6) {
	  n_hadr_jets++;
	  correct_dr += isMatched( calc->GetTTbarGen()->Wdecay1(), hadr_jets, index_h);
	  matched_hadr_jets.push_back(index_h);
	}
	if(abs(calc->GetTTbarGen()->Wdecay2().pdgId())<6) {
	  n_hadr_jets++;
	  correct_dr += isMatched( calc->GetTTbarGen()->Wdecay2(), hadr_jets, index_h);
	  matched_hadr_jets.push_back(index_h);
	}
	if(abs(calc->GetTTbarGen()->WMinusdecay1().pdgId())<6) {
	  n_hadr_jets++;
	  correct_dr += isMatched( calc->GetTTbarGen()->WMinusdecay1(), hadr_jets, index_h);
	  matched_hadr_jets.push_back(index_h);
	}
	if(abs(calc->GetTTbarGen()->WMinusdecay2().pdgId())<6) {
	  n_hadr_jets++;
	  correct_dr += isMatched( calc->GetTTbarGen()->WMinusdecay2(), hadr_jets, index_h);
	  matched_hadr_jets.push_back(index_h);
	}
	if(n_hadr_jets!=2){
	  std::cerr << "This should also not happen" << n_hadr_jets << std::endl;
	}

	//check that every jet in hadr_jets and lept_jets got matched
	for(int j=0; j<(int)hadr_jets.size(); ++j){
	  bool found=false;
	  for(unsigned int k=0; k<matched_hadr_jets.size(); ++k){
	    if(j==matched_hadr_jets[k]) {
	      found=true; 
	      break;
	    }
	  }
	  if(!found) correct_dr =9999.;
	}
	for(int j=0; j<(int)lept_jets.size(); ++j){
	  bool found=false;
	  for(unsigned int k=0; k<matched_lept_jets.size(); ++k){
	    if(j==matched_lept_jets[k]) {
	      found=true; 
	      break;
	    }
	  }
	  if(!found) correct_dr =9999.;
	}

	//add deltaR between reconstructed and true neutrino
	GenParticle neutrino = calc->GetTTbarGen()->Neutrino();
	correct_dr += deltaR(neutrino.v4(), hyp->neutrino_v4());
	if(correct_dr < bestdr) bestdr = correct_dr;
      }
    }
         
    hyp->add_qualityflag(m_label, correct_dr);

  }

//   int dummy[6]={-1,-1,-1,-1,-1,-1};
//   float alljets_dr = isMatched(calc->GetTTbarGen()->bTop(), *bcc->jets,dummy[0] ) + isMatched(calc->GetTTbarGen()->bAntitop(), *bcc->jets,dummy[1]); 
//   if(abs(calc->GetTTbarGen()->Wdecay1().pdgId())<6) {
//     alljets_dr += isMatched( calc->GetTTbarGen()->Wdecay1(), *bcc->jets,dummy[2]);
//   }
//   if(abs(calc->GetTTbarGen()->Wdecay2().pdgId())<6) {
//     alljets_dr += isMatched( calc->GetTTbarGen()->Wdecay2(), *bcc->jets,dummy[3]);
//   }
//   if(abs(calc->GetTTbarGen()->WMinusdecay1().pdgId())<6) {
//     alljets_dr += isMatched( calc->GetTTbarGen()->WMinusdecay1(), *bcc->jets,dummy[4]);
//   }
//   if(abs(calc->GetTTbarGen()->WMinusdecay2().pdgId())<6) {
//     alljets_dr += isMatched( calc->GetTTbarGen()->WMinusdecay2(), *bcc->jets,dummy[5]);
//   }

//   if(alljets_dr<9999 && bestdr>=9999){
//     if(calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_ehad || calc->GetTTbarGen()->DecayChannel() == TTbarGen::e_muhad){
//       std::cout << "sdiohsgdihsgdkihsgki: " << alljets_dr << "   " << bestdr <<std::endl;
//       std::cout << dummy[0] << "  " << dummy[1] << "  " << dummy[2] << "  " << dummy[3] << "  " << dummy[4] << "  " << dummy[5] << std::endl;
//     }
//   }
}

ReconstructionHypothesis* CorrectMatchDiscriminator::GetBestHypothesis(){

  ReconstructionHypothesis *hyp = GetHypWithSmallestDiscriminator();
  //return null pointer for non-matchable events
  if(hyp->discriminator(m_label)>=9999) return NULL;
  return hyp;

}

float CorrectMatchDiscriminator::isMatched(GenParticle p, std::vector<Jet> jets, int& index){

  float mindr = 9999.;
  index = -1;
  for(unsigned int i=0; i<jets.size(); ++i){
    float dR = deltaR( p.v4(), jets.at(i).v4());
    if( dR <0.3 && dR<mindr) {
      mindr=dR;
      index=i;
    }
  }
  return mindr;

}
