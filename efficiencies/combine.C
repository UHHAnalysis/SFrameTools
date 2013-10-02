void combine()
{
  // short macro to combine muon SFs from different files into one
  // all SFs are written out for 'tight' IDs!

  // target file
  TFile* tf = new TFile("muon_effs_2012_53x.root", "RECREATE");
  
  TDirectory* d_a = tf->mkdir("2012ABCD");
//   TDirectory* d_b = tf->mkdir("2012B");
//   TDirectory* d_c = tf->mkdir("2012C");
//   TDirectory* d_d = tf->mkdir("2012D");
  
  // eta binning: 
  // 0.0 - 0.9: eta bin 0
  // 0.9 - 1.2: eta bin 1
  // 1.2 - 2.1: eta bin 2
  
  // -------------------------------------------------------
  // ---------------- process 2012A ------------------------
  // -------------------------------------------------------

  TFile* file_ab = new TFile("MuonEfficiencies_RunABCD_53X.root", "READ");
  
  // tight ID
  TGraphAsymmErrors* sf_eta0_id = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_Tight_pt_abseta<0.9");
  TGraphAsymmErrors* sf_eta1_id = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_Tight_pt_abseta0.9-1.2");
  TGraphAsymmErrors* sf_eta2_id = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_Tight_pt_abseta1.2-2.1");

  sf_eta0_id->SetName("SF_eta0_ID_tight");
  sf_eta1_id->SetName("SF_eta1_ID_tight");
  sf_eta2_id->SetName("SF_eta2_ID_tight");

  // isolation
  TGraphAsymmErrors* sf_eta0_iso = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta<0.9");
  TGraphAsymmErrors* sf_eta1_iso = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta0.9-1.2");
  TGraphAsymmErrors* sf_eta2_iso = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta1.2-2.1");

  sf_eta0_iso->SetName("SF_eta0_ISO_tight");
  sf_eta1_iso->SetName("SF_eta1_ISO_tight");
  sf_eta2_iso->SetName("SF_eta2_ISO_tight");

  // trigger: IsoMu24
  TGraphAsymmErrors* sf_eta0_trig_isomu24 = (TGraphAsymmErrors*) file_ab->Get("IsoMu24_eta2p1_DATA_over_MC_TightID_IsodB_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD");
  TGraphAsymmErrors* sf_eta1_trig_isomu24 = (TGraphAsymmErrors*) file_ab->Get("IsoMu24_eta2p1_DATA_over_MC_TightID_IsodB_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD");
  TGraphAsymmErrors* sf_eta2_trig_isomu24 = (TGraphAsymmErrors*) file_ab->Get("IsoMu24_eta2p1_DATA_over_MC_TightID_IsodB_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD");

  sf_eta0_trig_isomu24->SetName("SF_eta0_TRIG_IsoMu24");
  sf_eta1_trig_isomu24->SetName("SF_eta1_TRIG_IsoMu24");
  sf_eta2_trig_isomu24->SetName("SF_eta2_TRIG_IsoMu24");

  // trigger: Mu40
  TGraphAsymmErrors* sf_eta0_trig_mu40 = (TGraphAsymmErrors*) file_ab->Get("Mu40_eta2p1_DATA_over_MC_TightID_PT_ABSETA_Barrel_0to0p9_pt45-500_2012ABCD");
  TGraphAsymmErrors* sf_eta1_trig_mu40 = (TGraphAsymmErrors*) file_ab->Get("Mu40_eta2p1_DATA_over_MC_TightID_PT_ABSETA_Transition_0p9to1p2_pt45-500_2012ABCD");
  TGraphAsymmErrors* sf_eta2_trig_mu40 = (TGraphAsymmErrors*) file_ab->Get("Mu40_eta2p1_DATA_over_MC_TightID_PT_ABSETA_Endcaps_1p2to2p1_pt45-500_2012ABCD");

  sf_eta0_trig_mu40->SetName("SF_eta0_TRIG_Mu40");
  sf_eta1_trig_mu40->SetName("SF_eta1_TRIG_Mu40");
  sf_eta2_trig_mu40->SetName("SF_eta2_TRIG_Mu40");

  d_a->cd();
  sf_eta0_id->Write();
  sf_eta1_id->Write();
  sf_eta2_id->Write();

  sf_eta0_iso->Write();
  sf_eta1_iso->Write();
  sf_eta2_iso->Write();

  sf_eta0_trig_isomu24->Write();
  sf_eta1_trig_isomu24->Write();
  sf_eta2_trig_isomu24->Write();

  sf_eta0_trig_mu40->Write();
  sf_eta1_trig_mu40->Write();
  sf_eta2_trig_mu40->Write();


  // -------------------------------------------------------
  // ---------------- process 2012B ------------------------
  // -------------------------------------------------------

  // TFile* file_ab = new TFile("MuonEfficiencies_Run_2012A_2012B_53X.root", "READ");
  
//   // tight ID
//   TGraphAsymmErrors* sf_eta0_id = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_Tight_pt_abseta<0.9");
//   TGraphAsymmErrors* sf_eta1_id = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_Tight_pt_abseta0.9-1.2");
//   TGraphAsymmErrors* sf_eta2_id = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_Tight_pt_abseta1.2-2.1");

//   sf_eta0_id->SetName("SF_eta0_ID_tight");
//   sf_eta1_id->SetName("SF_eta1_ID_tight");
//   sf_eta2_id->SetName("SF_eta2_ID_tight");

//   // isolation
//   TGraphAsymmErrors* sf_eta0_iso = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta<0.9");
//   TGraphAsymmErrors* sf_eta1_iso = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta0.9-1.2");
//   TGraphAsymmErrors* sf_eta2_iso = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta1.2-2.1");

//   sf_eta0_iso->SetName("SF_eta0_ISO_tight");
//   sf_eta1_iso->SetName("SF_eta1_ISO_tight");
//   sf_eta2_iso->SetName("SF_eta2_ISO_tight");

//   // trigger: IsoMu24
//   TGraphAsymmErrors* sf_eta0_trig_isomu24 = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_IsoMu24_eta2p1_TightIso_pt_abseta<0.9_2012B");
//   TGraphAsymmErrors* sf_eta1_trig_isomu24 = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_IsoMu24_eta2p1_TightIso_pt_abseta0.9-1.2_2012B");
//   TGraphAsymmErrors* sf_eta2_trig_isomu24 = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_IsoMu24_eta2p1_TightIso_pt_abseta1.2-2.1_2012B");

//   sf_eta0_trig_isomu24->SetName("SF_eta0_TRIG_IsoMu24");
//   sf_eta1_trig_isomu24->SetName("SF_eta1_TRIG_IsoMu24");
//   sf_eta2_trig_isomu24->SetName("SF_eta2_TRIG_IsoMu24");

//   // trigger: Mu40
//   TGraphAsymmErrors* sf_eta0_trig_mu40 = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_Mu40_eta2p1_Tight_pt_40_abseta<0.9_2012B");
//   TGraphAsymmErrors* sf_eta1_trig_mu40 = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_Mu40_eta2p1_Tight_pt_40_abseta0.9-1.2_2012B");
//   TGraphAsymmErrors* sf_eta2_trig_mu40 = (TGraphAsymmErrors*) file_ab->Get("DATA_over_MC_Mu40_eta2p1_Tight_pt_40_abseta1.2-2.1_2012B");

//   sf_eta0_trig_mu40->SetName("SF_eta0_TRIG_Mu40");
//   sf_eta1_trig_mu40->SetName("SF_eta1_TRIG_Mu40");
//   sf_eta2_trig_mu40->SetName("SF_eta2_TRIG_Mu40");

//   d_b->cd();
//   sf_eta0_id->Write();
//   sf_eta1_id->Write();
//   sf_eta2_id->Write();

//   sf_eta0_iso->Write();
//   sf_eta1_iso->Write();
//   sf_eta2_iso->Write();

//   sf_eta0_trig_isomu24->Write();
//   sf_eta1_trig_isomu24->Write();
//   sf_eta2_trig_isomu24->Write();

//   sf_eta0_trig_mu40->Write();
//   sf_eta1_trig_mu40->Write();
//   sf_eta2_trig_mu40->Write();


//   // -------------------------------------------------------
//   // ---------------- process 2012C ------------------------
//   // -------------------------------------------------------

//   TFile* file_c = new TFile("MuonEfficiencies_Run_2012C_53X.root", "READ");
  
//   // tight ID
//   TGraphAsymmErrors* sf_eta0_id = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_Tight_pt_abseta<0.9");
//   TGraphAsymmErrors* sf_eta1_id = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_Tight_pt_abseta0.9-1.2");
//   TGraphAsymmErrors* sf_eta2_id = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_Tight_pt_abseta1.2-2.1");

//   sf_eta0_id->SetName("SF_eta0_ID_tight");
//   sf_eta1_id->SetName("SF_eta1_ID_tight");
//   sf_eta2_id->SetName("SF_eta2_ID_tight");

//   // isolation
//   TGraphAsymmErrors* sf_eta0_iso = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta<0.9");
//   TGraphAsymmErrors* sf_eta1_iso = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta0.9-1.2");
//   TGraphAsymmErrors* sf_eta2_iso = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta1.2-2.1");

//   sf_eta0_iso->SetName("SF_eta0_ISO_tight");
//   sf_eta1_iso->SetName("SF_eta1_ISO_tight");
//   sf_eta2_iso->SetName("SF_eta2_ISO_tight");

//   // trigger: IsoMu24
//   TGraphAsymmErrors* sf_eta0_trig_isomu24 = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_IsoMu24_eta2p1_TightIso_pt_abseta<0.9");
//   TGraphAsymmErrors* sf_eta1_trig_isomu24 = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_IsoMu24_eta2p1_TightIso_pt_abseta0.9-1.2");
//   TGraphAsymmErrors* sf_eta2_trig_isomu24 = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_IsoMu24_eta2p1_TightIso_pt_abseta1.2-2.1");

//   sf_eta0_trig_isomu24->SetName("SF_eta0_TRIG_IsoMu24");
//   sf_eta1_trig_isomu24->SetName("SF_eta1_TRIG_IsoMu24");
//   sf_eta2_trig_isomu24->SetName("SF_eta2_TRIG_IsoMu24");

//   // trigger: Mu40
//   TGraphAsymmErrors* sf_eta0_trig_mu40 = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_Mu40_eta2p1_Tight_pt_40_abseta<0.9");
//   TGraphAsymmErrors* sf_eta1_trig_mu40 = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_Mu40_eta2p1_Tight_pt_40_abseta0.9-1.2");
//   TGraphAsymmErrors* sf_eta2_trig_mu40 = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_Mu40_eta2p1_Tight_pt_40_abseta1.2-2.1");

//   sf_eta0_trig_mu40->SetName("SF_eta0_TRIG_Mu40");
//   sf_eta1_trig_mu40->SetName("SF_eta1_TRIG_Mu40");
//   sf_eta2_trig_mu40->SetName("SF_eta2_TRIG_Mu40");

//   d_c->cd();
//   sf_eta0_id->Write();
//   sf_eta1_id->Write();
//   sf_eta2_id->Write();

//   sf_eta0_iso->Write();
//   sf_eta1_iso->Write();
//   sf_eta2_iso->Write();

//   sf_eta0_trig_isomu24->Write();
//   sf_eta1_trig_isomu24->Write();
//   sf_eta2_trig_isomu24->Write();

//   sf_eta0_trig_mu40->Write();
//   sf_eta1_trig_mu40->Write();
//   sf_eta2_trig_mu40->Write();


//   // -------------------------------------------------------
//   // ---------------- process 2012D ------------------------
//   // -------------------------------------------------------

//   TFile* file_d = new TFile("TriggerMuonEfficiencies_Run_2012D_53X.root", "READ");
//   TFile* file_c = new TFile("MuonEfficiencies_Run_2012C_53X.root", "READ");
  
//   // tight ID
//   TGraphAsymmErrors* sf_eta0_id = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_Tight_pt_abseta<0.9");
//   TGraphAsymmErrors* sf_eta1_id = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_Tight_pt_abseta0.9-1.2");
//   TGraphAsymmErrors* sf_eta2_id = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_Tight_pt_abseta1.2-2.1");

//   sf_eta0_id->SetName("SF_eta0_ID_tight");
//   sf_eta1_id->SetName("SF_eta1_ID_tight");
//   sf_eta2_id->SetName("SF_eta2_ID_tight");

//   // isolation
//   TGraphAsymmErrors* sf_eta0_iso = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta<0.9");
//   TGraphAsymmErrors* sf_eta1_iso = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta0.9-1.2");
//   TGraphAsymmErrors* sf_eta2_iso = (TGraphAsymmErrors*) file_c->Get("DATA_over_MC_combRelIsoPF04dBeta<02_Tight_pt_abseta1.2-2.1");

//   sf_eta0_iso->SetName("SF_eta0_ISO_tight");
//   sf_eta1_iso->SetName("SF_eta1_ISO_tight");
//   sf_eta2_iso->SetName("SF_eta2_ISO_tight");

//   // trigger: IsoMu24
//   TGraphAsymmErrors* sf_eta0_trig_isomu24 = (TGraphAsymmErrors*) file_d->Get("DATA_over_MC_IsoMu24_eta2p1_TightIso_pt_abseta<0.9_2012D");
//   TGraphAsymmErrors* sf_eta1_trig_isomu24 = (TGraphAsymmErrors*) file_d->Get("DATA_over_MC_IsoMu24_eta2p1_TightIso_pt_abseta0.9-1.2_2012D");
//   TGraphAsymmErrors* sf_eta2_trig_isomu24 = (TGraphAsymmErrors*) file_d->Get("DATA_over_MC_IsoMu24_eta2p1_TightIso_pt_abseta1.2-2.1_2012D");

//   // attention: correct unphysical error in bin 4:
//   sf_eta0_trig_isomu24->SetPointEYhigh(3, sf_eta0_trig_isomu24->GetErrorYlow(3));

//   sf_eta0_trig_isomu24->SetName("SF_eta0_TRIG_IsoMu24");
//   sf_eta1_trig_isomu24->SetName("SF_eta1_TRIG_IsoMu24");
//   sf_eta2_trig_isomu24->SetName("SF_eta2_TRIG_IsoMu24");

//   // trigger: Mu40
//   TGraphAsymmErrors* sf_eta0_trig_mu40 = (TGraphAsymmErrors*) file_d->Get("DATA_over_MC_Mu40_eta2p1_Tight_pt_40_abseta<0.9_2012D");
//   TGraphAsymmErrors* sf_eta1_trig_mu40 = (TGraphAsymmErrors*) file_d->Get("DATA_over_MC_Mu40_eta2p1_Tight_pt_40_abseta0.9-1.2_2012D");
//   TGraphAsymmErrors* sf_eta2_trig_mu40 = (TGraphAsymmErrors*) file_d->Get("DATA_over_MC_Mu40_eta2p1_Tight_pt_40_abseta1.2-2.1_2012D");

//   sf_eta0_trig_mu40->SetName("SF_eta0_TRIG_Mu40");
//   sf_eta1_trig_mu40->SetName("SF_eta1_TRIG_Mu40");
//   sf_eta2_trig_mu40->SetName("SF_eta2_TRIG_Mu40");

//   d_d->cd();
//   sf_eta0_id->Write();
//   sf_eta1_id->Write();
//   sf_eta2_id->Write();

//   sf_eta0_iso->Write();
//   sf_eta1_iso->Write();
//   sf_eta2_iso->Write();

//   sf_eta0_trig_isomu24->Write();
//   sf_eta1_trig_isomu24->Write();
//   sf_eta2_trig_isomu24->Write();

//   sf_eta0_trig_mu40->Write();
//   sf_eta1_trig_mu40->Write();
//   sf_eta2_trig_mu40->Write();


  tf->Write();
  tf->Close();
  

}
