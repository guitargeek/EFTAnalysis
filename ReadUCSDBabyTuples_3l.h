#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <TGraphAsymmErrors.h>
#include <TVector3.h>
#include <TRandom.h>

using std::string;

double MUON_MASS = 105.6583745*10e-03;
double ELECTRON_MASS = 511*10e-06;


bool sameValHighPrecision(double a, double b)
{
   return fabs(a - b) < 1.000e-08;
}

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double mass)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiM(pT, eta, phi, mass);
  return object_p4;
}

typedef struct
{
  float pt; 
  float phi;
  float eta;
  float pfiso;
  float trkiso;
  int charge;
  float dz;
  float d0;
  float sip3d;
  bool id;
  float recoEW;
  float triggerEW;
  float trigger;
  bool isTightUCSD;
  bool isTightMIT;
  bool isveryLooseUCSD;
  bool isHLTsafeMIT;
} leptonInfo;

typedef struct
{
  float jetPt;
  float jetEta;
  float jetPhi;
  float jetMass;
  float jetCSV;
} jetInfo;

typedef struct
{
  float ak8JetPrunmass;
  float ak8JetTrimmass;
  float ak8Jetsd0;
  float ak8JetPt;
  float ak8JetEta;
  float ak8JetPhi;
  float ak8JetTau1;
  float ak8JetTau2;
} fatJetInfo;

typedef struct
{
    TLorentzVector JetLV;
    float BTag_CSV;
} AnalysisJetInfo;

bool sortLeptonsInDescendingpT(leptonInfo lep1, leptonInfo lep2)
{
  return (lep1.pt > lep2.pt);
}

bool sortJetVectorsInDescendingpT(AnalysisJetInfo jet1, AnalysisJetInfo jet2)
{
  return (jet1.JetLV.Pt() > jet2.JetLV.Pt());
}

bool sortJetLorentzVectorsInDescendingpT(TLorentzVector jet1, TLorentzVector jet2)
{
  return (jet1.Pt() > jet2.Pt());
}


bool sortFatJetVectorsInDescendingpT(fatJetInfo fatjet1, fatJetInfo fatjet2)
{
  return (fatjet1.ak8JetPt > fatjet2.ak8JetPt);
}

typedef struct
{
  TH1F *h_mu_pt_leading;
  TH1F *h_mu_pt_trailing;
  TH1F *h_mu_eta_leading;
  TH1F *h_mu_eta_trailing;
  TH1F *h_mu_phi_leading;
  TH1F *h_mu_phi_trailing;
  TH1F *h_mu_pTRatio_leading;
  TH1F *h_mu_pTRatio_trailing;
  TH1F *h_mu_ip3d_leading;
  TH1F *h_mu_ip3d_trailing;
  TH1F *h_mu_d0_leading;
  TH1F *h_mu_d0_trailing;
  TH1F *h_mu_dz_leading;
  TH1F *h_mu_dz_trailing;
  TH1F *h_mu_pTrel_leading;
  TH1F *h_mu_pTrel_trailing;
  TH1F *h_mu_relIso_leading;
  TH1F *h_mu_relIso_trailing;
  TH1F *h_mu_deltaR_leading;
  TH1F *h_mu_deltaR_trailing;
  TH1F *h_InvariantMass;
  TH1F *h_MET;
  TH1F *h_nJets;
  TH1F *h_nbJets;
  TH1F *h_InvariantMassJJ;
  TH1F *h_deltaEtaJJ;
  TH1F *h_jet_pt_leading;
  TH1F *h_jet_pt_trailing;
  TH1F *h_jet_eta_leading;
  TH1F *h_jet_eta_trailing;
  TH1F *h_jet_phi_leading;
  TH1F *h_jet_phi_trailing;
  TH1F *h_jet_csv_leading;
  TH1F *h_jet_csv_trailing;
  TH1F *h_bjet_pt_leading;
  TH1F *h_bjet_eta_leading;
  TH1F *h_bjet_phi_leading;
  TH1F *h_bjet_pt_trailing;
  TH1F *h_bjet_eta_trailing;
  TH1F *h_bjet_phi_trailing;
  TH1F *h_el_pt_leading;
  TH1F *h_el_pt_trailing;
  TH1F *h_el_eta_leading;
  TH1F *h_el_eta_trailing;
  TH1F *h_el_phi_leading;
  TH1F *h_el_phi_trailing;
  TH1F *h_el_pTRatio_leading;
  TH1F *h_el_pTRatio_trailing;
  TH1F *h_el_ip3d_leading;
  TH1F *h_el_ip3d_trailing;
  TH1F *h_el_d0_leading;
  TH1F *h_el_d0_trailing;
  TH1F *h_el_dz_leading;
  TH1F *h_el_dz_trailing;
  TH1F *h_el_pTrel_leading;
  TH1F *h_el_pTrel_trailing;
  TH1F *h_el_relIso_leading;
  TH1F *h_el_relIso_trailing;
  TH1F *h_el_deltaR_leading;
  TH1F *h_el_deltaR_trailing;
  TH1F *h_nfatJets;
  TH1F *h_fatjet_pt_leading;
  TH1F *h_fatjet_eta_leading;
  TH1F *h_fatjet_phi_leading;
  TH1F *h_fatjet_prunMass_leading;
  TH1F *h_fatjet_trimMass_leading;
  TH1F *h_fatjet_sd0Mass_leading;
  TH1F *h_fatjet_tau21_leading;
  TH1F *h_ST;
} HistCollection;

void initializeHistCollection(HistCollection &histCol, std::string suffix)
{
   histCol.h_mu_pt_leading = new TH1F(("h_mu_pt_leading_"+suffix).c_str(), "Leading muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_pt_leading->Sumw2();
   histCol.h_mu_pt_trailing = new TH1F(("h_mu_pt_trailing_"+suffix).c_str(), "Trailing muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_pt_trailing->Sumw2();
   histCol.h_mu_eta_leading = new TH1F(("h_mu_eta_leading_"+suffix).c_str(), "Leading muon #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_mu_eta_leading->Sumw2();
   histCol.h_mu_eta_trailing = new TH1F(("h_mu_eta_trailing_"+suffix).c_str(), "Trailing muon #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_mu_eta_trailing->Sumw2();
   histCol.h_mu_phi_leading = new TH1F(("h_mu_phi_leading_"+suffix).c_str(), "Leading muon #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_mu_phi_leading->Sumw2();
   histCol.h_mu_phi_trailing = new TH1F(("h_mu_phi_trailing_"+suffix).c_str(), "Trailing muon #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_mu_phi_trailing->Sumw2();
   histCol.h_mu_pTRatio_leading = new TH1F(("h_mu_pTRatio_leading_"+suffix).c_str(), "Leading muon pT ratio ; pT ratio; Events", 150, 0.0, 1.5); histCol.h_mu_pTRatio_leading->Sumw2();
   histCol.h_mu_pTRatio_trailing = new TH1F(("h_mu_pTRatio_trailing_"+suffix).c_str(), "Trailing muon pT ratio ; pT ratio; Events", 150, 0.0, 1.5); histCol.h_mu_pTRatio_trailing->Sumw2();
   histCol.h_mu_ip3d_leading = new TH1F(("h_mu_ip3d_leading_"+suffix).c_str(), "Leading muon 3D impact parameter; Leading muon 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_ip3d_leading->Sumw2();
   histCol.h_mu_ip3d_trailing = new TH1F(("h_mu_ip3d_trailing_"+suffix).c_str(), "Trailing muon 3D impact parameter; Trailing muon 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_ip3d_trailing->Sumw2();
   histCol.h_mu_d0_leading = new TH1F(("h_mu_d0_leading_"+suffix).c_str(), "Leading muon d0 impact parameter; Leading muon d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_d0_leading->Sumw2();
   histCol.h_mu_d0_trailing = new TH1F(("h_mu_d0_trailing_"+suffix).c_str(), "Trailing muon d0 impact parameter; Trailing muon d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_d0_trailing->Sumw2();
   histCol.h_mu_dz_leading = new TH1F(("h_mu_dz_leading_"+suffix).c_str(), "Leading muon dz impact parameter; Leading muon dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_dz_leading->Sumw2();
   histCol.h_mu_dz_trailing = new TH1F(("h_mu_dz_trailing_"+suffix).c_str(), "Trailing muon dz impact parameter; Trailing muon dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_dz_trailing->Sumw2();
   histCol.h_mu_pTrel_leading = new TH1F(("h_mu_pTrel_leading_"+suffix).c_str(), "Leading muon pT_{rel}; Leading muon pT_{rel} [GeV]; Events", 100000, 0.0, 1000.0);histCol.h_mu_pTrel_leading->Sumw2();
   histCol.h_mu_pTrel_trailing = new TH1F(("h_mu_pTrel_trailing_"+suffix).c_str(), "Trailing muon pT_{rel}; Trailing muon pT_{rel} [GeV]; Events", 100000, 0.0, 1000.0);histCol.h_mu_pTrel_trailing->Sumw2();
   histCol.h_mu_relIso_leading = new TH1F(("h_mu_relIso_leading_"+suffix).c_str(), "Leading relative isolation; Leading relative isolation; Events", 10000, 0.0, 1.0);histCol.h_mu_relIso_leading->Sumw2();
   histCol.h_mu_relIso_trailing = new TH1F(("h_mu_relIso_trailing_"+suffix).c_str(), "Trailing relative isolation; Trailing relative isolation; Events", 10000, 0.0, 1.0);histCol.h_mu_relIso_trailing->Sumw2();
   histCol.h_mu_deltaR_leading = new TH1F(("h_mu_deltaR_leading_"+suffix).c_str(), "Leading muon-jet #DeltaR; Leading muon-jet #DeltaR; Events", 10000, 0.0, 10.0);histCol.h_mu_deltaR_leading->Sumw2();
   histCol.h_mu_deltaR_trailing = new TH1F(("h_mu_deltaR_trailing_"+suffix).c_str(), "Trailing muon-jet #DeltaR; Trailing muon-jet #DeltaR; Events", 10000, 0.0, 10.0);histCol.h_mu_deltaR_trailing->Sumw2();
   histCol.h_InvariantMass=new TH1F(("h_InvariantMass_"+suffix).c_str(), "Di-lepton invariant mass; M_{ll} [GeV]; Events/GeV", 10000, 0.0, 1000.0); histCol.h_InvariantMass->Sumw2();
   histCol.h_MET=new TH1F(("h_MET_"+suffix).c_str(), "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); histCol.h_MET->Sumw2();
   histCol.h_nJets = new TH1F(("h_nJets_"+suffix).c_str(), "Number of Jets; Number of Jets; Events", 20, -0.5, 19.5);histCol.h_nJets->Sumw2();
   histCol.h_nbJets = new TH1F(("h_nbJets_"+suffix).c_str(), "Number of b-Jets; Number of b-Jets; Events", 20, -0.5, 19.5);histCol.h_nbJets->Sumw2();
   histCol.h_InvariantMassJJ = new TH1F(("h_InvariantMassJJ_"+suffix).c_str(), "Di-jet invariant mass M_{jj} [GeV]; Events/GeV", 50000.0, 0.0, 5000.0);histCol.h_InvariantMassJJ->Sumw2();
   histCol.h_deltaEtaJJ = new TH1F(("h_deltaEtaJJ_"+suffix).c_str(), "#Delta #eta_{jj}; Events", 1200.0, -6.0, 6.0);histCol.h_deltaEtaJJ->Sumw2();
   histCol.h_jet_pt_leading=new TH1F(("h_jet_pt_leading_"+suffix).c_str(), "Leading jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_pt_leading->Sumw2();
   histCol.h_jet_pt_trailing=new TH1F(("h_jet_pt_trailing_"+suffix).c_str(), "Trailing jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_pt_trailing->Sumw2();
   histCol.h_jet_eta_leading=new TH1F(("h_jet_eta_leading_"+suffix).c_str(), "Leading jet #eta; #eta; Events", 1200, -6.0, 6.0); histCol.h_jet_eta_leading->Sumw2();
   histCol.h_jet_eta_trailing=new TH1F(("h_jet_eta_trailing_"+suffix).c_str(), "Trailing jet #eta; #eta; Events", 1200, -6.0, 6.0); histCol.h_jet_eta_trailing->Sumw2();
   histCol.h_jet_phi_leading=new TH1F(("h_jet_phi_leading_"+suffix).c_str(), "Leading jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_leading->Sumw2();
   histCol.h_jet_phi_trailing=new TH1F(("h_jet_phi_trailing_"+suffix).c_str(), "Trailing jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_trailing->Sumw2();
   histCol.h_jet_csv_leading=new TH1F(("h_jet_csv_leading_"+suffix).c_str(), "Leading jet CSV discriminator; CSV discriminator; Events", 100, 0, 1.0);histCol.h_jet_csv_leading->Sumw2();
   histCol.h_jet_csv_trailing=new TH1F(("h_jet_csv_trailing_"+suffix).c_str(), "Trailing jet CSV discriminator; CSV discriminator; Events", 100, 0, 1.0);histCol.h_jet_csv_trailing->Sumw2();
   histCol.h_bjet_pt_leading=new TH1F(("h_bjet_pt_leading_"+suffix).c_str(), "Leading b-jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000);histCol.h_bjet_pt_leading->Sumw2();
   histCol.h_bjet_eta_leading=new TH1F(("h_bjet_eta_leading_"+suffix).c_str(), "Leading b-jet |#eta|; |#eta|; Events", 1200, -6.0, 6.0);histCol.h_bjet_eta_leading->Sumw2();
   histCol.h_bjet_phi_leading=new TH1F(("h_bjet_phi_leading_"+suffix).c_str(), "Leading b-jet #phi;  #phi; Events", 800, -4.0, 4.0);histCol.h_bjet_phi_leading->Sumw2();
   histCol.h_bjet_pt_trailing=new TH1F(("h_bjet_pt_trailing_"+suffix).c_str(), "Trailing b-jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000);histCol.h_bjet_pt_trailing->Sumw2();
   histCol.h_bjet_eta_trailing=new TH1F(("h_bjet_eta_trailing_"+suffix).c_str(), "Trailing b-jet |#eta|; |#eta|; Events", 1200, -6.0, 6.0);histCol.h_bjet_eta_trailing->Sumw2();
   histCol.h_bjet_phi_trailing=new TH1F(("h_bjet_phi_trailing_"+suffix).c_str(), "Trailing b-jet #phi;  #phi; Events", 800, -4.0, 4.0);histCol.h_bjet_phi_trailing->Sumw2();
   histCol.h_el_pt_leading = new TH1F(("h_el_pt_leading_"+suffix).c_str(), "Leading electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_pt_leading->Sumw2();     histCol.h_el_pt_trailing = new TH1F(("h_el_pt_trailing_"+suffix).c_str(), "Trailing electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_pt_trailing->Sumw2();
   histCol.h_el_eta_leading = new TH1F(("h_el_eta_leading_"+suffix).c_str(), "Leading electron #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_el_eta_leading->Sumw2();
   histCol.h_el_eta_trailing = new TH1F(("h_el_eta_trailing_"+suffix).c_str(), "Trailing electron #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_el_eta_trailing->Sumw2();
   histCol.h_el_phi_leading = new TH1F(("h_el_phi_leading_"+suffix).c_str(), "Leading electron #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_el_phi_leading->Sumw2();
   histCol.h_el_phi_trailing = new TH1F(("h_el_phi_trailing_"+suffix).c_str(), "Trailing electron #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_el_phi_trailing->Sumw2();
   histCol.h_el_pTRatio_leading = new TH1F(("h_el_pTRatio_leading_"+suffix).c_str(), "Leading electron pT ratio ; pT ratio; Events", 150, 0.0, 1.5); histCol.h_el_pTRatio_leading->Sumw2();
   histCol.h_el_pTRatio_trailing = new TH1F(("h_el_pTRatio_trailing_"+suffix).c_str(), "Trailing electron pT ratio ; pT ratio; Events", 150, 0.0, 1.5); histCol.h_el_pTRatio_trailing->Sumw2();
   histCol.h_el_ip3d_leading = new TH1F(("h_el_ip3d_leading_"+suffix).c_str(), "Leading electron 3D impact parameter; Leading electron 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_ip3d_leading->Sumw2();
   histCol.h_el_ip3d_trailing = new TH1F(("h_el_ip3d_trailing_"+suffix).c_str(), "Trailing electron 3D impact parameter; Trailing electron 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_ip3d_trailing->Sumw2();
   histCol.h_el_d0_leading = new TH1F(("h_el_d0_leading_"+suffix).c_str(), "Leading electron d0 impact parameter; Leading electron d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_d0_leading->Sumw2();
   histCol.h_el_d0_trailing = new TH1F(("h_el_d0_trailing_"+suffix).c_str(), "Trailing electron d0 impact parameter; Trailing electron d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_d0_trailing->Sumw2();
   histCol.h_el_dz_leading = new TH1F(("h_el_dz_leading_"+suffix).c_str(), "Leading electron dz impact parameter; Leading electron dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_dz_leading->Sumw2();
   histCol.h_el_dz_trailing = new TH1F(("h_el_dz_trailing_"+suffix).c_str(), "Trailing electron dz impact parameter; Trailing electron dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_dz_trailing->Sumw2();
   histCol.h_el_pTrel_leading = new TH1F(("h_el_pTrel_leading_"+suffix).c_str(), "Leading electron pT_{rel}; Leading electron pT_{rel} [GeV]; Events", 100000, 0.0, 1000.0);histCol.h_el_pTrel_leading->Sumw2();
   histCol.h_el_pTrel_trailing = new TH1F(("h_el_pTrel_trailing_"+suffix).c_str(), "Trailing electron pT_{rel}; Trailing electron pT_{rel} [GeV]; Events", 100000, 0.0, 1000.0);histCol.h_el_pTrel_trailing->Sumw2();
   histCol.h_el_relIso_leading = new TH1F(("h_el_relIso_leading_"+suffix).c_str(), "Leading relative isolation; Leading relative isolation; Events", 10000, 0.0, 1.0);histCol.h_el_relIso_leading->Sumw2();
   histCol.h_el_relIso_trailing = new TH1F(("h_el_relIso_trailing_"+suffix).c_str(), "Trailing relative isolation; Trailing relative isolation; Events", 10000, 0.0, 1.0);histCol.h_el_relIso_trailing->Sumw2();
   histCol.h_el_deltaR_leading = new TH1F(("h_el_deltaR_leading_"+suffix).c_str(), "Leading electron-jet #DeltaR; Leading electron-jet #DeltaR; Events", 10000, 0.0, 10.0);histCol.h_el_deltaR_leading->Sumw2();
   histCol.h_el_deltaR_trailing = new TH1F(("h_el_deltaR_trailing_"+suffix).c_str(), "Trailing electron-jet #DeltaR; Trailing electron-jet #DeltaR; Events", 10000, 0.0, 10.0);histCol.h_el_deltaR_trailing->Sumw2();
   histCol.h_fatjet_pt_leading=new TH1F(("h_fatjet_pt_leading_"+suffix).c_str(), "Leading fatjet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_fatjet_pt_leading->Sumw2();
   histCol.h_fatjet_phi_leading=new TH1F(("h_fatjet_phi_leading_"+suffix).c_str(), "Leading fatjet #phi; #phi; Events", 800, -4.0, 4.0); histCol.h_fatjet_phi_leading->Sumw2();
   histCol.h_fatjet_eta_leading=new TH1F(("h_fatjet_eta_leading_"+suffix).c_str(), "Leading fatjet #eta; #eta; Events", 600, -3.0, 3.0); histCol.h_fatjet_eta_leading->Sumw2();
   histCol.h_nfatJets = new TH1F(("h_nfatJets_"+suffix).c_str(), "Number of fatJets; Number of fatJets; Events", 20, -0.5, 19.5);histCol.h_nfatJets->Sumw2();
   histCol.h_fatjet_prunMass_leading = new TH1F(("h_fatjet_prunMass_leading_"+suffix).c_str(), "Pruned Mass [GeV]; Pruned Mass [GeV]; Events", 10000, 0, 1000);histCol.h_fatjet_prunMass_leading->Sumw2();
   histCol.h_fatjet_trimMass_leading = new TH1F(("h_fatjet_trimMass_leading_"+suffix).c_str(), "Trimmed Mass [GeV]; Trimmed Mass [GeV]; Events", 10000, 0, 1000);histCol.h_fatjet_trimMass_leading->Sumw2();
   histCol.h_fatjet_sd0Mass_leading = new TH1F(("h_fatjet_sd0Mass_leading_"+suffix).c_str(), "Soft drop Mass [GeV]; soft drop Mass [GeV]; Events", 10000, 0, 1000);histCol.h_fatjet_sd0Mass_leading->Sumw2();
   histCol.h_fatjet_tau21_leading = new TH1F(("h_fatjet_tau21_leading_"+suffix).c_str(), "N-subjettiness; #tau2/#tau1; Events", 1000, 0, 1.0);histCol.h_fatjet_tau21_leading->Sumw2();
   histCol.h_ST = new TH1F(("h_ST_"+suffix).c_str(), "ST; ST (sum of lepton pT + jet pT + MET); Events/GeV", 5000, 0, 5000.0);histCol.h_ST->Sumw2();
}

void writeHistCollection(HistCollection &histCol)
{
   histCol.h_mu_pt_leading->Write();
   histCol.h_mu_pt_trailing->Write();
   histCol.h_mu_eta_leading->Write();
   histCol.h_mu_eta_trailing->Write();
   histCol.h_mu_phi_leading->Write();
   histCol.h_mu_phi_trailing->Write();
   histCol.h_mu_pTRatio_leading->Write();
   histCol.h_mu_pTRatio_trailing->Write();
   histCol.h_mu_ip3d_leading->Write();
   histCol.h_mu_ip3d_trailing->Write();
   histCol.h_mu_d0_leading->Write();
   histCol.h_mu_d0_trailing->Write();
   histCol.h_mu_dz_leading->Write();
   histCol.h_mu_dz_trailing->Write();
   histCol.h_mu_pTrel_leading->Write();
   histCol.h_mu_pTrel_trailing->Write();
   histCol.h_mu_relIso_leading->Write();
   histCol.h_mu_relIso_trailing->Write();
   histCol.h_mu_deltaR_leading->Write();
   histCol.h_mu_deltaR_trailing->Write();
   histCol.h_InvariantMass->Write();
   histCol.h_MET->Write();
   histCol.h_nJets->Write();
   histCol.h_nbJets->Write();
   histCol.h_InvariantMassJJ->Write();
   histCol.h_deltaEtaJJ->Write();
   histCol.h_jet_pt_leading->Write();
   histCol.h_jet_pt_trailing->Write();
   histCol.h_jet_eta_leading->Write();
   histCol.h_jet_eta_trailing->Write();
   histCol.h_jet_phi_leading->Write();
   histCol.h_jet_csv_leading->Write();
   histCol.h_jet_csv_trailing->Write();
   histCol.h_bjet_pt_leading->Write();
   histCol.h_bjet_eta_leading->Write();
   histCol.h_bjet_phi_leading->Write();
   histCol.h_bjet_pt_trailing->Write();
   histCol.h_bjet_eta_trailing->Write();
   histCol.h_bjet_phi_trailing->Write();
   histCol.h_el_pt_leading->Write();
   histCol.h_el_pt_trailing->Write();
   histCol.h_el_eta_leading->Write();
   histCol.h_el_eta_trailing->Write();
   histCol.h_el_phi_leading->Write();
   histCol.h_el_phi_trailing->Write();
   histCol.h_el_pTRatio_leading->Write();
   histCol.h_el_pTRatio_trailing->Write();
   histCol.h_el_ip3d_leading->Write();
   histCol.h_el_ip3d_trailing->Write();
   histCol.h_el_d0_leading->Write();
   histCol.h_el_d0_trailing->Write();
   histCol.h_el_dz_leading->Write();
   histCol.h_el_dz_trailing->Write();
   histCol.h_el_pTrel_leading->Write();
   histCol.h_el_pTrel_trailing->Write();
   histCol.h_el_relIso_leading->Write();
   histCol.h_el_relIso_trailing->Write();
   histCol.h_el_deltaR_leading->Write();
   histCol.h_el_deltaR_trailing->Write();
   histCol.h_nfatJets->Write();
   histCol.h_fatjet_pt_leading->Write();
   histCol.h_fatjet_phi_leading->Write();
   histCol.h_fatjet_eta_leading->Write();
   histCol.h_fatjet_prunMass_leading->Write();
   histCol.h_fatjet_trimMass_leading->Write();
   histCol.h_fatjet_sd0Mass_leading->Write();
   histCol.h_fatjet_tau21_leading->Write();
   histCol.h_ST->Write();
}

void fillMuHistCollection(HistCollection &histCol, double mu1pt, double mu2pt, double mu1eta,  double mu2eta, double mu1phi, double mu2phi, double mu1pTRatio, double mu2pTRatio, double mu1ip3D, double mu2ip3D, double mu1d0, double mu2d0, double mu1dz, double mu2dz, double mu1pTRel, double mu2pTRel, double mu1relIso, double mu2relIso, double mu1DeltaR, double mu2DeltaR, double invmass, double MET, int njets, int nbjets, double mjjL, double etajj, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double ST, double weight)
{
  histCol.h_mu_pt_leading->Fill(mu1pt, weight);
  histCol.h_mu_pt_trailing->Fill(mu2pt, weight);
  histCol.h_mu_eta_leading->Fill(mu1eta, weight);
  histCol.h_mu_eta_trailing->Fill(mu2eta, weight);
  histCol.h_mu_phi_leading->Fill(mu1phi, weight);
  histCol.h_mu_phi_trailing->Fill(mu2phi, weight);
  histCol.h_mu_pTRatio_leading->Fill(mu1pTRatio, weight);
  histCol.h_mu_pTRatio_trailing->Fill(mu2pTRatio, weight);
  histCol.h_mu_ip3d_leading->Fill(mu1ip3D, weight);
  histCol.h_mu_ip3d_trailing->Fill(mu2ip3D, weight);
  histCol.h_mu_d0_leading->Fill(mu1d0, weight);
  histCol.h_mu_d0_trailing->Fill(mu2d0, weight);
  histCol.h_mu_dz_leading->Fill(mu1dz, weight);
  histCol.h_mu_dz_trailing->Fill(mu2dz, weight);
  histCol.h_mu_pTrel_leading->Fill(mu1pTRel, weight);
  histCol.h_mu_pTrel_trailing->Fill(mu2pTRel, weight);
  histCol.h_mu_relIso_leading->Fill(mu1relIso, weight);
  histCol.h_mu_relIso_trailing->Fill(mu2relIso, weight);
  histCol.h_mu_deltaR_leading->Fill(mu1DeltaR, weight);
  histCol.h_mu_deltaR_trailing->Fill(mu2DeltaR, weight);
  histCol.h_InvariantMass->Fill(invmass, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  histCol.h_InvariantMassJJ->Fill(mjjL, weight);
  histCol.h_deltaEtaJJ->Fill(etajj, weight);
  if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
  if(njets >= 2) histCol.h_jet_pt_trailing->Fill(jet2pt, weight);     
  if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
  if(njets >= 2) histCol.h_jet_eta_trailing->Fill(jet2eta, weight);
  if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
  if(njets >= 2) histCol.h_jet_phi_trailing->Fill(jet2phi, weight);
  if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
  if(njets >= 2) histCol.h_jet_csv_trailing->Fill(jet2csv, weight);
  if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
  if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
  if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight); 
  if(nbjets >= 2) histCol.h_bjet_pt_trailing->Fill(bjet2pt, weight);
  if(nbjets >= 2) histCol.h_bjet_eta_trailing->Fill(bjet2eta, weight);
  if(nbjets >= 2) histCol.h_bjet_phi_trailing->Fill(bjet2phi, weight);
  histCol.h_ST->Fill(ST, weight);
}

void fillElHistCollection(HistCollection &histCol, double el1pt, double el2pt, double el1eta,  double el2eta, double el1phi, double el2phi,  double el1pTRatio, double el2pTRatio, double el1ip3D, double el2ip3D, double el1d0, double el2d0, double el1dz, double el2dz, double invmass, double MET, int njets, int nbjets, double mjjL, double etajj, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double ST, double weight)
{
  histCol.h_el_pt_leading->Fill(el1pt, weight);
  histCol.h_el_pt_trailing->Fill(el2pt, weight);
  histCol.h_el_eta_leading->Fill(el1eta, weight);
  histCol.h_el_eta_trailing->Fill(el2eta, weight);
  histCol.h_el_phi_leading->Fill(el1phi, weight);
  histCol.h_el_phi_trailing->Fill(el2phi, weight);
  histCol.h_el_pTRatio_leading->Fill(el1pTRatio, weight);
  histCol.h_el_pTRatio_trailing->Fill(el2pTRatio, weight); 
  histCol.h_el_ip3d_leading->Fill(el1ip3D, weight);
  histCol.h_el_ip3d_trailing->Fill(el2ip3D, weight);
  histCol.h_el_d0_leading->Fill(el1d0, weight);
  histCol.h_el_d0_trailing->Fill(el2d0, weight);
  histCol.h_el_dz_leading->Fill(el1dz, weight);
  histCol.h_el_dz_trailing->Fill(el2dz, weight);
  histCol.h_InvariantMass->Fill(invmass, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  histCol.h_InvariantMassJJ->Fill(mjjL, weight);
  histCol.h_deltaEtaJJ->Fill(etajj, weight);
  if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
  if(njets >= 2) histCol.h_jet_pt_trailing->Fill(jet2pt, weight);
  if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
  if(njets >= 2) histCol.h_jet_eta_trailing->Fill(jet2eta, weight);
  if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
  if(njets >= 2) histCol.h_jet_phi_trailing->Fill(jet2phi, weight);
  if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
  if(njets >= 2) histCol.h_jet_csv_trailing->Fill(jet2csv, weight);
  if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
  if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
  if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight);
  if(nbjets >= 2) histCol.h_bjet_pt_trailing->Fill(bjet2pt, weight);
  if(nbjets >= 2) histCol.h_bjet_eta_trailing->Fill(bjet2eta, weight);
  if(nbjets >= 2) histCol.h_bjet_phi_trailing->Fill(bjet2phi, weight);
  histCol.h_ST->Fill(ST, weight);
}

void fillElMuHistCollection(HistCollection &histCol, double el1pt, double mu1pt, double el1eta,  double mu1eta, double el1phi, double mu1phi, double el1pTRatio, double mu1pTRatio, double el1ip3D, double mu1ip3D, double el1d0, double mu1d0, double el1dz, double mu1dz, double invmass, double MET, int njets, int nbjets, double mjjL, double etajj, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double ST, double weight)
{
  histCol.h_el_pt_leading->Fill(el1pt, weight);
  histCol.h_mu_pt_leading->Fill(mu1pt, weight);
  histCol.h_el_eta_leading->Fill(el1eta, weight);
  histCol.h_mu_eta_leading->Fill(mu1eta, weight);
  histCol.h_el_phi_leading->Fill(el1phi, weight);
  histCol.h_mu_phi_leading->Fill(mu1phi, weight);
  histCol.h_el_pTRatio_leading->Fill(el1pTRatio, weight);
  histCol.h_mu_pTRatio_leading->Fill(mu1pTRatio, weight);
  histCol.h_el_ip3d_leading->Fill(el1ip3D, weight);
  histCol.h_mu_ip3d_leading->Fill(mu1ip3D, weight);
  histCol.h_el_d0_leading->Fill(el1d0, weight);
  histCol.h_mu_d0_leading->Fill(mu1d0, weight);
  histCol.h_el_dz_leading->Fill(el1dz, weight);
  histCol.h_mu_dz_leading->Fill(mu1dz, weight);
  histCol.h_InvariantMass->Fill(invmass, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  histCol.h_InvariantMassJJ->Fill(mjjL, weight);
  histCol.h_deltaEtaJJ->Fill(etajj, weight);
  if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
  if(njets >= 2) histCol.h_jet_pt_trailing->Fill(jet2pt, weight);
  if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
  if(njets >= 2) histCol.h_jet_eta_trailing->Fill(jet2eta, weight);
  if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
  if(njets >= 2) histCol.h_jet_phi_trailing->Fill(jet2phi, weight);
  if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
  if(njets >= 2) histCol.h_jet_csv_trailing->Fill(jet2csv, weight);
  if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
  if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
  if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight);
  if(nbjets >= 2) histCol.h_bjet_pt_trailing->Fill(bjet2pt, weight);
  if(nbjets >= 2) histCol.h_bjet_eta_trailing->Fill(bjet2eta, weight);
  if(nbjets >= 2) histCol.h_bjet_phi_trailing->Fill(bjet2phi, weight);
  histCol.h_ST->Fill(ST, weight);
}

void fillMuHistCollectionWithFatJet(HistCollection &histCol, double mu1pt, double mu2pt, double mu1eta,  double mu2eta, double mu1phi, double mu2phi, double invmass, double MET, int njets, int nbjets, int nfatjets, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double fatJetPt, double fatJetEta, double fatJetPhi, double prunMass, double trimMass, double sd0, double tau21, double weight)
{
  histCol.h_mu_pt_leading->Fill(mu1pt, weight);
  histCol.h_mu_pt_trailing->Fill(mu2pt, weight);
  histCol.h_mu_eta_leading->Fill(mu1eta, weight);
  histCol.h_mu_eta_trailing->Fill(mu2eta, weight);
  histCol.h_mu_phi_leading->Fill(mu1phi, weight);
  histCol.h_mu_phi_trailing->Fill(mu2phi, weight);
  histCol.h_InvariantMass->Fill(invmass, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  histCol.h_nfatJets->Fill(nfatjets, weight);
  if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
  if(njets >= 2) histCol.h_jet_pt_trailing->Fill(jet2pt, weight); 
  if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
  if(njets >= 2) histCol.h_jet_eta_trailing->Fill(jet2eta, weight);
  if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
  if(njets >= 2) histCol.h_jet_phi_trailing->Fill(jet2phi, weight);
  if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
  if(njets >= 2) histCol.h_jet_csv_trailing->Fill(jet2csv, weight);
  if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
  if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
  if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight); 
  if(nbjets >= 2) histCol.h_bjet_pt_trailing->Fill(bjet2pt, weight);
  if(nbjets >= 2) histCol.h_bjet_eta_trailing->Fill(bjet2eta, weight);
  if(nbjets >= 2) histCol.h_bjet_phi_trailing->Fill(bjet2phi, weight); 
  if(nfatjets > 0) histCol.h_fatjet_pt_leading->Fill(fatJetPt, weight);   
  if(nfatjets > 0) histCol.h_fatjet_eta_leading->Fill(fatJetEta, weight);
  if(nfatjets > 0) histCol.h_fatjet_phi_leading->Fill(fatJetPhi, weight);
  if(nfatjets > 0) histCol.h_fatjet_prunMass_leading->Fill(prunMass, weight);
  if(nfatjets > 0) histCol.h_fatjet_trimMass_leading->Fill(trimMass, weight);
  if(nfatjets > 0) histCol.h_fatjet_sd0Mass_leading->Fill(sd0, weight);
  if(nfatjets > 0) histCol.h_fatjet_tau21_leading->Fill(tau21, weight); 
}

void fillElHistCollectionWithFatJet(HistCollection &histCol, double el1pt, double el2pt, double el1eta,  double el2eta, double el1phi, double el2phi, double invmass, double MET, int njets, int nbjets, int nfatjets, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double fatJetPt, double fatJetEta, double fatJetPhi, double prunMass, double trimMass, double sd0, double tau21, double weight)
{
  histCol.h_el_pt_leading->Fill(el1pt, weight);
  histCol.h_el_pt_trailing->Fill(el2pt, weight);
  histCol.h_el_eta_leading->Fill(el1eta, weight);
  histCol.h_el_eta_trailing->Fill(el2eta, weight);
  histCol.h_el_phi_leading->Fill(el1phi, weight);
  histCol.h_el_phi_trailing->Fill(el2phi, weight);
  histCol.h_InvariantMass->Fill(invmass, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  histCol.h_nfatJets->Fill(nfatjets, weight);
  if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
  if(njets >= 2) histCol.h_jet_pt_trailing->Fill(jet2pt, weight);
  if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
  if(njets >= 2) histCol.h_jet_eta_trailing->Fill(jet2eta, weight);
  if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
  if(njets >= 2) histCol.h_jet_phi_trailing->Fill(jet2phi, weight);
  if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
  if(njets >= 2) histCol.h_jet_csv_trailing->Fill(jet2csv, weight);
  if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
  if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
  if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight);
  if(nbjets >= 2) histCol.h_bjet_pt_trailing->Fill(bjet2pt, weight);
  if(nbjets >= 2) histCol.h_bjet_eta_trailing->Fill(bjet2eta, weight);
  if(nbjets >= 2) histCol.h_bjet_phi_trailing->Fill(bjet2phi, weight);
  if(nfatjets >= 1) histCol.h_fatjet_pt_leading->Fill(fatJetPt, weight);
  if(nfatjets >= 1) histCol.h_fatjet_eta_leading->Fill(fatJetEta, weight);
  if(nfatjets >= 1) histCol.h_fatjet_phi_leading->Fill(fatJetPhi, weight);
  if(nfatjets > 0) histCol.h_fatjet_prunMass_leading->Fill(prunMass, weight);
  if(nfatjets > 0) histCol.h_fatjet_trimMass_leading->Fill(trimMass, weight);
  if(nfatjets > 0) histCol.h_fatjet_sd0Mass_leading->Fill(sd0, weight);
  if(nfatjets > 0) histCol.h_fatjet_tau21_leading->Fill(tau21, weight);
}

void fillElMuHistCollectionWithFatJet(HistCollection &histCol, double el1pt, double mu1pt, double el1eta,  double mu1eta, double el1phi, double mu1phi, double invmass, double MET, int njets, int nbjets, int nfatjets, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double fatJetPt, double fatJetEta, double fatJetPhi,  double prunMass, double trimMass, double sd0, double tau21, double weight)
{
  histCol.h_el_pt_leading->Fill(el1pt, weight);
  histCol.h_mu_pt_leading->Fill(mu1pt, weight);
  histCol.h_el_eta_leading->Fill(el1eta, weight);
  histCol.h_mu_eta_leading->Fill(mu1eta, weight);
  histCol.h_el_phi_leading->Fill(el1phi, weight);
  histCol.h_mu_phi_leading->Fill(mu1phi, weight);
  histCol.h_InvariantMass->Fill(invmass, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  histCol.h_nfatJets->Fill(nfatjets, weight);
  if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
  if(njets >= 2) histCol.h_jet_pt_trailing->Fill(jet2pt, weight);
  if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
  if(njets >= 2) histCol.h_jet_eta_trailing->Fill(jet2eta, weight);
  if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
  if(njets >= 2) histCol.h_jet_phi_trailing->Fill(jet2phi, weight);
  if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
  if(njets >= 2) histCol.h_jet_csv_trailing->Fill(jet2csv, weight);
  if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
  if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
  if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight);
  if(nbjets >= 2) histCol.h_bjet_pt_trailing->Fill(bjet2pt, weight);
  if(nbjets >= 2) histCol.h_bjet_eta_trailing->Fill(bjet2eta, weight);
  if(nbjets >= 2) histCol.h_bjet_phi_trailing->Fill(bjet2phi, weight);
  if(nfatjets >= 1) histCol.h_fatjet_pt_leading->Fill(fatJetPt, weight);
  if(nfatjets >= 1) histCol.h_fatjet_eta_leading->Fill(fatJetEta, weight);
  if(nfatjets >= 1) histCol.h_fatjet_phi_leading->Fill(fatJetPhi, weight);
  if(nfatjets > 0) histCol.h_fatjet_prunMass_leading->Fill(prunMass, weight);
  if(nfatjets > 0) histCol.h_fatjet_trimMass_leading->Fill(trimMass, weight);
  if(nfatjets > 0) histCol.h_fatjet_sd0Mass_leading->Fill(sd0, weight);
  if(nfatjets > 0) histCol.h_fatjet_tau21_leading->Fill(tau21, weight);
}

double weightMu(std::vector<leptonInfo> v_mu)
{
  return v_mu.at(0).recoEW*v_mu.at(1).recoEW*(1 - (1 - v_mu.at(0).triggerEW)*(1 - v_mu.at(1).triggerEW));
}

double weightEl(std::vector<leptonInfo> v_el)
{
  return v_el.at(0).recoEW*v_el.at(1).recoEW*(1 - (1 - v_el.at(0).triggerEW)*(1 - v_el.at(1).triggerEW));
}

double weightElMu(std::vector<leptonInfo> v_el, std::vector<leptonInfo> v_mu)
{
  return v_el.at(0).recoEW*v_mu.at(0).recoEW*(1 - (1 - v_el.at(0).triggerEW)*(1 - v_mu.at(0).triggerEW));
}

double trigSFMuLead(TFile* fileData, double pt, double eta)
{
   TH2F *h_dimu_trig_SF = (TH2F*) fileData->Get("mu_lead_leg_eta_v_pt_sf");
   TAxis *xaxis = h_dimu_trig_SF->GetXaxis();
   TAxis *yaxis = h_dimu_trig_SF->GetYaxis();
   Int_t binx = xaxis->FindBin(eta);
   Int_t biny = yaxis->FindBin(pt);
   return h_dimu_trig_SF->GetBinContent(binx, biny);
}

double trigSFMuTrail(TFile* fileData, double pt, double eta)
{
   TH2F *h_dimu_trig_SF = (TH2F*) fileData->Get("mu_trail_leg_eta_v_pt_sf");
   TAxis *xaxis = h_dimu_trig_SF->GetXaxis();
   TAxis *yaxis = h_dimu_trig_SF->GetYaxis();
   Int_t binx = xaxis->FindBin(eta);
   Int_t biny = yaxis->FindBin(pt);
   return h_dimu_trig_SF->GetBinContent(binx, biny);
}

double trigSFElLead(TFile* fileData, double pt, double eta)
//double trigSFElLead(TFile* fileData, TH2F* hist2D, double pt, double eta)
{
   TH2F *h_diel_trig_SF = (TH2F*) fileData->Get("el_lead_leg_eta_v_pt_sf");
   TAxis *xaxis = h_diel_trig_SF->GetXaxis();
   TAxis *yaxis = h_diel_trig_SF->GetYaxis();
   Int_t binx = xaxis->FindBin(eta);
   Int_t biny = yaxis->FindBin(pt);
   return h_diel_trig_SF->GetBinContent(binx, biny);
}

double idSFMuTrk(TFile* fileData, double eta)
{
   TH1F *h_mu_trk_SF = (TH1F*) fileData->Get("muon_trk_sf");
   TAxis *xaxis = h_mu_trk_SF->GetXaxis();
   Int_t binx = xaxis->FindBin(eta);
   return h_mu_trk_SF->GetBinContent(binx);
}

double idSFMu(TFile* fileData, double pt, double eta)
{
   TH2F *h_mu_id_SF = (TH2F*) fileData->Get("muon_id_sf");
   TAxis *xaxis = h_mu_id_SF->GetXaxis();
   TAxis *yaxis = h_mu_id_SF->GetYaxis();
   Int_t binx = xaxis->FindBin(eta);
   Int_t biny = yaxis->FindBin(pt);
   return h_mu_id_SF->GetBinContent(binx, biny);
}

double sfMu(TFile* fileData, double pt, double eta)
{
   TH2F *h_mu_SF = (TH2F*) fileData->Get("sf_pt_vs_eta");
   TAxis *xaxis = h_mu_SF->GetXaxis();
   TAxis *yaxis = h_mu_SF->GetYaxis();
   Int_t binx = xaxis->FindBin(eta);
   Int_t biny = yaxis->FindBin(pt);
   return h_mu_SF->GetBinContent(binx, biny);
}

double idSFEl(TFile* fileData, TH2F* hist2D, double pt, double eta)
{
   //TH2F *h_el_id_SF = (TH2F*) fileData->Get("EGamma_SF2D");
   TAxis *xaxis = hist2D->GetXaxis();
   TAxis *yaxis = hist2D->GetYaxis();
   Int_t binx = xaxis->FindBin(eta);
   Int_t biny = yaxis->FindBin(pt);
   if(pt < 50.0 or pt > 500) return 1.0;
   else return hist2D->GetBinContent(binx, biny);
}

double idSFElMVA(TFile* fileData, TH2F* hist2D, double pt, double eta)
{
   //TH2F *h_el_idMVA_SF = (TH2F*) fileData->Get("EGamma_SF2D");
   TAxis *xaxis = hist2D->GetXaxis();
   TAxis *yaxis = hist2D->GetYaxis();
   Int_t binx = xaxis->FindBin(eta);
   Int_t biny = yaxis->FindBin(pt);
   return hist2D->GetBinContent(binx, biny);
}

double sfEl(TFile* fileData, TH2F* hist2D, double pt, double eta)
{
   //TH2F *h_el_SF = (TH2F*) fileData->Get("sf_pt_vs_eta");
   TAxis *xaxis = hist2D->GetXaxis();
   TAxis *yaxis = hist2D->GetYaxis();
   Int_t binx = xaxis->FindBin(eta);
   Int_t biny = yaxis->FindBin(pt);
   return hist2D->GetBinContent(binx, biny);
}
/*
double idSFEl(double pt, double eta)
{

  double sf = 1.0;
  if(eta > 0.0 and eta < 1.0)
  {
    if(pt > 25 and pt < 40)


  } 

}*/
