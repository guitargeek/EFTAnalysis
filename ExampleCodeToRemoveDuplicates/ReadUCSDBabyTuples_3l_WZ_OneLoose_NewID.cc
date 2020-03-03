#include "ReadUCSDBabyTuples_2018.h"
#include "Math/LorentzVector.h"
#include "Math/GenVector/LorentzVector.h"
#include "TROOT.h"

using std::string;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

int closestJet(const TLorentzVector& lep_p4, std::vector<TLorentzVector>& v_Jets, float dRmin, float maxAbsEta, TLorentzVector& jet_p4)
{
  bool returnValue = -1;
  for (unsigned int pfjidx=0; pfjidx<v_Jets.size(); ++pfjidx)
  {
    TLorentzVector tmp_jet_p4 = v_Jets.at(pfjidx);
    if (fabs(tmp_jet_p4.Eta()) < maxAbsEta)
    {
      float tmp_dRmin = tmp_jet_p4.DeltaR(lep_p4);
      if (tmp_dRmin < dRmin)
      {
        returnValue = pfjidx;
	dRmin = tmp_dRmin;
	jet_p4 = tmp_jet_p4;
      }
    }
  }
  return  returnValue;
}

bool CutHLT(std::vector<leptonInfo> v_leptons, int HLT_DoubleEl_temp, int HLT_MuEG_temp, int HLT_DoubleMu_temp)
{
  bool passTrigger = false;
  for (unsigned int i=0; i<v_leptons.size(); i++)
  {
    for (unsigned int j=0; j<v_leptons.size(); j++)
    {
      if (i==j) continue;
      // Check if any of the combination of leptons pass the trigger thresholds
      // Ele 23 12
      // El23 Mu8
      // Mu23 El12
      // Mu 17 8
      // The thresholds are rounded up to 25, 15, or 10
      if (abs(v_leptons.at(i).id) == 11 and abs(v_leptons.at(j).id) == 11)
        passTrigger |= (HLT_DoubleEl_temp and v_leptons.at(i).pt > 25 and v_leptons.at(j).pt > 15);
      else if (abs(v_leptons.at(i).id) == 13 and abs(v_leptons.at(j).id) == 11)
        passTrigger |= (HLT_MuEG_temp and v_leptons.at(i).pt > 25 and v_leptons.at(j).pt > 15);
      else if (abs(v_leptons.at(i).id) == 11 and abs(v_leptons.at(j).id) == 13)
        passTrigger |= (HLT_MuEG_temp and v_leptons.at(i).pt > 25 and v_leptons.at(j).pt > 10);
      else if (abs(v_leptons.at(i).id) == 13 and abs(v_leptons.at(j).id) == 13)
        passTrigger |= (HLT_DoubleMu_temp and v_leptons.at(i).pt > 20 and v_leptons.at(j).pt > 10);
    }
  }
  return passTrigger;
}

void is5leptonZandWtag(std::vector<leptonInfo> v_leptons, double met_pt, double met_phi, int &z_lep1, int &z_lep2, int &z_lep3, int &z_lep4, int &w_lep)
{
  double chi1_sq, chi2_sq, mT;
  z_lep1=z_lep2=z_lep3=z_lep4=w_lep=-999;
  chi1_sq=chi2_sq=0.0;
  double Mz = 91.1876;
  double Mw = 80.379;
  double pair1massDiff, pair2massDiff;
  pair1massDiff=pair2massDiff=0.0;
  double compare1 = 10;
  double compare2 = 10;
  for(unsigned int i=0; i<v_leptons.size(); i++)
  {
    for(unsigned int j=0; j<v_leptons.size(); j++)
    {
      if(i!=j)//make sure not checking same lepton
      {
        if(v_leptons.at(i).id*v_leptons.at(j).id==-121 or v_leptons.at(i).id*v_leptons.at(j).id==-169) //check opposite sign pair 
        {
          //chi1_sq = pow(((v_leptons.at(i).lep_lv+v_leptons.at(j).lep_lv).M() - Mz), 2);
          pair1massDiff = ((v_leptons.at(i).lep_lv+v_leptons.at(j).lep_lv).M() - Mz);
          for(unsigned int k=0; k<v_leptons.size(); k++)
          {
            for(unsigned int l=0; l<v_leptons.size(); l++)
            {
              if(j!=l and j!=k and i!=k and i!=l)//make sure not checking same lepton
              {
                if(v_leptons.at(k).id*v_leptons.at(l).id==-121 or v_leptons.at(k).id*v_leptons.at(l).id==-169) //check second opposite sign pair 
                {
                  //chi2_sq = pow(((v_leptons.at(k).lep_lv+v_leptons.at(l).lep_lv).M() - Mz), 2);
                  pair2massDiff = ((v_leptons.at(k).lep_lv+v_leptons.at(l).lep_lv).M() - Mz);
                  if(fabs(pair1massDiff) < compare1 and fabs(pair2massDiff) < compare2) 
                  {
                    compare1 = fabs(pair1massDiff);
                    compare2 = fabs(pair2massDiff);
                    z_lep1=i;
                    z_lep2=j;
                    z_lep3=k;
                    z_lep4=l;
                    for(unsigned int m=0; m<v_leptons.size(); m++)
                    {
                      if(i!=m and j!=m and k!=m and l!=m)//make sure not checking same lepton
                      {
                        mT = sqrt(2*met_pt*v_leptons.at(m).pt*(1.0 - cos(v_leptons.at(m).phi - met_phi)));
                        //w_lep=m;
                        if(mT > 50.0) w_lep=m;
                        //if(fabs(mT-Mw)<=20) w_lep=m;
                      }//w-tag
                    }//m-loop  
                  }//zz-tag
                }//second opposite sign pair
              }//lepton index check
            }//l-loop 
          }//k-loop
        }//first opposite sign pair
      }//lepton index check
    }//j-loop
  }//i-loop
}

void is4leptonZandWtag(std::vector<leptonInfo> v_leptons, double met_pt, double met_phi, int &z_lep1, int &z_lep2, int &w_lep)
{
  double mT = 0.0;
  z_lep1=z_lep2=w_lep=-999;
  double Mz = 91.1876;
  double Mw = 80.379;
  double pair1massDiff, pair2massDiff;
  pair1massDiff=pair2massDiff=0.0;
  double compare1 = 10;
  for(unsigned int i=0; i<v_leptons.size(); i++)
  {
    for(unsigned int j=0; j<v_leptons.size(); j++)
    {
      if(i!=j)//make sure not checking same lepton
      {
        if(v_leptons.at(i).id*v_leptons.at(j).id==-121 or v_leptons.at(i).id*v_leptons.at(j).id==-169) //check opposite sign pair 
        {
          pair1massDiff = ((v_leptons.at(i).lep_lv+v_leptons.at(j).lep_lv).M() - Mz);
          if(fabs(pair1massDiff) < compare1)
          {
            compare1 = fabs(pair1massDiff);
            z_lep1=i;
            z_lep2=j;
            for(unsigned int k=0; k<v_leptons.size(); k++)
            {
              if(i!=k and j!=k)//make sure not checking same lepton
              {
                mT = sqrt(2*met_pt*v_leptons.at(k).pt*(1.0 - cos(v_leptons.at(k).phi - met_phi)));
                w_lep=k;
                //if(mT > 50.0) w_lep=k;
                //if(fabs(mT-Mw)<=20) w_lep=m;
              }//w-tag
            }//k-loop 
          }//first opposite sign pair
        }//lepton index check
      }//not checking same lepton 
    }//j-loop
  }//i-loop
}

bool elecIsoCut(float eta, float pt, float iso)
{
  if(fabs(eta) <= 1.479)
  {
    if(iso >= 0.198 + 0.506/pt) return false;
  }
  else if ((fabs(eta) > 1.479) && (fabs(eta) < 2.5))
  {
    if (iso >= 0.203 + 0.963/pt) return false;
  }
  else return false;
  return true;
}

float wMassConstraintPlus(float E_lep, float plep_x, float plep_y, float plep_z, float pv_T, float pv_x, float pv_y, float Mw, float Ml)
{
  float A = 4*(E_lep*E_lep - plep_z*plep_z);
  float a = Mw*Mw - Ml*Ml + 2*(plep_x*pv_x + plep_y*pv_y);
  float b = -4*a*plep_z;
  float C = 4*(E_lep*E_lep)*(pv_T*pv_T) - a*a;
  if(b*b - 4*A*C < 0.0) return (1.0/(2*A))*(-b);
  return (1.0/(2*A))*(-b + sqrt(b*b - 4*A*C));

}

float wMassConstraintMinus(float E_lep, float plep_x, float plep_y, float plep_z, float pv_T, float pv_x, float pv_y, float Mw, float Ml)
{
  float A = 4*(E_lep*E_lep - plep_z*plep_z);
  float a = Mw*Mw - Ml*Ml + 2*(plep_x*pv_x + plep_y*pv_y);
  float b = -4*a*plep_z;
  float C = 4*(E_lep*E_lep)*(pv_T*pv_T) - a*a;
  if(b*b - 4*A*C < 0.0) return (1.0/(2*A))*(-b);
  else return (1.0/(2*A))*(-b - sqrt(b*b - 4*A*C));

}

int ReadUCSDBabyTuples_3l_WZ_OneLoose_NewID(std::string infile, std::string treeStr, std::string Sample, std::string Trigger="None")
{

  std::string inputfilename=(infile+".root").c_str();
  TFile *inputFile = new TFile((inputfilename).c_str());
  TChain *tree=new TChain(treeStr.c_str());
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  TFile *sfElRecoFile = new TFile("egammaEffi.txt_EGM2D_updatedAll.root");
  TH2F  *idSFElRecoHist = (TH2F*) sfElRecoFile->Get("EGamma_SF2D");
  TFile *sfElVetoFile = new TFile("2018_ElectronWPVeto_Fall17V2.root");
  TH2F  *idSFElVetoHist = (TH2F*) sfElVetoFile->Get("EGamma_SF2D");

  TFile *sfMuRecoFile = new TFile("EfficiencyStudies_2018_rootfiles_RunABCD_SF_ID.root");
  TH2F  *sfMuRecoHist = (TH2F*) sfMuRecoFile->Get("NUM_MediumID_DEN_TrackerMuons_pt_abseta");
  TFile *sfMuRecoLowPtFile = new TFile("EfficiencyStudies_2018_rootfiles_lowpt_RunABCD_SF_ID.root");
  TH2F  *sfMuRecoLowPtHist = (TH2F*) sfMuRecoLowPtFile->Get("NUM_MediumID_DEN_genTracks_pt_abseta");
  TFile *sfMuIsoFile = new TFile("EfficiencyStudies_2018_rootfiles_RunABCD_SF_ISO.root");
  TH2F  *sfMuIsoHist = (TH2F*) sfMuIsoFile->Get("NUM_LooseRelIso_DEN_MediumID_pt_abseta");

  bool debug=false;

  Int_t           run;
  Int_t           lumi;
  ULong64_t       evt;
  Int_t           isData;
  Float_t         evt_scale1fb;
  Float_t         genps_weight;
  Float_t         xsec_br;
  Int_t           evt_passgoodrunlist;
  TString         *CMS4path;
  Int_t           CMS4index;
  Float_t         weight_fr_r1_f1;
  Float_t         weight_fr_r1_f2;
  Float_t         weight_fr_r1_f0p5;
  Float_t         weight_fr_r2_f1;
  Float_t         weight_fr_r2_f2;
  Float_t         weight_fr_r2_f0p5;
  Float_t         weight_fr_r0p5_f1;
  Float_t         weight_fr_r0p5_f2;
  Float_t         weight_fr_r0p5_f0p5;
  Float_t         weight_pdf_up;
  Float_t         weight_pdf_down;
  Float_t         weight_alphas_down;
  Float_t         weight_alphas_up;
  Int_t           HLT_DoubleMu;
  Int_t           HLT_DoubleEl;
  Int_t           HLT_MuEG;
  Int_t           pass_duplicate_ee_em_mm;
  Int_t           pass_duplicate_mm_em_ee;
  Float_t         gen_ht;
  std::vector<LorentzVector>  *gen_V_p4;
  std::vector<float>   *gen_V_pt;
  std::vector<float>   *gen_V_eta;
  std::vector<float>   *gen_V_phi;
  std::vector<float>   *gen_V_mass;
  std::vector<int>     *gen_V_id;
  std::vector<LorentzVector>   *gen_lep_p4;
  std::vector<float>   *gen_lep_pt;
  std::vector<float>   *gen_lep_eta;
  std::vector<float>   *gen_lep_phi;
  std::vector<float>   *gen_lep_mass;
  std::vector<int>     *gen_lep_id;
  Int_t           VHchannel;
  Int_t           Higgschannel;
  Int_t           firstgoodvertex;
  Int_t           nvtx;
  Int_t           nTrueInt;
  std::vector<LorentzVector>  *lep_p4;
  std::vector<float>   *lep_pt;
  std::vector<float>   *lep_eta;
  std::vector<float>   *lep_phi;
  std::vector<float>   *lep_energy;
  std::vector<float>   *lep_mva;
  std::vector<float>   *lep_relIso04DB;
  std::vector<float>   *lep_relIso03EA;
  std::vector<float>   *lep_relIso03EAwLep;
  std::vector<float>   *lep_ip3d;
  std::vector<float>   *lep_sip3d;
  std::vector<float>   *lep_dxy;
  std::vector<float>   *lep_dz;
  std::vector<int>     *lep_mc_motherid;
  std::vector<int>     *lep_mc_id;
  std::vector<int>     *lep_motherIdv2;
  std::vector<int>     *lep_idx;
  std::vector<int>     *lep_id;
  std::vector<int>     *lep_isTightPOG;
  std::vector<int>     *lep_isMediumPOG;
  std::vector<int>     *lep_isCutBasedNoIsoVetoPOG;
  std::vector<int>     *lep_isCutBasedNoIsoLoosePOG;
  std::vector<int>     *lep_isCutBasedNoIsoMediumPOG;
  std::vector<int>     *lep_isCutBasedNoIsoTightPOG;
  std::vector<int>     *lep_isCutBasedIsoVetoPOG;
  std::vector<int>     *lep_isCutBasedIsoLoosePOG;
  std::vector<int>     *lep_isCutBasedIsoMediumPOG;
  std::vector<int>     *lep_isCutBasedIsoTightPOG;  
  Float_t         met_pt;
  Float_t         met_phi;
  Float_t         met_up_pt;
  Float_t         met_up_phi;
  Float_t         met_dn_pt;
  Float_t         met_dn_phi;
  Float_t         met_gen_pt;
  Float_t         met_gen_phi;
  std::vector<LorentzVector>  *jets_p4;
  std::vector<float>   *jets_pt;
  std::vector<float>   *jets_eta;
  std::vector<float>   *jets_phi;
  std::vector<float>   *jets_mass;
  std::vector<LorentzVector>  *jets_cen_p4;
  std::vector<float>   *jets_cen_pt;
  std::vector<float>   *jets_cen_eta;
  std::vector<float>   *jets_cen_phi;
  std::vector<float>   *jets_cen_mass;
  std::vector<int>     *lep_isVVVVeto;
  std::vector<int>     *lep_isMVAwp90IsoPOG;
  Int_t           nj;
  Int_t           nb;
  Int_t           nbmed;
  Float_t         ht;
  Int_t           nj_cen;
  Float_t         weight_btagsf;
  Float_t         weight_btagsf_heavy_DN;
  Float_t         weight_btagsf_heavy_UP;
  Float_t         weight_btagsf_light_DN;
  Float_t         weight_btagsf_light_UP;

  CMS4path     = 0; 
  gen_V_p4     = 0;
  gen_V_pt     = 0;
  gen_V_eta    = 0;
  gen_V_phi    = 0;
  gen_V_mass   = 0;
  gen_V_id     = 0;
  gen_lep_p4  = 0;
  gen_lep_pt   = 0;
  gen_lep_eta  = 0;
  gen_lep_phi  = 0;
  gen_lep_mass = 0;
  gen_lep_id   = 0;
  lep_p4       = 0;
  lep_pt       = 0;
  lep_eta      = 0;
  lep_phi      = 0;
  lep_energy   = 0;
  lep_mva      = 0;
  lep_relIso04DB = 0;
  lep_relIso03EA = 0;
  lep_relIso03EAwLep = 0;
  lep_ip3d = 0;
  lep_sip3d = 0;
  lep_dxy = 0;
  lep_dz  = 0;
  lep_mc_id = 0;
  lep_motherIdv2 = 0;
  lep_idx = 0;
  lep_id  = 0;
  lep_mc_motherid = 0;
  lep_isTightPOG = 0;
  lep_isMediumPOG = 0;
  lep_isCutBasedNoIsoVetoPOG = 0;
  lep_isCutBasedNoIsoLoosePOG = 0;
  lep_isCutBasedNoIsoMediumPOG = 0;
  lep_isCutBasedNoIsoTightPOG = 0;
  lep_isCutBasedIsoVetoPOG = 0;
  lep_isCutBasedIsoLoosePOG = 0;
  lep_isCutBasedIsoMediumPOG = 0;
  lep_isCutBasedIsoTightPOG = 0;
  lep_isVVVVeto = 0;
  lep_isMVAwp90IsoPOG = 0;
  jets_p4 = 0;
  jets_pt = 0;
  jets_eta = 0;
  jets_phi = 0;
  jets_mass = 0;
  jets_cen_p4 = 0;
  jets_cen_pt = 0;
  jets_cen_eta = 0;
  jets_cen_phi = 0;
  jets_cen_mass = 0;

  tree->SetBranchAddress("run", &(run));
  tree->SetBranchAddress("lumi", &(lumi));
  tree->SetBranchAddress("evt", &(evt));
  tree->SetBranchAddress("isData", &(isData));
  tree->SetBranchAddress("evt_scale1fb", &(evt_scale1fb));
  tree->SetBranchAddress("genps_weight", &(genps_weight));
  tree->SetBranchAddress("xsec_br", &(xsec_br));
  tree->SetBranchAddress("evt_passgoodrunlist", &(evt_passgoodrunlist));
  tree->SetBranchAddress("CMS4path", &(CMS4path));
  tree->SetBranchAddress("CMS4index", &(CMS4index));
  tree->SetBranchAddress("weight_fr_r1_f1", &(weight_fr_r1_f1));
  tree->SetBranchAddress("weight_fr_r1_f2", &(weight_fr_r1_f2));
  tree->SetBranchAddress("weight_fr_r1_f0p5", &(weight_fr_r1_f0p5));
  tree->SetBranchAddress("weight_fr_r2_f1", &(weight_fr_r2_f1));
  tree->SetBranchAddress("weight_fr_r2_f2", &(weight_fr_r2_f2));
  tree->SetBranchAddress("weight_fr_r2_f0p5", &(weight_fr_r2_f0p5));
  tree->SetBranchAddress("weight_fr_r0p5_f1", &(weight_fr_r0p5_f1));
  tree->SetBranchAddress("weight_fr_r0p5_f2", &(weight_fr_r0p5_f2));
  tree->SetBranchAddress("weight_fr_r0p5_f0p5", &(weight_fr_r0p5_f0p5));
  tree->SetBranchAddress("weight_pdf_up", &(weight_pdf_up));
  tree->SetBranchAddress("weight_pdf_down", &(weight_pdf_down));
  tree->SetBranchAddress("weight_alphas_down", &(weight_alphas_down));
  tree->SetBranchAddress("weight_alphas_up", &(weight_alphas_up));
  tree->SetBranchAddress("HLT_DoubleMu", &(HLT_DoubleMu));
  tree->SetBranchAddress("HLT_DoubleEl", &(HLT_DoubleEl));
  tree->SetBranchAddress("HLT_MuEG", &(HLT_MuEG));
  tree->SetBranchAddress("pass_duplicate_ee_em_mm", &(pass_duplicate_ee_em_mm));
  tree->SetBranchAddress("pass_duplicate_mm_em_ee", &(pass_duplicate_mm_em_ee));
  tree->SetBranchAddress("gen_ht", &(gen_ht));
  tree->SetBranchAddress("gen_V_p4", &(gen_V_p4));
  tree->SetBranchAddress("gen_V_pt", &(gen_V_pt));
  tree->SetBranchAddress("gen_V_eta", &(gen_V_eta));
  tree->SetBranchAddress("gen_V_phi", &(gen_V_phi));
  tree->SetBranchAddress("gen_V_mass", &(gen_V_mass));
  tree->SetBranchAddress("gen_V_id", &(gen_V_id));
  tree->SetBranchAddress("gen_lep_p4", &(gen_lep_p4));
  tree->SetBranchAddress("gen_lep_pt", &(gen_lep_pt));
  tree->SetBranchAddress("gen_lep_eta", &(gen_lep_eta));
  tree->SetBranchAddress("gen_lep_phi", &(gen_lep_phi));
  tree->SetBranchAddress("gen_lep_mass", &(gen_lep_mass));
  tree->SetBranchAddress("gen_lep_id", &(gen_lep_id));
  tree->SetBranchAddress("VHchannel", &(VHchannel));
  tree->SetBranchAddress("Higgschannel", &(Higgschannel));
  tree->SetBranchAddress("firstgoodvertex", &(firstgoodvertex));
  tree->SetBranchAddress("nvtx", &(nvtx));
  tree->SetBranchAddress("nTrueInt", &(nTrueInt));
  tree->SetBranchAddress("lep_p4", &(lep_p4));
  tree->SetBranchAddress("lep_pt", &(lep_pt));
  tree->SetBranchAddress("lep_eta", &(lep_eta));
  tree->SetBranchAddress("lep_phi", &(lep_phi));
  tree->SetBranchAddress("lep_energy", &(lep_energy));
  tree->SetBranchAddress("lep_mva", &(lep_mva));
  tree->SetBranchAddress("lep_relIso04DB", &(lep_relIso04DB));
  tree->SetBranchAddress("lep_relIso03EA", &(lep_relIso03EA));
  tree->SetBranchAddress("lep_relIso03EAwLep", &(lep_relIso03EAwLep));
  tree->SetBranchAddress("lep_ip3d", &(lep_ip3d));
  tree->SetBranchAddress("lep_sip3d", &(lep_sip3d));
  tree->SetBranchAddress("lep_dxy", &(lep_dxy));
  tree->SetBranchAddress("lep_dz", &(lep_dz));
  tree->SetBranchAddress("lep_mc_id", &(lep_mc_id));
  tree->SetBranchAddress("lep_motherIdv2", &(lep_motherIdv2));
  tree->SetBranchAddress("lep_idx", &(lep_idx));
  tree->SetBranchAddress("lep_id", &(lep_id));
  tree->SetBranchAddress("lep_mc_motherid", &(lep_mc_motherid));
  tree->SetBranchAddress("lep_isTightPOG", &(lep_isTightPOG));
  tree->SetBranchAddress("lep_isMediumPOG", &(lep_isMediumPOG));
  tree->SetBranchAddress("lep_isCutBasedNoIsoVetoPOG", &(lep_isCutBasedNoIsoVetoPOG));
  tree->SetBranchAddress("lep_isCutBasedNoIsoLoosePOG", &(lep_isCutBasedNoIsoLoosePOG));
  tree->SetBranchAddress("lep_isCutBasedNoIsoMediumPOG", &(lep_isCutBasedNoIsoMediumPOG));
  tree->SetBranchAddress("lep_isCutBasedNoIsoTightPOG", &(lep_isCutBasedNoIsoTightPOG));
  tree->SetBranchAddress("lep_isCutBasedIsoVetoPOG", &(lep_isCutBasedIsoVetoPOG));
  tree->SetBranchAddress("lep_isCutBasedIsoLoosePOG", &(lep_isCutBasedIsoLoosePOG));
  tree->SetBranchAddress("lep_isCutBasedIsoMediumPOG", &(lep_isCutBasedIsoMediumPOG));
  tree->SetBranchAddress("lep_isCutBasedIsoTightPOG", &(lep_isCutBasedIsoTightPOG));
  tree->SetBranchAddress("lep_isVVVVeto", &(lep_isVVVVeto));
  tree->SetBranchAddress("lep_isMVAwp90IsoPOG", &(lep_isMVAwp90IsoPOG));
  tree->SetBranchAddress("met_pt", &(met_pt));
  tree->SetBranchAddress("met_phi", &(met_phi));
  tree->SetBranchAddress("met_up_pt", &(met_up_pt));
  tree->SetBranchAddress("met_up_phi", &(met_up_phi));
  tree->SetBranchAddress("met_dn_pt", &(met_dn_pt));
  tree->SetBranchAddress("met_dn_phi", &(met_dn_phi));
  tree->SetBranchAddress("met_gen_pt", &(met_gen_pt));
  tree->SetBranchAddress("met_gen_phi", &(met_gen_phi));
  tree->SetBranchAddress("jets_p4", &(jets_p4));
  tree->SetBranchAddress("jets_pt", &(jets_pt));
  tree->SetBranchAddress("jets_eta", &(jets_eta));
  tree->SetBranchAddress("jets_phi", &(jets_phi));
  tree->SetBranchAddress("jets_mass", &(jets_mass));
  tree->SetBranchAddress("jets_cen_p4", &(jets_cen_p4));
  tree->SetBranchAddress("jets_cen_pt", &(jets_cen_pt));
  tree->SetBranchAddress("jets_cen_eta", &(jets_cen_eta));
  tree->SetBranchAddress("jets_cen_phi", &(jets_cen_phi));
  tree->SetBranchAddress("jets_cen_mass", &(jets_cen_mass));
  tree->SetBranchAddress("nj", &(nj));
  tree->SetBranchAddress("nb", &(nb));
  tree->SetBranchAddress("nbmed", &(nbmed));
  tree->SetBranchAddress("ht", &(ht));
  tree->SetBranchAddress("nj_cen", &(nj_cen));
  tree->SetBranchAddress("weight_btagsf", &(weight_btagsf));
  tree->SetBranchAddress("weight_btagsf_heavy_DN", &(weight_btagsf_heavy_DN));
  tree->SetBranchAddress("weight_btagsf_heavy_UP", &(weight_btagsf_heavy_UP));
  tree->SetBranchAddress("weight_btagsf_light_DN", &(weight_btagsf_light_DN));
  tree->SetBranchAddress("weight_btagsf_light_UP", &(weight_btagsf_light_UP));

  TH1D *h_nLeptons = new TH1D("h_nLeptons", "h_nLeptons", 10, -0.5, 9.5);h_nLeptons->Sumw2();

  TH2D *h_ee_mm_ElMu = new TH2D("h_ee_mm_ElMu", "h_ee_mm_ElMu", 300, 0.0, 300.0, 300, 0.0, 300.0); h_ee_mm_ElMu->Sumw2();
  TH2D *h_ee_ee_4El = new TH2D("h_ee_ee_4El", "h_ee_ee_4El", 300, 0.0, 300.0, 300, 0.0, 300.0); h_ee_ee_4El->Sumw2();
  TH2D *h_mm_mm_4Mu = new TH2D("h_mm_mm_4Mu", "h_mm_mm_4Mu", 300, 0.0, 300.0, 300, 0.0, 300.0); h_mm_mm_4Mu->Sumw2();

  HistCollection Cut1_3Mu1El;//no cuts
  initializeHistCollection(Cut1_3Mu1El, "Cut1_3Mu1El");

  HistCollection Cut2_3Mu1El;//b-tag
  initializeHistCollection(Cut2_3Mu1El, "Cut2_3Mu1El");

  HistCollection Cut3_3Mu1El;//z-tag
  initializeHistCollection(Cut3_3Mu1El, "Cut3_3Mu1El");

  HistCollection Cut4_3Mu1El;//MET > 40
  initializeHistCollection(Cut4_3Mu1El, "Cut4_3Mu1El");

  HistCollection Cut1_2Mu2El_LooseMu;//no cuts
  initializeHistCollection(Cut1_2Mu2El_LooseMu, "Cut1_2Mu2El_LooseMu");

  HistCollection Cut2_2Mu2El_LooseMu;//b-tag
  initializeHistCollection(Cut2_2Mu2El_LooseMu, "Cut2_2Mu2El_LooseMu");

  HistCollection Cut3_2MuPlus2El_LooseMu;//z-tag
  initializeHistCollection(Cut3_2MuPlus2El_LooseMu, "Cut3_2MuPlus2El_LooseMu");

  HistCollection Cut4_2MuPlus2El_LooseMu;//MET > 40
  initializeHistCollection(Cut4_2MuPlus2El_LooseMu, "Cut4_2MuPlus2El_LooseMu"); 

  HistCollection Cut3_2MuMinus2El_LooseMu;//z-tag
  initializeHistCollection(Cut3_2MuMinus2El_LooseMu, "Cut3_2MuMinus2El_LooseMu");

  HistCollection Cut4_2MuMinus2El_LooseMu;//MET > 40
  initializeHistCollection(Cut4_2MuMinus2El_LooseMu, "Cut4_2MuMinus2El_LooseMu");

  HistCollection Cut1_2Mu2El_LooseEl;//no cuts
  initializeHistCollection(Cut1_2Mu2El_LooseEl, "Cut1_2Mu2El_LooseEl");

  HistCollection Cut2_2Mu2El_LooseEl;//b-tag
  initializeHistCollection(Cut2_2Mu2El_LooseEl, "Cut2_2Mu2El_LooseEl");

  HistCollection Cut3_2Mu2ElPlus_LooseEl;//z-tag
  initializeHistCollection(Cut3_2Mu2ElPlus_LooseEl, "Cut3_2Mu2ElPlus_LooseEl");

  HistCollection Cut4_2Mu2ElPlus_LooseEl;//MET > 40
  initializeHistCollection(Cut4_2Mu2ElPlus_LooseEl, "Cut4_2Mu2ElPlus_LooseEl");

  HistCollection Cut3_2Mu2ElMinus_LooseEl;//z-tag
  initializeHistCollection(Cut3_2Mu2ElMinus_LooseEl, "Cut3_2Mu2ElMinus_LooseEl");

  HistCollection Cut4_2Mu2ElMinus_LooseEl;//MET > 40
  initializeHistCollection(Cut4_2Mu2ElMinus_LooseEl, "Cut4_2Mu2ElMinus_LooseEl");


  TH1D *h_TotalEvents_4Mu = new TH1D("h_TotalEvents_4Mu", "h_TotalEvents_4Mu", 15, -0.5, 14.5); h_TotalEvents_4Mu->Sumw2();
  TH1D *h_TotalEvents_1El2Mu = new TH1D("h_TotalEvents_1El2Mu", "h_TotalEvents_1El2Mu", 15, -0.5, 14.5); h_TotalEvents_1El2Mu->Sumw2();
  TH1D *h_TotalEvents_2El1Mu = new TH1D("h_TotalEvents_2El1Mu", "h_TotalEvents_2El1Mu", 15, -0.5, 14.5); h_TotalEvents_2El1Mu->Sumw2();
  TH1D *h_TotalEvents_4El = new TH1D("h_TotalEvents_4El", "h_TotalEvents_4El", 15, -0.5, 14.5); h_TotalEvents_4El->Sumw2();

  TH1D *h_TotalEvents_6l = new TH1D("h_TotalEvents_6l", "h_TotalEvents_6l", 15, -0.5, 14.5); h_TotalEvents_6l->Sumw2();

  TH1D *h_zmass_1 = new TH1D("h_zmass_1", "h_zmass_1", 500, 0.0, 500.0);h_zmass_1->Sumw2();
  TH1D *h_zmass_2 = new TH1D("h_zmass_2", "h_zmass_2", 500, 0.0, 500.0);h_zmass_2->Sumw2();
  TH1D *h_mT = new TH1D("h_mT", "h_mT", 500, 0.0, 500.0);h_mT->Sumw2();
  TH1D *h_nJets = new TH1D("h_nJets", "h_nJets", 15, -0.5, 14.5); h_nJets->Sumw2();

  int nEvents=tree->GetEntries();
  std::cout << nEvents << std::endl;
  //nEvents=10;
  int lep_ZCand_idx1 = 0;
  unsigned int lep_ZCand_idx2 = 0;
  int lep_WCand_idx1 = 0;
  int lep_2ndZCand_idx1 = 0;
  int lep_2ndZCand_idx2 = 0;
  vector <long long> checkDuplicates;
  int nDup = 0;
  long long RUNPREF = 1000 * 1000;
  RUNPREF *= 1000 * 1000;
  int n_5l, n_pos_5l, n_neg_5l;
  n_5l=n_pos_5l=n_neg_5l=0;
  int mother_heavy_flavor, mother_pi0, mother_pi, mother_lep, mother_other;
  mother_heavy_flavor=mother_pi0=mother_pi=mother_lep=mother_other=0;
  double n_4mu=0;
  //nEvents=10;
  int n_1el, n_3mu;
  n_1el=n_3mu=0;
  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);  
   
    double weight = evt_scale1fb*getTruePUw2018(nTrueInt);
    if(evt_passgoodrunlist==0) continue;
    if(firstgoodvertex!=0) continue;

    std::vector<TLorentzVector> v_selectedJets;
    for(unsigned int iselJet=0; iselJet<jets_p4->size(); ++iselJet)
    {
      TLorentzVector Jet;

      if(fabs(jets_p4->at(iselJet).Eta())<5.0 and jets_p4->at(iselJet).Pt()>30.0)
      {
        Jet.SetPtEtaPhiE(jets_p4->at(iselJet).Pt(), jets_p4->at(iselJet).Eta(), jets_p4->at(iselJet).Phi(), jets_p4->at(iselJet).E());
        v_selectedJets.push_back(Jet);
      }
    }

    std::vector<leptonInfo> v_leptons;
    v_leptons.clear();

    std::vector<leptonInfo> v_muons;
    v_muons.clear();

    for (unsigned int imuon=0; imuon<lep_p4->size(); imuon++)
    {
       leptonInfo muon;
       muon.pt = lep_pt->at(imuon); 
       muon.phi = lep_phi->at(imuon);
       muon.eta = lep_eta->at(imuon);
       muon.sip3d = lep_sip3d->at(imuon);
       muon.dxy  = lep_dxy->at(imuon);
       muon.dz = lep_dz->at(imuon);
       muon.id = lep_id->at(imuon);
       muon.lep_lv.SetPtEtaPhiM(lep_pt->at(imuon), lep_eta->at(imuon), lep_phi->at(imuon), MUON_MASS);
       muon.iso = lep_relIso04DB->at(imuon);
       //muon.sf = recoSFMu(sfMuRecoFile, sfMuRecoHist, lep_pt->at(imuon), lep_eta->at(imuon))*recoLowPtSFMu(sfMuRecoLowPtFile, sfMuRecoLowPtHist, lep_pt->at(imuon), lep_eta->at(imuon))*isoSFMu(sfMuIsoFile, sfMuIsoHist, lep_pt->at(imuon), lep_eta->at(imuon));
       muon.charge = lep_id->at(imuon)/13;
       muon.sf = 1.0;
       //if(abs(lep_id->at(imuon))==13 and muon.pt>10.0 and fabs(lep_eta->at(imuon)) < 2.4 and fabs(lep_relIso04DB->at(imuon)) > 0.15 and fabs(lep_sip3d->at(imuon)) < 4 and lep_isMediumPOG->at(imuon)==1)
       //if(abs(lep_id->at(imuon))==13 and muon.pt>10.0 and fabs(lep_eta->at(imuon)) < 2.4 and fabs(lep_relIso04DB->at(imuon)) > 0.25 and fabs(lep_sip3d->at(imuon)) < 4 and lep_isMediumPOG->at(imuon)==1) 
       if(abs(lep_id->at(imuon))==13 and muon.pt>10.0 and fabs(lep_eta->at(imuon)) < 2.4 and lep_isVVVVeto->at(imuon)==1 and fabs(lep_sip3d->at(imuon)) < 4 and lep_isMediumPOG->at(imuon)==1 and fabs(lep_relIso04DB->at(imuon)) < 0.25)// and fabs(lep_sip3d->at(imuon)) < 4 and lep_isMediumPOG->at(imuon)==1)
       //if(abs(lep_id->at(imuon))==13 and muon.pt>10.0 and fabs(lep_eta->at(imuon)) < 2.4 and fabs(lep_sip3d->at(imuon)) < 4 and lep_isMediumPOG->at(imuon)==1) 
       //if(abs(lep_id->at(imuon))==13 and lep_isVVVVeto->at(imuon)==1 and (fabs(lep_sip3d->at(imuon)) < 4))
       //if((lep_isVVVVeto->at(imuon)==1) and (fabs(lep_sip3d->at(imuon)) < 4))// and (not (lep_isMVAwp90IsoPOG->at(imuon)==1)))
       {
         //std::cout << "fabs(lep_relIso04DB->at(imuon)) = " << fabs(lep_relIso04DB->at(imuon)) << std::endl;
         v_muons.push_back(muon);
         v_leptons.push_back(muon);
       }
       //if(abs(lep_id->at(imuon))==13 and muon.pt<15.0 and fabs(lep_eta->at(imuon)) < 2.4 and fabs(lep_sip3d->at(imuon)) < 4 and lep_isMediumPOG->at(imuon)==1) std::cout << "Low pt muon found" << std::endl; 
    }
    std::sort (v_muons.begin(), v_muons.end(), sortLeptonsInDescendingpT);

    std::vector<leptonInfo> v_electrons;
    v_electrons.clear();    

    for (unsigned int ielectron=0; ielectron<lep_p4->size(); ielectron++)
    {
       leptonInfo electron;
       electron.pt = lep_pt->at(ielectron);
       electron.phi = lep_phi->at(ielectron);
       electron.eta = lep_eta->at(ielectron);
       electron.sip3d = lep_sip3d->at(ielectron);
       electron.dxy  = lep_dxy->at(ielectron);
       electron.dz = lep_dz->at(ielectron);
       electron.mva = lep_mva->at(ielectron);
       electron.id = lep_id->at(ielectron);
       electron.lep_lv.SetPtEtaPhiM(lep_pt->at(ielectron), lep_eta->at(ielectron), lep_phi->at(ielectron), ELECTRON_MASS);
       electron.iso = lep_relIso03EA->at(ielectron);
       electron.charge = lep_id->at(ielectron)/11; 
       //electron.sf = recoSFEl(sfElRecoFile, idSFElRecoHist, lep_pt->at(ielectron), lep_eta->at(ielectron))*vetoSFEl(sfElVetoFile, idSFElVetoHist, lep_pt->at(ielectron), lep_eta->at(ielectron));
       electron.sf = 1.0;
       //if(abs(lep_id->at(ielectron))==11 and electron.pt>10.0 and fabs(lep_eta->at(ielectron)) < 2.5 and (elecIsoCut(lep_eta->at(ielectron), lep_pt->at(ielectron), lep_relIso03EA->at(ielectron))==false))
       //if(abs(lep_id->at(ielectron))==11 and electron.pt>10.0 and fabs(lep_eta->at(ielectron)) < 2.5 and lep_isCutBasedNoIsoVetoPOG->at(ielectron)==1 and fabs(lep_sip3d->at(ielectron)) < 4)  
       //if(abs(lep_id->at(ielectron))==11 and electron.pt>10.0 and fabs(lep_eta->at(ielectron)) < 2.5 and lep_isCutBasedIsoVetoPOG->at(ielectron)==1 and fabs(lep_sip3d->at(ielectron)) < 4) 
       //if(abs(lep_id->at(ielectron))==11 and electron.pt>10.0 and fabs(lep_eta->at(ielectron)) < 2.5)
       //if(abs(lep_id->at(ielectron))==11 and lep_isVVVVeto->at(ielectron)==1 and (fabs(lep_sip3d->at(ielectron)) < 4) and (not (lep_isMVAwp90IsoPOG->at(ielectron)==1)) and fabs(lep_sip3d->at(ielectron)) < 4)
       //if(abs(lep_id->at(ielectron))==11 and fabs(lep_relIso03EAwLep->at(ielectron)) > 0.2 and (fabs(lep_sip3d->at(ielectron)) < 4))
       //if(abs(lep_id->at(ielectron))==11 and lep_isVVVVeto->at(ielectron)==1 and fabs(lep_relIso03EAwLep->at(ielectron)) > 0.2)
       if(abs(lep_id->at(ielectron))==11 and lep_isVVVVeto->at(ielectron)==1 and (fabs(lep_sip3d->at(ielectron)) < 4) and fabs(lep_relIso03EAwLep->at(ielectron)) > 0.2)
       //if((lep_isVVVVeto->at(ielectron)==1) and (fabs(lep_sip3d->at(ielectron)) < 4) and (not (lep_isMVAwp90IsoPOG->at(ielectron)==1)))
       {
         //std::cout << "fabs(lep_relIso03EAwLep->at(ielectron)) = " << fabs(lep_relIso03EAwLep->at(ielectron)) << std::endl;
         v_electrons.push_back(electron);
         v_leptons.push_back(electron);
       }
       //if(fabs(lep_eta->at(ielectron)) < 2.5 and lep_isCutBasedNoIsoVetoPOG->at(ielectron)==1 and fabs(lep_sip3d->at(ielectron)) < 4 and electron.pt<15.0) std::cout << "Low pt electron found" << std::endl; 
    }
    std::sort (v_electrons.begin(), v_electrons.end(), sortLeptonsInDescendingpT);
    
    std::sort (v_leptons.begin(), v_leptons.end(), sortLeptonsInDescendingpT);

    unsigned int nv_leptons = v_muons.size() + v_electrons.size();
    if(nv_leptons != v_leptons.size()) std::cout << "size mismatch" << std::endl;   

    if(not CutHLT(v_leptons, HLT_DoubleEl, HLT_MuEG, HLT_DoubleMu)) continue;
 
    if(Sample=="MC") h_nLeptons->Fill(nv_leptons, weight);
    else if(Sample=="Data") 
    {
      long long dupCheck = run*RUNPREF + evt;
      bool bDuplicate = false;
      for (unsigned int uid = 0; uid < checkDuplicates.size(); uid++)
      {
        if (checkDuplicates[uid] == dupCheck)
        {
          cout<<dupCheck<<endl;
          bDuplicate = true;
          nDup++;
          break;
        }
      }
      if (bDuplicate) continue;
      else checkDuplicates.push_back(dupCheck);
      h_nLeptons->Fill(nv_leptons, 1.0); 
    }

    if(v_electrons.size()==1) n_1el++;
    if(v_muons.size()==3) n_3mu++;


    //if(v_electrons.size()==1)
    if(v_muons.size()==3 and v_electrons.size()==1)
    //if(v_muons.size()==3 and v_electrons.size()==1 and v_muons.at(0).iso < 0.25 and v_muons.at(1).iso < 0.25 and v_muons.at(2).iso < 0.25 and (elecIsoCut(v_electrons.at(0).eta, v_electrons.at(0).pt, v_electrons.at(0).iso)==false))
    {
      if(Sample=="MC") weight *= v_muons.at(0).sf*v_muons.at(1).sf*v_muons.at(2).sf*v_electrons.at(0).sf; 
      else weight = 1.0;
      //std::cout << "weight = " << weight << std::endl;
      //std::cout << " sf = " << v_muons.at(0).sf*v_muons.at(1).sf*v_muons.at(2).sf*v_electrons.at(0).sf << std::endl;
      int nlep = v_muons.size() + v_electrons.size();

      h_nJets->Fill(nj, weight);
      fillElMuHistCollection4l(Cut1_3Mu1El, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
      h_TotalEvents_1El2Mu->Fill(1, weight);
      if(nb==0)
      {
        fillElMuHistCollection4l(Cut2_3Mu1El, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
        h_TotalEvents_1El2Mu->Fill(2, weight);
        //if((fabs((v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() - 91.1875) < 10 and v_muons.at(0).charge*v_muons.at(1).charge==-1) or (fabs((v_muons.at(1).lep_lv+v_muons.at(2).lep_lv).M() - 91.1875) < 10 and v_muons.at(1).charge*v_muons.at(2).charge==-1) or (fabs((v_muons.at(0).lep_lv+v_muons.at(2).lep_lv).M() - 91.1875) < 10 and v_muons.at(0).charge*v_muons.at(2).charge==-1))
        //if(fabs((v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() - 91.1875) < 10 or fabs((v_muons.at(1).lep_lv+v_muons.at(2).lep_lv).M() - 91.1875) < 10 or fabs((v_muons.at(0).lep_lv+v_muons.at(2).lep_lv).M() - 91.1875) < 10)
        //{
          fillElMuHistCollection4l(Cut3_3Mu1El, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
          h_TotalEvents_1El2Mu->Fill(3, weight);
          if(v_electrons.at(0).lep_lv.Pt() > 15.0 and v_muons.at(0).lep_lv.Pt() > 10.0 and v_muons.at(1).lep_lv.Pt() > 10.0)
          //if(v_electrons.at(0).lep_lv.Pt() > 15.0 and v_muons.at(0).lep_lv.Pt() > 10.0 and v_muons.at(1).lep_lv.Pt() > 10.0 and (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() > 100.0)
          //if(met_pt > 30.0 and (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() > 100.0 and v_electrons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(0).lep_lv.Pt() > 10.0 and v_muons.at(1).lep_lv.Pt() > 10.0)
          //if((fabs((v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M() - 91.1875) > 10 and v_electrons.at(0).charge*v_muons.at(0).charge==-1) or (fabs((v_electrons.at(0).lep_lv+v_muons.at(1).lep_lv).M() - 91.1875) > 10 and v_electrons.at(0).charge*v_muons.at(1).charge==-1) or (fabs((v_electrons.at(0).lep_lv+v_muons.at(2).lep_lv).M() - 91.1875) > 10 and v_electrons.at(0).charge*v_muons.at(2).charge==-1))
          {
            fillElMuHistCollection4l(Cut4_3Mu1El, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
            if(Sample=="Data") std::cout << "run, evt, lumi = " << run << " , " << evt << " , " << lumi << std::endl;
            h_TotalEvents_1El2Mu->Fill(4, weight);
          }//met cut
        //}//z-tag
      }//b-tag
    }//no cuts

    if(v_muons.size()==2 and v_electrons.size()==2 and v_muons.at(0).iso < 0.25 and v_muons.at(1).iso > 0.25 and (elecIsoCut(v_electrons.at(0).eta, v_electrons.at(0).pt, v_electrons.at(0).iso)==true) and (elecIsoCut(v_electrons.at(1).eta, v_electrons.at(1).pt, v_electrons.at(1).iso)==true))
    { 
      if(Sample=="MC") weight *= v_muons.at(0).sf*v_muons.at(1).sf*v_electrons.at(0).sf*v_electrons.at(1).sf;
      else weight = 1.0;
      int nlep = v_muons.size() + v_electrons.size();
      
      h_nJets->Fill(nj, weight);
      fillElMuHistCollection4l(Cut1_2Mu2El_LooseMu, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
      h_TotalEvents_1El2Mu->Fill(1, weight);
      if(nb==0)
      {
        fillElMuHistCollection4l(Cut2_2Mu2El_LooseMu, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight); 
        h_TotalEvents_1El2Mu->Fill(2, weight);
        if(fabs((v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv).M() - 91.1875) < 10 and (v_electrons.at(0).charge*v_electrons.at(1).charge==-1) and fabs((v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() - 91.1875) > 20 and (v_muons.at(0).charge*v_muons.at(1).charge==-1))
        {
          fillElMuHistCollection4l(Cut3_2MuMinus2El_LooseMu, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
          h_TotalEvents_1El2Mu->Fill(3, weight);
          if(v_electrons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(1).lep_lv.Pt() > 10.0 and v_electrons.at(1).lep_lv.Pt() > 10.0)
          //if(met_pt > 30.0 and (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() > 100.0 and v_electrons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(1).lep_lv.Pt() > 10.0 and v_electrons.at(1).lep_lv.Pt() > 10.0)
          {
            double mt1 = sqrt(2*met_pt*v_electrons.at(0).lep_lv.Pt()*(1.0 - cos(v_electrons.at(0).lep_lv.Phi() - met_phi))); 
            double mt2 = sqrt(2*met_pt*v_electrons.at(1).lep_lv.Pt()*(1.0 - cos(v_electrons.at(1).lep_lv.Phi() - met_phi)));
            fillElMuHistCollection4l(Cut4_2MuMinus2El_LooseMu, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, mt1, mt2, (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
            h_TotalEvents_1El2Mu->Fill(4, weight);
          }//met cut
        }//z-tag
        if(fabs((v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv).M() - 91.1875) < 10 and (v_electrons.at(0).charge*v_electrons.at(1).charge==-1) and fabs((v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() - 91.1875) > 20 and (v_muons.at(0).charge*v_muons.at(1).charge==1))
        {
          fillElMuHistCollection4l(Cut3_2MuPlus2El_LooseMu, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
          h_TotalEvents_1El2Mu->Fill(3, weight);
          if(v_electrons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(1).lep_lv.Pt() > 10.0 and v_electrons.at(1).lep_lv.Pt() > 10.0)
          //if(met_pt > 30.0 and (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() > 100.0 and v_electrons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(1).lep_lv.Pt() > 10.0 and v_electrons.at(1).lep_lv.Pt() > 10.0)
          {
            double mt1 = sqrt(2*met_pt*v_electrons.at(0).lep_lv.Pt()*(1.0 - cos(v_electrons.at(0).lep_lv.Phi() - met_phi)));
            double mt2 = sqrt(2*met_pt*v_electrons.at(1).lep_lv.Pt()*(1.0 - cos(v_electrons.at(1).lep_lv.Phi() - met_phi)));
            fillElMuHistCollection4l(Cut4_2MuPlus2El_LooseMu, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, mt1, mt2, (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
            h_TotalEvents_1El2Mu->Fill(4, weight);
          }//met cut
        }//z-tag
      }//b-tag
    }//no cuts

    if(v_muons.size()==2 and v_electrons.size()==2 and v_muons.at(0).iso < 0.25 and v_muons.at(1).iso < 0.25 and (elecIsoCut(v_electrons.at(0).eta, v_electrons.at(0).pt, v_electrons.at(0).iso)==true) and (elecIsoCut(v_electrons.at(1).eta, v_electrons.at(1).pt, v_electrons.at(1).iso)==false))
    {
      if(Sample=="MC") weight *= v_muons.at(0).sf*v_muons.at(1).sf*v_electrons.at(0).sf*v_electrons.at(1).sf;
      else weight = 1.0;
      int nlep = v_muons.size() + v_electrons.size();

      h_nJets->Fill(nj, weight);
      fillElMuHistCollection4l(Cut1_2Mu2El_LooseEl, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
      h_TotalEvents_1El2Mu->Fill(1, weight);
      if(nb==0)
      {
        fillElMuHistCollection4l(Cut2_2Mu2El_LooseEl, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
        h_TotalEvents_1El2Mu->Fill(2, weight);
        if(fabs((v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv).M() - 91.1875) > 20 and (v_electrons.at(0).charge*v_electrons.at(1).charge==-1) and fabs((v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() - 91.1875) < 10 and (v_muons.at(0).charge*v_muons.at(1).charge==-1))
        {
          fillElMuHistCollection4l(Cut3_2Mu2ElMinus_LooseEl, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
          h_TotalEvents_1El2Mu->Fill(3, weight);
          if(met_pt > 60.0 and v_electrons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(1).lep_lv.Pt() > 10.0 and v_electrons.at(1).lep_lv.Pt() > 10.0)
          //if(met_pt > 30.0 and (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() > 100.0 and v_electrons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(1).lep_lv.Pt() > 10.0 and v_electrons.at(1).lep_lv.Pt() > 10.0)
          {
            double mt1 = sqrt(2*met_pt*v_electrons.at(0).lep_lv.Pt()*(1.0 - cos(v_electrons.at(0).lep_lv.Phi() - met_phi)));
            double mt2 = sqrt(2*met_pt*v_electrons.at(1).lep_lv.Pt()*(1.0 - cos(v_electrons.at(1).lep_lv.Phi() - met_phi)));
            fillElMuHistCollection4l(Cut4_2Mu2ElMinus_LooseEl, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, mt1, mt2, (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
            h_TotalEvents_1El2Mu->Fill(4, weight);
          }//met cut
        }//z-tag
        if(fabs((v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv).M() - 91.1875) > 20 and (v_electrons.at(0).charge*v_electrons.at(1).charge==1) and fabs((v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() - 91.1875) < 10 and (v_muons.at(0).charge*v_muons.at(1).charge==-1))
        {
          fillElMuHistCollection4l(Cut3_2Mu2ElPlus_LooseEl, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
          h_TotalEvents_1El2Mu->Fill(3, weight);
          if(met_pt > 60.0 and v_electrons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(1).lep_lv.Pt() > 10.0 and v_electrons.at(1).lep_lv.Pt() > 10.0)
          //if(met_pt > 30.0 and (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() > 100.0 and v_electrons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(1).lep_lv.Pt() > 10.0 and v_electrons.at(1).lep_lv.Pt() > 10.0)
          {
            double mt1 = sqrt(2*met_pt*v_electrons.at(0).lep_lv.Pt()*(1.0 - cos(v_electrons.at(0).lep_lv.Phi() - met_phi)));
            double mt2 = sqrt(2*met_pt*v_electrons.at(1).lep_lv.Pt()*(1.0 - cos(v_electrons.at(1).lep_lv.Phi() - met_phi)));
            fillElMuHistCollection4l(Cut4_2Mu2ElPlus_LooseEl, v_electrons.at(0).pt, 0.0, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, 0.0, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, 0.0, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, 0.0, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, 0.0, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, 0.0, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, 0.0, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, 0.0, 0.0, mt1, mt2, (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
            h_TotalEvents_1El2Mu->Fill(4, weight);
          }//met cut
        }//z-tag
      }//b-tag
    }//no cuts  

    //if(v_muons.size()==1 and v_electrons.size()==3 and v_muons.at(0).iso > 0.25 and (elecIsoCut(v_electrons.at(0).eta, v_electrons.at(0).pt, v_electrons.at(0).iso)==true) and (elecIsoCut(v_electrons.at(1).eta, v_electrons.at(1).pt, v_electrons.at(1).iso)==true) and (elecIsoCut(v_electrons.at(2).eta, v_electrons.at(2).pt, v_electrons.at(2).iso)==true))
    /*if(v_muons.size()==2 and v_electrons.size()==2 and v_muons.at(0).iso < 0.25 and v_muons.at(1).iso < 0.25 and (elecIsoCut(v_electrons.at(0).eta, v_electrons.at(0).pt, v_electrons.at(0).iso)==true) and (elecIsoCut(v_electrons.at(1).eta, v_electrons.at(1).pt, v_electrons.at(1).iso)==false))
    {
      if(Sample=="MC") weight *= v_muons.at(0).sf*v_electrons.at(0).sf*v_electrons.at(1).sf*v_electrons.at(2).sf;
      else weight = 1.0;
      int nlep = v_muons.size() + v_electrons.size();

      h_nJets->Fill(nj, weight);
      fillElMuHistCollection4l(Cut1_2Mu2El_LooseEl, v_electrons.at(0).pt, v_electrons.at(1).pt, v_electrons.at(2).pt, v_muons.at(0).pt, v_electrons.at(0).eta, v_electrons.at(1).eta, v_electrons.at(2).eta, v_muons.at(0).eta, v_electrons.at(0).phi, v_electrons.at(1).phi, v_electrons.at(2).phi, v_muons.at(0).phi, v_electrons.at(0).sip3d, v_electrons.at(1).sip3d, v_electrons.at(2).sip3d, v_muons.at(0).sip3d, v_electrons.at(0).dxy, v_electrons.at(1).dxy, v_electrons.at(2).dxy, v_muons.at(0).dxy, v_electrons.at(0).dz, v_electrons.at(1).dz, v_electrons.at(2).dz, v_muons.at(0).dz, v_electrons.at(0).iso, v_electrons.at(1).iso, v_electrons.at(2).iso, v_muons.at(0).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv+v_electrons.at(2).lep_lv).M(), (v_muons.at(0).lep_lv).M(), weight);
      h_TotalEvents_1El2Mu->Fill(1, weight);
      if(nb==0)
      {
        fillElMuHistCollection4l(Cut2_2Mu2El_LooseEl, v_electrons.at(0).pt, v_electrons.at(1).pt, v_electrons.at(2).pt, v_muons.at(0).pt, v_electrons.at(0).eta, v_electrons.at(1).eta, v_electrons.at(2).eta, v_muons.at(0).eta, v_electrons.at(0).phi, v_electrons.at(1).phi, v_electrons.at(2).phi, v_muons.at(0).phi, v_electrons.at(0).sip3d, v_electrons.at(1).sip3d, v_electrons.at(2).sip3d, v_muons.at(0).sip3d, v_electrons.at(0).dxy, v_electrons.at(1).dxy, v_electrons.at(2).dxy, v_muons.at(0).dxy, v_electrons.at(0).dz, v_electrons.at(1).dz, v_electrons.at(2).dz, v_muons.at(0).dz, v_electrons.at(0).iso, v_electrons.at(1).iso, v_electrons.at(2).iso, v_muons.at(0).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv+v_electrons.at(2).lep_lv).M(), (v_muons.at(0).lep_lv).M(), weight);
        h_TotalEvents_1El2Mu->Fill(2, weight);
        //if(fabs((v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv).M() - 91.1875) < 10 or fabs((v_electrons.at(1).lep_lv+v_electrons.at(2).lep_lv).M() - 91.1875) < 10 or fabs((v_electrons.at(0).lep_lv+v_electrons.at(2).lep_lv).M() - 91.1875) < 10)
        //if(fabs((v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv).M() - 91.1875) < 10 and fabs((v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() - 91.1875) > 20)
        //if(fabs((v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv).M() - 91.1875) > 20 and fabs((v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() - 91.1875) < 10)
        {
          fillElMuHistCollection4l(Cut3_2Mu2El_LooseEl, v_electrons.at(0).pt, v_electrons.at(1).pt, v_electrons.at(2).pt, v_muons.at(0).pt, v_electrons.at(0).eta, v_electrons.at(1).eta, v_electrons.at(2).eta, v_muons.at(0).eta, v_electrons.at(0).phi, v_electrons.at(1).phi, v_electrons.at(2).phi, v_muons.at(0).phi, v_electrons.at(0).sip3d, v_electrons.at(1).sip3d, v_electrons.at(2).sip3d, v_muons.at(0).sip3d, v_electrons.at(0).dxy, v_electrons.at(1).dxy, v_electrons.at(2).dxy, v_muons.at(0).dxy, v_electrons.at(0).dz, v_electrons.at(1).dz, v_electrons.at(2).dz, v_muons.at(0).dz, v_electrons.at(0).iso, v_electrons.at(1).iso, v_electrons.at(2).iso, v_muons.at(0).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv+v_electrons.at(2).lep_lv).M(), (v_muons.at(0).lep_lv).M(), weight);
          h_TotalEvents_1El2Mu->Fill(3, weight);
          if(met_pt > 30.0 and (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv+v_electrons.at(1).lep_lv).M() > 100.0 and v_electrons.at(0).lep_lv.Pt() > 25.0 and v_muons.at(0).lep_lv.Pt() > 25.0 and v_electrons.at(1).lep_lv.Pt() > 10.0 and v_electrons.at(2).lep_lv.Pt() > 10.0)
          {
            double mt1 = sqrt(2*met_pt*v_electrons.at(0).lep_lv.Pt()*(1.0 - cos(v_electrons.at(0).lep_lv.Phi() - met_phi)));
            double mt2 = sqrt(2*met_pt*v_electrons.at(1).lep_lv.Pt()*(1.0 - cos(v_electrons.at(1).lep_lv.Phi() - met_phi)));
            fillElMuHistCollection4l(Cut4_2Mu2El_LooseEl, v_electrons.at(0).pt, v_electrons.at(1).pt, v_electrons.at(2).pt, v_muons.at(0).pt, v_electrons.at(0).eta, v_electrons.at(1).eta, v_electrons.at(2).eta, v_muons.at(0).eta, v_electrons.at(0).phi, v_electrons.at(1).phi, v_electrons.at(2).phi, v_muons.at(0).phi, v_electrons.at(0).sip3d, v_electrons.at(1).sip3d, v_electrons.at(2).sip3d, v_muons.at(0).sip3d, v_electrons.at(0).dxy, v_electrons.at(1).dxy, v_electrons.at(2).dxy, v_muons.at(0).dxy, v_electrons.at(0).dz, v_electrons.at(1).dz, v_electrons.at(2).dz, v_muons.at(0).dz, v_electrons.at(0).iso, v_electrons.at(1).iso, v_electrons.at(2).iso, v_muons.at(0).iso, met_pt, nj, nb, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_muons.at(0).lep_lv).M(), (v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv+v_electrons.at(2).lep_lv).M(), (v_muons.at(0).lep_lv).M(), weight);
            h_TotalEvents_1El2Mu->Fill(4, weight);
          }//met cut
        }//z-tag
      }//b-tag
    }//no cuts 
*/
    if(v_leptons.size()>=6) h_TotalEvents_6l->Fill(1, weight); 
    
  }//event loop
  std::cout << "n_4mu = " << n_4mu << std::endl;
  std::cout << "n_3mu = " << n_3mu << std::endl;
  std::cout << "n_1el = " << n_1el << std::endl;
  std::string histfilename=("output_"+infile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  tFile->cd();
  tFile->mkdir("Cut1_3Mu1El");
  tFile->cd("Cut1_3Mu1El");
  writeHistCollection(Cut1_3Mu1El);
  tFile->cd();
  tFile->mkdir("Cut2_3Mu1El");
  tFile->cd("Cut2_3Mu1El");
  writeHistCollection(Cut2_3Mu1El);
  tFile->cd();
  tFile->mkdir("Cut3_3Mu1El");
  tFile->cd("Cut3_3Mu1El");
  writeHistCollection(Cut3_3Mu1El);
  tFile->cd();
  tFile->mkdir("Cut4_3Mu1El");
  tFile->cd("Cut4_3Mu1El");
  writeHistCollection(Cut4_3Mu1El);
  tFile->cd();
  tFile->mkdir("Cut1_2Mu2El_LooseMu");
  tFile->cd("Cut1_2Mu2El_LooseMu");
  writeHistCollection(Cut1_2Mu2El_LooseMu);
  tFile->cd();
  tFile->mkdir("Cut2_2Mu2El_LooseMu");
  tFile->cd("Cut2_2Mu2El_LooseMu");
  writeHistCollection(Cut2_2Mu2El_LooseMu);
  tFile->cd();
  tFile->mkdir("Cut3_2MuPlus2El_LooseMu");
  tFile->cd("Cut3_2MuPlus2El_LooseMu");
  writeHistCollection(Cut3_2MuPlus2El_LooseMu);
  tFile->cd();
  tFile->mkdir("Cut4_2MuPlus2El_LooseMu");
  tFile->cd("Cut4_2MuPlus2El_LooseMu");
  writeHistCollection(Cut4_2MuPlus2El_LooseMu);
  tFile->cd();
  tFile->mkdir("Cut3_2MuMinus2El_LooseMu");
  tFile->cd("Cut3_2MuMinus2El_LooseMu");
  writeHistCollection(Cut3_2MuMinus2El_LooseMu);
  tFile->cd();
  tFile->mkdir("Cut4_2MuMinus2El_LooseMu");
  tFile->cd("Cut4_2MuMinus2El_LooseMu");
  writeHistCollection(Cut4_2MuMinus2El_LooseMu);
  tFile->cd();
  tFile->mkdir("Cut1_2Mu2El_LooseEl");
  tFile->cd("Cut1_2Mu2El_LooseEl");
  writeHistCollection(Cut1_2Mu2El_LooseEl);
  tFile->cd();
  tFile->mkdir("Cut2_2Mu2El_LooseEl");
  tFile->cd("Cut2_2Mu2El_LooseEl");
  writeHistCollection(Cut2_2Mu2El_LooseEl);
  tFile->cd();
  tFile->mkdir("Cut3_2Mu2ElPlus_LooseEl");
  tFile->cd("Cut3_2Mu2ElPlus_LooseEl");
  writeHistCollection(Cut3_2Mu2ElPlus_LooseEl);
  tFile->cd();
  tFile->mkdir("Cut4_2Mu2ElPlus_LooseEl");
  tFile->cd("Cut4_2Mu2ElPlus_LooseEl");
  writeHistCollection(Cut4_2Mu2ElPlus_LooseEl);
  tFile->cd();
  tFile->mkdir("Cut3_2Mu2ElMinus_LooseEl");
  tFile->cd("Cut3_2Mu2ElMinus_LooseEl");
  writeHistCollection(Cut3_2Mu2ElMinus_LooseEl);
  tFile->cd();
  tFile->mkdir("Cut4_2Mu2ElMinus_LooseEl");
  tFile->cd("Cut4_2Mu2ElMinus_LooseEl");
  writeHistCollection(Cut4_2Mu2ElMinus_LooseEl);
  tFile->cd();
  h_zmass_1->Write();
  h_zmass_2->Write();
  h_mT->Write();
  h_nLeptons->Write();
  h_ee_mm_ElMu->Write(); 
  h_TotalEvents_4Mu->Write();
  h_TotalEvents_1El2Mu->Write(); 
  h_TotalEvents_2El1Mu->Write();
  h_TotalEvents_4El->Write(); 
  h_TotalEvents_6l->Write();
  h_ee_ee_4El->Write();
  h_mm_mm_4Mu->Write();
  h_nJets->Write();
  tFile->Close();
  sfElRecoFile->Close();
  sfElVetoFile->Close();
  sfMuRecoFile->Close();
  sfMuRecoLowPtFile->Close();
  sfMuIsoFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;
}
