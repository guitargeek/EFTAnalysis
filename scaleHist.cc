#include <TF1.h>
#include <TH1D.h>
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
#include <TGraph.h>
#include <TRandom.h>
#include <TMath.h>
#include <TKey.h>
#include <TDirectory.h>

using std::string;
/*
void scaleAllHistos()
{
   TFile *infile = new TFile("zzz.root", "UPDATE");
   //TFile *outfile = new TFile("zzz_output.root", "RECREATE");
   int i=0;

   TIter nextkey(gDirectory->GetListOfKeys());
   while (TKey *key = (TKey*)nextkey()) 
   {
      //std::cout << "i: " << i << std::endl;
      ++i;
      TObject* obj = key->ReadObj();
      std::string name = obj->GetName();
      TH1D *hist = (TH1D*)infile->Get(name.c_str());
      TH1D *newHist=(TH1D*)hist->Clone(name.c_str());
      newHist->Scale(35.9/59.7);
      newHist->Write();
   }
   //outfile->Write();
}*/

void scaleAllHistos()
{
  TFile *infile = new TFile("zzz.root");

  int i=0;
  std::vector<TH1D*> v_hist;

  TIter nextkey(gDirectory->GetListOfKeys());
  while (TKey *key = (TKey*)nextkey()) 
  {
    ++i;
    TObject* obj = key->ReadObj();
    std::string name = obj->GetName();
    TH1D *hist = (TH1D*)infile->Get(name.c_str());
    hist->Scale(35.9/59.7);
    v_hist.push_back((TH1D*)hist->Clone(name.c_str()));
  }
  TFile *outfile = new TFile("zzz_2016.root", "RECREATE");
  for(unsigned int i=0; i<v_hist.size(); i++)
  {
    v_hist.at(i)->Write();
  }
  outfile->Write();
} 
