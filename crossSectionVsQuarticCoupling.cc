#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include <TCanvas.h>
#include <TLatex.h>
#include "TGraphErrors.h"
#include "TLegend.h"
#include <TPad.h>
#include <sstream>
#include "TVectorD.h"
#include "TGraph.h"
#include "TFile.h"
#include "THStack.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TColor.h"

using std::string;
using std::cout;
using std::endl;
using std::istringstream;

double BR = (0.33*0.33*0.67*3 + 0.33*0.33*0.33);

void crossSectionVsQuarticCoupling()
{

  TCanvas c1("c1","kLong energy versus hits", 10, 10, 600, 400);
  gStyle->SetOptFit(1);
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  c1.SetGrid();

  double quartic[7] = {1.00, 0.995, 0.99, 0.985, 0.98, 0.975, 0.97}; 
  double crosssection[7] = {0.1125, 0.1127, 0.116, 0.1179, 0.1203, 0.1262, 0.1341};


  TGraph *gr1obs = new TGraph(7, quartic, crosssection);
  gr1obs->SetLineStyle(5);
  gr1obs->SetLineColor(kRed);
  gr1obs->SetLineWidth(2);
  gr1obs->SetMarkerColor(kRed);
  gr1obs->SetTitle("");
  gr1obs->GetXaxis()->SetTitle("Strength of quartic coupling");
  gr1obs->GetYaxis()->CenterTitle();
  gr1obs->GetYaxis()->SetTitleOffset(1.3);
  gr1obs->GetYaxis()->SetTitle("Cross section [pb]");
  gr1obs->SetMaximum(0.145);
  gr1obs->SetMinimum(0.105);
  gr1obs->Draw("APL* e");

  TLegend *leg = new TLegend(0.60,0.75,0.85,0.85,NULL,"brNDC");
  leg->AddEntry(gr1obs, "Cross section", "l");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  TLatex TL;
  TL.SetTextAlign(11);
  TL.SetTextSize(0.05);
  TL.SetTextFont(22);
  TL.DrawLatexNDC(0.35,0.94,"Strength of quartic coupling versus cross section");
  

  c1.SaveAs("crosssection_quartic.pdf");
  c1.SaveAs("crosssection_quartic.png");  
   

} 
