#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TProfile.h"
#include "TF1.h"
#include "TKey.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooBinning.h"

//default drawing options
void Plot(){
 TFile* file0= new TFile ("/afs/cern.ch/user/m/mshi/CMSSW_7_4_12_patch4/src/GGHAA2Mu2TauAnalysis/MuMuTauTauGenAnalyzer/python/histodemo_heavyHiggs_125_light19.root"); 
 TFile* file1= new TFile ("/afs/cern.ch/user/m/mshi/CMSSW_7_4_12_patch4/src/GGHAA2Mu2TauAnalysis/MuMuTauTauGenAnalyzer/python/Drell_Yan.root");
  THStack *hs= new THStack("hs","DeltaR between generation level muons decayed from same light higgs(signal, with red) or z(DrellYan, blue color) ");
  TH1F *h1st=(TH1F*)file0->Get("demo/DeltaR_of_2_mu_decayed_from_same_higgs");
  Double_t norm=8;
  h1st->Scale(norm/h1st->Integral(), "width");
   h1st->SetFillColor(kRed);
   h1st->SetMarkerStyle(21);
   h1st->SetMarkerColor(kRed);
  std::cout<<"test"<<std::endl;
   hs->Add(h1st);
  std::cout<<"test1"<<std::endl;
  TH1F *h2st=(TH1F*)file1->Get("DrellYan/DeltaR of 2 mu decayed from z");
  Double_t norm1=4000000;
  h2st->Scale(norm1/h2st->Integral(),"width");
  h2st->SetFillColor(kBlue);
  h2st->SetMarkerStyle(21);
  h2st->SetMarkerColor(kBlue);
  hs->Add(h2st);
  
  TCanvas* cst=new TCanvas("cst","stacked hists",10,10,700,700);
   cst->SetLogy();
   cst->SetFillColor(41);
   cst->cd();
   hs->SetMaximum(10000000);
   hs->SetMinimum(1);
   hs->Draw();
   leg = new TLegend(0.1,0.7,0.48,0.9);
   leg->SetHeader("The Legend Title");
   leg->AddEntry(h1st,"ggh-> aa, with a->mumu, another a->tautau","f");
   leg->AddEntry(h2st,"DrellYan with z->mumu","f");
   leg->Draw();

}
