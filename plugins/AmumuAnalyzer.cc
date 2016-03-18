// -*- C++ -*-
//
// Package:    Amumu/AmumuAnalyzer
// Class:      AmumuAnalyzer
// 
/**\class AmumuAnalyzer AmumuAnalyzer.cc Amumu/AmumuAnalyzer/plugins/AmumuAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mengyao Shi
//         Created:  Mon, 19 Oct 2015 10:59:36 GMT
//
//


// system include files
#include <memory>
#include <map>
#include <string>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include <math.h>
#include <sstream>
#include <typeinfo>
#include "TCanvas.h"
#include "TFile.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
 #include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "GGHAA2Mu2TauAnalysis/VariousFunctions/interface/VariousFunctions.h"
//
// class declaration
//
using namespace std;
using namespace edm;
using namespace reco;

class AmumuAnalyzer : public edm::EDAnalyzer {
   public:
      explicit AmumuAnalyzer(const edm::ParameterSet&);
      ~AmumuAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void reset(const bool);
      double mass;
      double pt_of_muon_from_tau;
      double pt_of_higher_muon;
      double pt_of_lower_muon;
      TH2F* DeltaR_a_and_tau_VS_pt_of_a_;
      TH2F* DeltaR_two_tau_VS_pt_of_a_;
      TH2F* DeltaR_two_mu_VS_pt_of_a_;
      TH2F* pt_of_muon_from_tau_VS_pt_of_higher_muon_;
      TH2F* pt_of_muon_from_tau_VS_pt_of_lower_muon_;
      double eta_of_higher_pt_A_mu;
      double eta_of_lower_pt_A_mu;
      TH1F* DR_Signal_;
      TFile* out_;
      std::string outFileName_;
      edm::EDGetTokenT<reco::GenParticleCollection> genParticleTag_;
      std::map<std::string, TH1D*> histos1D_;
      std::map<std::string, TH2D*> histos2D_;

      //std::string jetOutputFileName_;
      //ofstream jetOutput_;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
//
AmumuAnalyzer::AmumuAnalyzer(const edm::ParameterSet& iConfig):
 outFileName_(iConfig.getParameter<std::string>("outFileName")),
 //jetOutputFileName_(iConfig.getParameter<std::string>("jetOutputFileName")),
 genParticleTag_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
 histos1D_(),
 histos2D_()
{
  reset(false);
}
AmumuAnalyzer::~AmumuAnalyzer()
{
   reset(true);
   // do anything here that needs to be done at desctruction time
   //    // (e.g. close files, deallocate resources etc.)
   //
}

void AmumuAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //now do what ever initialization is needed
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByToken(genParticleTag_, pGenParticles);
  bool tau_mu=false;
  for(reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); iGenParticle != pGenParticles->end(); ++iGenParticle)
  {
    //this is from bottom to top
    if((fabs(iGenParticle->pdgId())==13) && (fabs(iGenParticle->mother()->pdgId())!=13))
    {
      if(fabs(iGenParticle->mother()->pdgId())==15)
      { 
        if(fabs(iGenParticle->mother()->mother()->pdgId())==36)
        {
          pt_of_muon_from_tau=iGenParticle->pt();
          if(abs(iGenParticle->eta())<2.1)
          histos1D_["pt_of_muon_from_tau"]->Fill(pt_of_muon_from_tau);
          tau_mu=true;
        }
      }
    }

    //this is from top to bottom
    if((*iGenParticle).pdgId() == 35 && ((*iGenParticle).numberOfDaughters()==2) )
    {
      reco::GenParticleRef child0 = iGenParticle->daughterRef(0);
      reco::GenParticleRef child1 = iGenParticle->daughterRef(1);
      reco::GenParticleRef childchild00;
      reco::GenParticleRef childchild01;
      reco::GenParticleRef childchild10;
      reco::GenParticleRef childchild11;
 
      if(child0->pdgId() == 36 && child1->pdgId()==36 && (child0->numberOfDaughters()==2)&&(child1->numberOfDaughters()==2)&&(fabs(child0->daughterRef(0)->pdgId())==13)&&(fabs(child1->daughterRef(0)->pdgId())==15) )
      { 
        childchild00=child0->daughterRef(0);//always make sure childchild0* is muon
        childchild01=child0->daughterRef(1);
        childchild10=child1->daughterRef(0);
        childchild11=child1->daughterRef(1);
      }
      else
      {
        childchild00=child1->daughterRef(0);
        childchild01=child1->daughterRef(1);
        childchild10=child0->daughterRef(0);
        childchild11=child0->daughterRef(1); 
      }
          
      mass=0.0;
      mass=child0->mass();     
      double DR_tau=VariousFunctions::getDiThingDR_1(childchild10, childchild11);
      double DR_mu=VariousFunctions::getDiThingDR_1(childchild00, childchild01);
      double DR_a_and_tau=VariousFunctions::getDiThingDR_1(child1, childchild11);
      double p_of_a=child1->pt();
      reco::GenParticleRef higherPtObj=VariousFunctions::getHigherPtObj(childchild00,childchild01);
      eta_of_higher_pt_A_mu=higherPtObj->eta();
      double DR_tau_mu0=VariousFunctions::getMuTauDR(higherPtObj, childchild10, true);
      double DR_tau_mu1=VariousFunctions::getMuTauDR(higherPtObj, childchild11, true);
      double DR_tau_mu=0;

      reco::GenParticleRef lowerPtObj=VariousFunctions::getLowerPtObj(childchild00, childchild01);
      eta_of_lower_pt_A_mu=lowerPtObj->eta();
      pt_of_higher_muon=VariousFunctions::getHigherPt( childchild00, childchild01);
      pt_of_lower_muon=VariousFunctions::getLowerPt( childchild00, childchild01);  
      if(DR_tau > 2.0)
      {
        std::cout<<"eta" << childchild10->eta() << "phi"<< childchild10->phi()<< "pt of a"<< p_of_a <<std::endl;
        std::cout<<"eta" << childchild11->eta() << "phi"<< childchild11->phi() <<std::endl;
      }
      if(DR_tau_mu0<DR_tau_mu1)
        DR_tau_mu=DR_tau_mu0;
      else
        DR_tau_mu=DR_tau_mu1;
  
      histos1D_[ "test" ]->Fill(child0->mass());
      histos1D_[ "DeltaR_of_2_tau_decayed_from_same_higgs" ]->Fill(DR_tau);
      histos1D_[ "DeltaR_of_2_mu_decayed_from_same_higgs" ]->Fill(DR_mu,1);
      histos1D_[ "DeltaR of highest pt muon and nearest tau"]->Fill(DR_tau_mu);
      histos2D_[ "DeltaR_a_and_tau_VS_p_of_a" ]->Fill(DR_a_and_tau, p_of_a);
      if(abs(higherPtObj->eta())<2.1)
      histos1D_[ "pt_of_higher_muon" ]->Fill(pt_of_higher_muon);

      histos1D_[ "eta_of_higher_pt_A_mu"]->Fill(eta_of_higher_pt_A_mu);
      if(abs(lowerPtObj->eta())<2.1)
      histos1D_[ "pt_of_lower_muon" ]->Fill(pt_of_lower_muon); 

      histos1D_[ "eta_of_lower_pt_A_mu"]->Fill(eta_of_lower_pt_A_mu);
      DeltaR_a_and_tau_VS_pt_of_a_->Fill(DR_a_and_tau, p_of_a);
      DeltaR_two_tau_VS_pt_of_a_->Fill(DR_tau, p_of_a);
      DeltaR_two_mu_VS_pt_of_a_->Fill(DR_mu, p_of_a);
      DR_Signal_->Fill(DR_mu,1);
      std::cout <<"second test of bool"  << tau_mu << std::endl;
      if(tau_mu)
      {
        pt_of_muon_from_tau_VS_pt_of_higher_muon_->Fill(pt_of_muon_from_tau, pt_of_higher_muon); 
        pt_of_muon_from_tau_VS_pt_of_lower_muon_->Fill(pt_of_muon_from_tau, pt_of_lower_muon);
      }
    }
  }
}


void
AmumuAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (DeltaR_a_and_tau_VS_pt_of_a_ != NULL)) delete DeltaR_a_and_tau_VS_pt_of_a_;
  DeltaR_a_and_tau_VS_pt_of_a_ = NULL;
  if ((doDelete) && (DeltaR_two_tau_VS_pt_of_a_ != NULL)) delete DeltaR_two_tau_VS_pt_of_a_;
  DeltaR_two_tau_VS_pt_of_a_ = NULL;
  if ((doDelete) && (DeltaR_two_mu_VS_pt_of_a_ != NULL)) delete DeltaR_two_mu_VS_pt_of_a_;
  DeltaR_two_mu_VS_pt_of_a_ = NULL;
  if ((doDelete) && (pt_of_muon_from_tau_VS_pt_of_higher_muon_ != NULL)) delete pt_of_muon_from_tau_VS_pt_of_higher_muon_;
  pt_of_muon_from_tau_VS_pt_of_higher_muon_ = NULL;
  if ((doDelete) && (pt_of_muon_from_tau_VS_pt_of_lower_muon_ != NULL)) delete pt_of_muon_from_tau_VS_pt_of_lower_muon_;
  pt_of_muon_from_tau_VS_pt_of_lower_muon_ = NULL;


 
}
//
// member functions
//

// ------------ method called for each event  ------------


// ------------ method called once each job just before starting event loop  ------------
void 
AmumuAnalyzer::beginJob()
{
  out_=new TFile(outFileName_.c_str(), "RECREATE");
  //jetOutput_.open(jetOutputFileName_.c_str());
  edm::Service< TFileService > fileService;
  histos1D_[ "test" ]=fileService->make< TH1D >("test", "invariant mass of a from DiMuon", 30, 0, 10);
  histos1D_[ "DeltaR_of_2_tau_decayed_from_same_higgs" ]=fileService->make<TH1D>("DeltaR_of_2_tau_decayed_from_same_higgs", "DeltaR_of_2_tau_decayed_from_same_higgs H300a9", 100, 0.0, 0.5);
  histos1D_[ "DeltaR_of_2_mu_decayed_from_same_higgs"]=fileService->make<TH1D>("DeltaR_of_2_mu_decayed_from_same_higgs","DeltaR_of_2_mu_decayed_from_same_higgs DrellY", 100, 0, 6.0);
  histos1D_[ "DeltaR of highest pt muon and nearest tau"]=fileService->make<TH1D>("DeltaR of highest pt muon and nearest tau", "DeltaR of highest pt muon and nearest tau", 100, 0, 5.0);
  histos2D_[ "DeltaR_a_and_tau_VS_p_of_a" ]=fileService->make<TH2D>("DeltaR_a_and_tau_VS_p_of_a", "DeltaR_a_and_tau_VS_p_of_a DrellY", 100, 0, 0.5, 100, 0, 500);
  histos1D_[ "pt_of_higher_muon"]=fileService->make<TH1D>("pt_of_higher_muon","pt_of_higher_muon DrellY",100,0,500); 
  histos1D_[ "pt_of_muon_from_tau" ]=fileService->make<TH1D>("pt_of_muon_from_tau", "pt_of_muon_from_tau DrellY",100,0,500);
  histos1D_[ "pt_of_lower_muon"]=fileService->make<TH1D>("pt_of_lower_muon","pt_of_lower_muon DrellY",100,0,500);
  histos1D_[ "eta_of_higher_pt_A_mu" ] =fileService->make<TH1D>("eta_of_higher_pt_A_mu", "eta_of_higher_pt_A_mu DrellY", 100, -5.0, 5.0);
  histos1D_[ "eta_of_lower_pt_A_mu" ] =fileService->make<TH1D>("eta_of_lower_pt_A_mu", "eta_of_lower_pt_A_mu DrellY", 100, -5.0, 5.0);

  DeltaR_a_and_tau_VS_pt_of_a_ = new TH2F("DeltaR_a_and_tau_VS_pt_of_a", "heavy300GeV_light9GeV", 100, 0, 0.5, 100, 0, 500);
  DeltaR_two_tau_VS_pt_of_a_ = new TH2F("DeltaR_two_tau_VS_pt_of_a", "heavy300GeV_light9GeV", 100, 0, 0.5, 100, 0, 500);
  DeltaR_two_mu_VS_pt_of_a_ = new TH2F("DeltaR_two_mu_VS_pt_of_a", "heavy300GeV_light9GeV", 100, 0, 0.5, 100, 0, 500);
  pt_of_muon_from_tau_VS_pt_of_higher_muon_ = new TH2F(" pt_of_muon_from_tau_VS_pt_of_higher_muon", "DrellY", 100, 0, 300, 100, 0, 500);
  pt_of_muon_from_tau_VS_pt_of_lower_muon_ = new TH2F(" pt_of_muon_from_tau_VS_pt_of_lower_muon", "DrellY", 100, 0, 300, 100, 0, 500);
  DR_Signal_=new TH1F("DR_Signal", "DR_Signal",100,0.0,4.0);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AmumuAnalyzer::endJob() 
{
 // jetOutput_.close();
  TCanvas DeltaR_a_and_tau_VS_pt_of_a_Canvas("CanvasName","",600,600);
  TCanvas DeltaR_two_tau_VS_pt_of_a_Canvas("Canvas2tau","",600,600);
  TCanvas DeltaR_two_mu_VS_pt_of_a_Canvas("Canvas2mu","",600,600);
  TCanvas pt_of_muon_from_tau_VS_pt_of_higher_muon_Canvas("Canvas_muon_pts","",600,600);
  TCanvas pt_of_muon_from_tau_VS_pt_of_lower_muon_Canvas("Canvas_muonl_pts","",600,600);
  TCanvas Hstack("cst","",10,10,600,600);
  VariousFunctions::formatAndDrawCanvasAndHist2D(DeltaR_a_and_tau_VS_pt_of_a_Canvas, DeltaR_a_and_tau_VS_pt_of_a_, 0, 0, 0, kBlack, 7, 20, "Delta_R_of_a_and_tau", .04, .04, 1.1, "pt of a", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(DeltaR_two_tau_VS_pt_of_a_Canvas, DeltaR_two_tau_VS_pt_of_a_, 0, 0, 0, kBlack, 7, 20, "Delta_R_of_two_tau", .04, .04, 1.1, "pt of a", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(DeltaR_two_mu_VS_pt_of_a_Canvas, DeltaR_two_mu_VS_pt_of_a_, 0, 0, 0, kBlack, 7, 20, "Delta_R_of_two_mu", .04, .04, 1.1, "pt of a", .04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(pt_of_muon_from_tau_VS_pt_of_higher_muon_Canvas, pt_of_muon_from_tau_VS_pt_of_higher_muon_, 0, 0, 0, kBlack, 7, 20, "pt_of_muon_from_tau", .04, .04, 1.1, "pt_of_higher_muon_from_a",.04, .04, 1.6, "", .04, .04, 1.0);
  VariousFunctions::formatAndDrawCanvasAndHist2D(pt_of_muon_from_tau_VS_pt_of_lower_muon_Canvas, pt_of_muon_from_tau_VS_pt_of_lower_muon_, 0, 0, 0, kBlack, 7, 20, "pt_of_muon_from_tau", .04, .04, 1.1, "pt_of_lower_muon_from_a",.04, .04, 1.6, "", .04, .04, 1.0);
  out_->cd();
  DeltaR_a_and_tau_VS_pt_of_a_Canvas.Write();
  DeltaR_two_tau_VS_pt_of_a_Canvas.Write();
  DeltaR_two_mu_VS_pt_of_a_Canvas.Write();
  pt_of_muon_from_tau_VS_pt_of_higher_muon_Canvas.Write();
  pt_of_muon_from_tau_VS_pt_of_lower_muon_Canvas.Write(); 
  Hstack.Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
AmumuAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
AmumuAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
AmumuAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
AmumuAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AmumuAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AmumuAnalyzer);
