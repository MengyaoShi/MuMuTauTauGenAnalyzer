// -*- C++ -*-
//
// Package:    Amumu/DrellYanAnalyzer
// Class:      DrellYanAnalyzer
// 
/**\class DrellYanAnalyzer DrellYanAnalyzer.cc Amumu/DrellYanAnalyzer/plugins/DrellYanAnalyzer.cc

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
#include "TH1F.h"
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
#include "Tools/Common/interface/Common.h"
//
// class declaration
//
using namespace std;
using namespace edm;
using namespace reco;

class DrellYanAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DrellYanAnalyzer(const edm::ParameterSet&);
      ~DrellYanAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      edm::EDGetTokenT<reco::GenParticleCollection> genParticleTag_;
      std::map<std::string, TH1D*> histos1D_;
      TH1F* DR_DrellYan_;
      TFile* out_;
      std::string outFileName_; 
      TH1F * HighestPTMuonPT_;
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
DrellYanAnalyzer::DrellYanAnalyzer(const edm::ParameterSet& iConfig):
 //jetOutputFileName_(iConfig.getParameter<std::string>("jetOutputFileName")),
genParticleTag_(consumes<reco::GenParticleCollection >(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
 histos1D_(),
outFileName_(iConfig.getParameter<std::string>("outFileName"))
{
}
DrellYanAnalyzer::~DrellYanAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   //    // (e.g. close files, deallocate resources etc.)
   //
}

void DrellYanAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //now do what ever initialization is needed
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByToken(genParticleTag_, pGenParticles);
  double HighestPTMuonPT=(pGenParticles->begin())->pt();
  double HighestPTMuonEta=(pGenParticles->begin())->eta();
  double Temp=0;
  for(reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); iGenParticle != pGenParticles->end(); ++iGenParticle){
    if(fabs((*iGenParticle).pdgId())==13)
    {
      Temp=(*iGenParticle).pt();
      if(HighestPTMuonPT<Temp)
        {
          HighestPTMuonPT=Temp;
          HighestPTMuonEta=((*iGenParticle).eta());
        }
      std::cout<<"You got a mu! it's mother pdg id  =="<<iGenParticle->motherRef()->pdgId()<<std::endl;
    }

    //this is from top to bottom
    if(fabs((*iGenParticle).pdgId()) == 23&& iGenParticle->numberOfDaughters()==2 )
    {
      if(fabs(iGenParticle->daughterRef(0)->pdgId())==13) 
      {
         reco::GenParticleRef child0 = iGenParticle->daughterRef(0);
         reco::GenParticleRef child1 = iGenParticle->daughterRef(1);
        double DR_mu=VariousFunctions::getDiThingDR_1(child0, child1);
        double InvMass=VariousFunctions::getInvMass(child0,child1);
        histos1D_[ "DeltaR of 2 mu decayed from z" ]->Fill(DR_mu,1);
        histos1D_[ "Invariant mass of 2 mu decayed from z" ]->Fill(InvMass);
        histos1D_[ "Invariant mass of 2mu from z or gamma"]->Fill(InvMass);
        DR_DrellYan_->Fill(DR_mu);
      }
    }
    if(fabs((*iGenParticle).pdgId())==22 &&iGenParticle->numberOfDaughters()==2)
    {
      if(fabs(iGenParticle->daughterRef(0)->pdgId())==13)
      {
         reco::GenParticleRef child0 = iGenParticle->daughterRef(0);
         reco::GenParticleRef child1 = iGenParticle->daughterRef(1);
        double InvMass=VariousFunctions::getInvMass(child0,child1);
        histos1D_[ "Invariant mass of 2mu from z or gamma" ]->Fill(InvMass);
      }
    }
  }
  HighestPTMuonPT_->Fill(HighestPTMuonPT);
  histos1D_["HighestPTMuonPT"]->Fill(HighestPTMuonPT);
  histos1D_["HighestPTMuonEta"]->Fill(HighestPTMuonEta); 
}


//
// member functions
//

// ------------ method called for each event  ------------


// ------------ method called once each job just before starting event loop  ------------
void 
DrellYanAnalyzer::beginJob()
{
  out_=new TFile(outFileName_.c_str(), "RECREATE");
  HighestPTMuonPT_=new TH1F("HighestPTMuonPT","HighestPTMuon PT of DrellYan", 100,0,200);
  edm::Service< TFileService > fileService;
 
  histos1D_[ "DeltaR of 2 mu decayed from z"]=fileService->make<TH1D>("DeltaR of 2 mu decayed from z","DeltaR_of_2_mu_decayed_from_z", 100, 0, 5.5); 
  histos1D_[ "Invariant mass of 2 mu decayed from z"]=fileService->make<TH1D>("Invariant mass of 2 mu decayed from z","Invariant mass of 2 mu decayed from z", 100, 0, 300); 
  histos1D_[ "Invariant mass of 2mu from z or gamma"]=fileService->make<TH1D>("Invariant mass of 2mu decayed from z or gamma","Invariant mass of 2 mu decayed from z or gamma", 100, 0, 300);
  histos1D_[ "HighestPTMuonPT"]=fileService->make<TH1D>("HighestPTMuonPT", "HighestPTMuonPT", 100,0,300);
  histos1D_["HighestPTMuonEta"]=fileService->make<TH1D>("HighestPTMuonEta","HighestPTMuonEta", 100, -5.0, 5.0);
  DR_DrellYan_=new TH1F("DR_DrellYan", "",100,0,6.0);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DrellYanAnalyzer::endJob() 
{
 out_->cd();
 // jetOutput_.close();
 TCanvas HighestPTMuonPTcanvas_("canvas", "", 600,600);
 HighestPTMuonPTcanvas_.SetLogy();
 HighestPTMuonPT_->Draw();
// Common::draw1DHistograms(HighestPTMuonPTcanvas_,HighestPTMuonPT_ );
 HighestPTMuonPTcanvas_.Write();
 out_->Write();
 out_->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
DrellYanAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DrellYanAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DrellYanAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DrellYanAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DrellYanAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DrellYanAnalyzer);
