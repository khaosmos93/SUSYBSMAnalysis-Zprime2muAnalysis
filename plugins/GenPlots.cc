#include "TH1F.h"
#include "TMath.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"

class GenPlots : public edm::EDAnalyzer {
 public:
  explicit GenPlots(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  const bool isMC;
  HardInteraction* hardInteraction;
  double eventWeight;
  bool useMadgraphWeight;
  double madgraphWeight;
  const bool both_in_acc;

  TH1F* Gen_Weight;
  TH1F* res_mass;
  TH1F* dil_mass;
  TH1F* res_pt;
  TH1F* dil_pt;
  TH1F* res_rap;
  TH1F* dil_rap;
  TH1F* res_eta;
  TH1F* dil_eta;
  TH1F* res_phi;
  TH1F* dil_phi;
};

GenPlots::GenPlots(const edm::ParameterSet& cfg)
  : isMC(cfg.getParameter<bool>("isMC")),
    hardInteraction(isMC ? new HardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")) : 0),
    //hardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")),
    useMadgraphWeight(cfg.getParameter<bool>("useMadgraphWeight")),
    madgraphWeight(1.0),
    both_in_acc(cfg.getParameter<bool>("both_in_acc"))
{
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));
  if (isMC) consumes<std::vector<reco::GenParticle>>(hardInteraction->src);

  edm::Service<TFileService> fs;
  Gen_Weight = fs->make<TH1F>("Gen_Weight", "", 4, -2, 2);
  res_mass = fs->make<TH1F>("res_mass", "", 20000, 0, 20000);
  dil_mass = fs->make<TH1F>("dil_mass", "", 20000, 0, 20000);
  res_pt = fs->make<TH1F>("res_pt", "", 2000, 0, 2000);
  dil_pt = fs->make<TH1F>("dil_pt", "", 2000, 0, 2000);
  res_rap = fs->make<TH1F>("res_rap", "", 100, -5, 5);
  dil_rap = fs->make<TH1F>("dil_rap", "", 100, -5, 5);
  res_eta = fs->make<TH1F>("res_eta", "", 100, -5, 5);
  dil_eta = fs->make<TH1F>("dil_eta", "", 100, -5, 5);
  res_phi = fs->make<TH1F>("res_phi", "", 100, -TMath::Pi(), TMath::Pi());
  dil_phi = fs->make<TH1F>("dil_phi", "", 100, -TMath::Pi(), TMath::Pi());
}

void GenPlots::analyze(const edm::Event& event, const edm::EventSetup& setup) {

  if(isMC) {

    eventWeight = 1.0;
    madgraphWeight = 1.0;

    if (useMadgraphWeight) {
      edm::Handle<GenEventInfoProduct> gen_ev_info;
      event.getByLabel(edm::InputTag("generator"), gen_ev_info);
      if (gen_ev_info.isValid() ){
        eventWeight = gen_ev_info->weight();
        madgraphWeight = ( eventWeight > 0 ) ? 1.0 : -1.0;
      }
    }
    else {
      eventWeight = 1.0;
      madgraphWeight = 1.0;
    }

    Gen_Weight->Fill(madgraphWeight);


    hardInteraction->Fill(event);

    if (!hardInteraction->IsValid() ) {
      //edm::LogWarning("GenPlots") << "!hardInteraction->IsValid()";
      return;
    }

    if (both_in_acc && !(fabs(hardInteraction->lepMinus->eta()) < 2.4 && fabs(hardInteraction->lepPlus->eta()) < 2.4))
      return;

    res_mass->Fill(hardInteraction->resonance->mass(), madgraphWeight);
    res_pt->Fill(hardInteraction->resonance->pt(), madgraphWeight);
    res_eta->Fill(hardInteraction->resonance->eta(), madgraphWeight);
    res_phi->Fill(hardInteraction->resonance->phi(), madgraphWeight);
    res_rap->Fill(hardInteraction->resonance->rapidity(), madgraphWeight);

    reco::Particle::LorentzVector dil = hardInteraction->dilepton();
    dil_mass->Fill(dil.mass(), madgraphWeight);
    dil_pt->Fill(dil.pt(), madgraphWeight);
    dil_eta->Fill(dil.eta(), madgraphWeight);
    dil_phi->Fill(dil.phi(), madgraphWeight);
    dil_rap->Fill(dil.Rapidity(), madgraphWeight);
  }
}

DEFINE_FWK_MODULE(GenPlots);
