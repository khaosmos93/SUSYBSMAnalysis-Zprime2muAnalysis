#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"

class MuonsForNminus1 : public edm::EDProducer {
public:
  explicit MuonsForNminus1(const edm::ParameterSet&);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::InputTag muon_src;
  edm::InputTag dimuon_src;
};

MuonsForNminus1::MuonsForNminus1(const edm::ParameterSet& cfg)
  : muon_src(cfg.getParameter<edm::InputTag>("muon_src")),
    dimuon_src(cfg.getParameter<edm::InputTag>("dimuon_src"))
{
  consumes<pat::MuonCollection>(muon_src);
  consumes<pat::CompositeCandidateCollection>(dimuon_src);
  produces<pat::MuonCollection>("muons");
}

void MuonsForNminus1::produce(edm::Event& event, const edm::EventSetup& setup) {

  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  event.getByLabel(dimuon_src, dimuons);

  if(!dimuons.isValid())
    edm::LogWarning("DileptonHandleInvalid") << "tried to get " << dimuon_src << " and failed!";

  edm::Handle<pat::MuonCollection> pre_muons;
  event.getByLabel(muon_src, pre_muons);

  if(!pre_muons.isValid())
    edm::LogWarning("LeptonHandleInvalid") << "tried to get " << muon_src << " with edm::Handle<edm::View<reco::Candidate> > and failed!";

  std::unique_ptr<pat::MuonCollection> new_muons(new pat::MuonCollection);

  // if there is no dimoun pair passing full selection, store all the muons.
  if(dimuons->size() < 1) {
    for(pat::MuonCollection::const_iterator mu = pre_muons->begin(), mue = pre_muons->end(); mu != mue; ++mu) {
      pat::Muon* new_mu = mu->clone();
      new_muons->push_back(*new_mu);
      delete new_mu;
    }
  }

  // if there are dimuon pairs passing full selection, store only those muons.
  else {
    for(pat::CompositeCandidateCollection::const_iterator di = dimuons->begin(), die = dimuons->end(); di != die; ++di) {
      new_muons->push_back(toConcrete<pat::Muon>(dileptonDaughter(*di, 0)));
      new_muons->push_back(toConcrete<pat::Muon>(dileptonDaughter(*di, 1)));
    }
  }

  event.put(std::move(new_muons), "muons");
}

DEFINE_FWK_MODULE(MuonsForNminus1);
