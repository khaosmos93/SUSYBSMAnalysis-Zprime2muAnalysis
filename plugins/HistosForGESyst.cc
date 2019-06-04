#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GEScaleSyst.h"

using namespace std;

class Zprime2muHistosForGESyst : public edm::EDAnalyzer {
 public:
  explicit Zprime2muHistosForGESyst(const edm::ParameterSet&);
  ~Zprime2muHistosForGESyst();
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  edm::InputTag dilepton_src;

  double _eventWeight;
  bool   _useMadgraphWeight;
  double _madgraphWeight;

  TH1F *NDileptons;
  TH1F *WeightMadGraph;

  int n_copies;
  std::vector<int> copies;

  std::vector<TH1F *> h_DimuonMassBB;
  std::vector<TH1F *> h_DimuonMassBE;

  std::vector<TH1F *> h_MuonEta;
  std::vector<TH1F *> h_MuonPhi;
  std::vector<TH1F *> h_MuonPtB;
  std::vector<TH1F *> h_MuonPtE;

  std::vector<TH1F *> h_MuonPEta;
  std::vector<TH1F *> h_MuonPPhi;
  std::vector<TH1F *> h_MuonPPtB;
  std::vector<TH1F *> h_MuonPPtE;

  std::vector<TH1F *> h_MuonMEta;
  std::vector<TH1F *> h_MuonMPhi;
  std::vector<TH1F *> h_MuonMPtB;
  std::vector<TH1F *> h_MuonMPtE;

  std::vector<TProfile *> p_DimuonMassBB;
  std::vector<TProfile *> p_DimuonMassBE;

  std::vector<TProfile *> p_MuonEta;
  std::vector<TProfile *> p_MuonPhi;
  std::vector<TProfile *> p_MuonPtB;
  std::vector<TProfile *> p_MuonPtE;

  std::vector<TProfile *> p_MuonPEta;
  std::vector<TProfile *> p_MuonPPhi;
  std::vector<TProfile *> p_MuonPPtB;
  std::vector<TProfile *> p_MuonPPtE;

  std::vector<TProfile *> p_MuonMEta;
  std::vector<TProfile *> p_MuonMPhi;
  std::vector<TProfile *> p_MuonMPtB;
  std::vector<TProfile *> p_MuonMPtE;

  const bool ShutUp;
};

Zprime2muHistosForGESyst::Zprime2muHistosForGESyst(const edm::ParameterSet& cfg)
  : dilepton_src(cfg.getParameter<edm::InputTag>("dilepton_src")),

    _useMadgraphWeight(cfg.getParameter<bool>("useMadgraphWeight")),
    _madgraphWeight(1.),

    copies(cfg.getParameter<std::vector<int>>("copies")),

    ShutUp(cfg.getParameter<bool>("ShutUp"))
{
  consumes<pat::CompositeCandidateCollection>(dilepton_src);
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));

  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  n_copies = copies.size();

  // Dilepton multiplicity.
  NDileptons = fs->make<TH1F>("NDileptons", "# dileptons/event", 10, 0, 10);

  // Generator weights
  WeightMadGraph = fs->make<TH1F>("WeightMadGraph", "weight per event", 4, -2,2);

  int n_eta_bins = 6;
  double eta_bins[7] = { -2.4, -2.1, -1.2, 0.0, 1.2, 2.1, 2.4 };

  int n_phi_bins = 6;
  double phi_bins[7] = {
    -TMath::Pi()*(3./3.), -TMath::Pi()*(2./3.), -TMath::Pi()*(1./3.), 0.0,
     TMath::Pi()*(1./3.), TMath::Pi()*(2./3.), TMath::Pi()*(3./3.)
  };

  for(int i=0; i<n_copies; ++i) {
    TString tag = TString::Format("_%d", copies[i]);

    TProfile *_p_DimuonMassBB = fs->make<TProfile>("p_DimuonMassBB"+tag, "", 100, 0, 10000 );
    TProfile *_p_DimuonMassBE = fs->make<TProfile>("p_DimuonMassBE"+tag, "", 100, 0, 10000 );
    TProfile *_p_MuonEta      = fs->make<TProfile>("p_MuonEta"+tag,      "", n_eta_bins, eta_bins );
    TProfile *_p_MuonPhi      = fs->make<TProfile>("p_MuonPhi"+tag,      "", n_phi_bins, phi_bins );
    TProfile *_p_MuonPtB      = fs->make<TProfile>("p_MuonPtB"+tag,      "", 100, 0, 10000 );
    TProfile *_p_MuonPtE      = fs->make<TProfile>("p_MuonPtE"+tag,      "", 100, 0, 10000 );
    TProfile *_p_MuonPEta     = fs->make<TProfile>("p_MuonPEta"+tag,     "", n_eta_bins, eta_bins );
    TProfile *_p_MuonPPhi     = fs->make<TProfile>("p_MuonPPhi"+tag,     "", n_phi_bins, phi_bins );
    TProfile *_p_MuonPPtB     = fs->make<TProfile>("p_MuonPPtB"+tag,     "", 100, 0, 10000 );
    TProfile *_p_MuonPPtE     = fs->make<TProfile>("p_MuonPPtE"+tag,     "", 100, 0, 10000 );
    TProfile *_p_MuonMEta     = fs->make<TProfile>("p_MuonMEta"+tag,     "", n_eta_bins, eta_bins );
    TProfile *_p_MuonMPhi     = fs->make<TProfile>("p_MuonMPhi"+tag,     "", n_phi_bins, phi_bins );
    TProfile *_p_MuonMPtB     = fs->make<TProfile>("p_MuonMPtB"+tag,     "", 100, 0, 10000 );
    TProfile *_p_MuonMPtE     = fs->make<TProfile>("p_MuonMPtE"+tag,     "", 100, 0, 10000 );

    p_DimuonMassBB.push_back(_p_DimuonMassBB);
    p_DimuonMassBE.push_back(_p_DimuonMassBE);
    p_MuonEta.push_back(_p_MuonEta);
    p_MuonPhi.push_back(_p_MuonPhi);
    p_MuonPtB.push_back(_p_MuonPtB);
    p_MuonPtE.push_back(_p_MuonPtE);
    p_MuonPEta.push_back(_p_MuonPEta);
    p_MuonPPhi.push_back(_p_MuonPPhi);
    p_MuonPPtB.push_back(_p_MuonPPtB);
    p_MuonPPtE.push_back(_p_MuonPPtE);
    p_MuonMEta.push_back(_p_MuonMEta);
    p_MuonMPhi.push_back(_p_MuonMPhi);
    p_MuonMPtB.push_back(_p_MuonMPtB);
    p_MuonMPtE.push_back(_p_MuonMPtE);


    TH1F *_h_DimuonMassBB = fs->make<TH1F>("h_DimuonMassBB"+tag, "", 100, 0, 10000 );
    TH1F *_h_DimuonMassBE = fs->make<TH1F>("h_DimuonMassBE"+tag, "", 100, 0, 10000 );
    TH1F *_h_MuonEta      = fs->make<TH1F>("h_MuonEta"+tag,      "", n_eta_bins, eta_bins );
    TH1F *_h_MuonPhi      = fs->make<TH1F>("h_MuonPhi"+tag,      "", n_phi_bins, phi_bins );
    TH1F *_h_MuonPtB      = fs->make<TH1F>("h_MuonPtB"+tag,      "", 100, 0, 10000 );
    TH1F *_h_MuonPtE      = fs->make<TH1F>("h_MuonPtE"+tag,      "", 100, 0, 10000 );
    TH1F *_h_MuonPEta     = fs->make<TH1F>("h_MuonPEta"+tag,     "", n_eta_bins, eta_bins );
    TH1F *_h_MuonPPhi     = fs->make<TH1F>("h_MuonPPhi"+tag,     "", n_phi_bins, phi_bins );
    TH1F *_h_MuonPPtB     = fs->make<TH1F>("h_MuonPPtB"+tag,     "", 100, 0, 10000 );
    TH1F *_h_MuonPPtE     = fs->make<TH1F>("h_MuonPPtE"+tag,     "", 100, 0, 10000 );
    TH1F *_h_MuonMEta     = fs->make<TH1F>("h_MuonMEta"+tag,     "", n_eta_bins, eta_bins );
    TH1F *_h_MuonMPhi     = fs->make<TH1F>("h_MuonMPhi"+tag,     "", n_phi_bins, phi_bins );
    TH1F *_h_MuonMPtB     = fs->make<TH1F>("h_MuonMPtB"+tag,     "", 100, 0, 10000 );
    TH1F *_h_MuonMPtE     = fs->make<TH1F>("h_MuonMPtE"+tag,     "", 100, 0, 10000 );

    h_DimuonMassBB.push_back(_h_DimuonMassBB);
    h_DimuonMassBE.push_back(_h_DimuonMassBE);
    h_MuonEta.push_back(_h_MuonEta);
    h_MuonPhi.push_back(_h_MuonPhi);
    h_MuonPtB.push_back(_h_MuonPtB);
    h_MuonPtE.push_back(_h_MuonPtE);
    h_MuonPEta.push_back(_h_MuonPEta);
    h_MuonPPhi.push_back(_h_MuonPPhi);
    h_MuonPPtB.push_back(_h_MuonPPtB);
    h_MuonPPtE.push_back(_h_MuonPPtE);
    h_MuonMEta.push_back(_h_MuonMEta);
    h_MuonMPhi.push_back(_h_MuonMPhi);
    h_MuonMPtB.push_back(_h_MuonMPtB);
    h_MuonMPtE.push_back(_h_MuonMPtE);
  }

}

Zprime2muHistosForGESyst::~Zprime2muHistosForGESyst() {}

void Zprime2muHistosForGESyst::analyze(const edm::Event& event, const edm::EventSetup& setup) {

  //---- Generator weights
  if (_useMadgraphWeight) {
    _eventWeight = 1.;
    _madgraphWeight = 1.;

    edm::Handle<GenEventInfoProduct> gen_ev_info;
    event.getByLabel(edm::InputTag("generator"), gen_ev_info);
    if (gen_ev_info.isValid()){
      _eventWeight = gen_ev_info->weight();
      _madgraphWeight = ( _eventWeight > 0 ) ? 1.0 : -1.0;
    }
    WeightMadGraph->Fill( _madgraphWeight );
  }

  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  event.getByLabel(dilepton_src, dileptons);

  if( !dileptons.isValid() ) {
    std::cout << "Zprime2muHistosForGESyst::analyze : !dileptons.isValid() ---> return" << std::endl;
    return;
  }

  NDileptons->Fill( dileptons->size(), _madgraphWeight );

  std::unique_ptr<GEScaleSyst> GE( new GEScaleSyst() );
  if(!ShutUp)  GE->SetVerbose(1);

  pat::CompositeCandidateCollection::const_iterator dil = dileptons->begin(), dile = dileptons->end();
  for ( ; dil != dile; ++dil) {

    double mass = dil->mass();

    const reco::CandidateBaseRef& lep0 = dileptonDaughter(*dil, 0);
    const reco::CandidateBaseRef& lep1 = dileptonDaughter(*dil, 1);
    if( !lep0.isNonnull() || !lep1.isNonnull())  continue;

    const reco::CandidateBaseRef& lepP = lep0->charge() > 0 ? lep0 : lep1;
    const reco::CandidateBaseRef& lepM = lep0->charge() > 0 ? lep1 : lep0;

    bool isbb = ( fabs(lepP->eta()) < 1.2 && fabs(lepM->eta()) < 1.2 );

    for(int i=0; i<n_copies; ++i) {
      TLorentzVector LvecP = GE->GEScaleCorrLvec( copies[i], lepP->pt(), lepP->eta(), lepP->phi(), lepP->charge() );
      TLorentzVector LvecM = GE->GEScaleCorrLvec( copies[i], lepM->pt(), lepM->eta(), lepM->phi(), lepM->charge() );

      double mass_corr = (LvecP+LvecM).M();

      if(copies[i] == 0) {
        double diff = (mass - mass_corr);
        if(diff > 1e-3) {
          cout << "WARNING: mass - mass_corr > 1e-3" << endl;
          cout << "         \t mass =     " << mass << endl;
          cout << "         \t mass_corr  = " << mass_corr << endl;
        }
      }

      //-- fill mass
      if(isbb) {
        h_DimuonMassBB[i]->Fill( mass_corr, _madgraphWeight );
        p_DimuonMassBB[i]->Fill( mass, mass_corr, _madgraphWeight );
      }
      else {
        h_DimuonMassBE[i]->Fill( mass_corr, _madgraphWeight );
        p_DimuonMassBE[i]->Fill( mass, mass_corr, _madgraphWeight );
      }

      //-- fill sigle muon variables
      h_MuonEta[i]->Fill(  LvecP.Eta(), _madgraphWeight );
      h_MuonEta[i]->Fill(  LvecM.Eta(), _madgraphWeight );
      h_MuonPEta[i]->Fill( LvecP.Eta(), _madgraphWeight );
      h_MuonMEta[i]->Fill( LvecM.Eta(), _madgraphWeight );
      h_MuonPhi[i]->Fill(  LvecP.Phi(), _madgraphWeight );
      h_MuonPhi[i]->Fill(  LvecM.Phi(), _madgraphWeight );
      h_MuonPPhi[i]->Fill( LvecP.Phi(), _madgraphWeight );
      h_MuonMPhi[i]->Fill( LvecM.Phi(), _madgraphWeight );

      p_MuonEta[i]->Fill(  lepP->eta(), LvecP.Pt(), _madgraphWeight );
      p_MuonEta[i]->Fill(  lepM->eta(), LvecM.Pt(), _madgraphWeight );
      p_MuonPEta[i]->Fill( lepP->eta(), LvecP.Pt(), _madgraphWeight );
      p_MuonMEta[i]->Fill( lepM->eta(), LvecM.Pt(), _madgraphWeight );
      p_MuonPhi[i]->Fill(  lepP->phi(), LvecP.Pt(), _madgraphWeight );
      p_MuonPhi[i]->Fill(  lepM->phi(), LvecM.Pt(), _madgraphWeight );
      p_MuonPPhi[i]->Fill( lepP->phi(), LvecP.Pt(), _madgraphWeight );
      p_MuonMPhi[i]->Fill( lepM->phi(), LvecM.Pt(), _madgraphWeight );

      if( fabs(lepP->eta()) < 1.2 ) {
        h_MuonPtB[i]->Fill(  LvecP.Pt(), _madgraphWeight );
        h_MuonPPtB[i]->Fill( LvecP.Pt(), _madgraphWeight );
        p_MuonPtB[i]->Fill(  lepP->pt(), LvecP.Pt(), _madgraphWeight );
        p_MuonPPtB[i]->Fill( lepP->pt(), LvecP.Pt(), _madgraphWeight );
      }
      else {
        h_MuonPtE[i]->Fill(  LvecP.Pt(), _madgraphWeight );
        h_MuonPPtE[i]->Fill( LvecP.Pt(), _madgraphWeight );
        p_MuonPtE[i]->Fill(  lepP->pt(), LvecP.Pt(), _madgraphWeight );
        p_MuonPPtE[i]->Fill( lepP->pt(), LvecP.Pt(), _madgraphWeight );
      }

      if( fabs(lepM->eta()) < 1.2 ) {
        h_MuonPtB[i]->Fill(  LvecM.Pt(), _madgraphWeight );
        h_MuonMPtB[i]->Fill( LvecM.Pt(), _madgraphWeight );
        p_MuonPtB[i]->Fill(  lepM->pt(), LvecM.Pt(), _madgraphWeight );
        p_MuonMPtB[i]->Fill( lepM->pt(), LvecM.Pt(), _madgraphWeight );
      }
      else {
        h_MuonPtE[i]->Fill(  LvecM.Pt(), _madgraphWeight );
        h_MuonMPtE[i]->Fill( LvecM.Pt(), _madgraphWeight );
        p_MuonPtE[i]->Fill(  lepM->pt(), LvecM.Pt(), _madgraphWeight );
        p_MuonMPtE[i]->Fill( lepM->pt(), LvecM.Pt(), _madgraphWeight );
      }

    }

  } // for ( ; dil != dile; ++dil)

}

DEFINE_FWK_MODULE(Zprime2muHistosForGESyst);
