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

//using namespace std;

class Zprime2muHistosForDimuTrigEff : public edm::EDAnalyzer {
 public:
  explicit Zprime2muHistosForDimuTrigEff(const edm::ParameterSet&);
  ~Zprime2muHistosForDimuTrigEff();
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  void getBSandPV(const edm::Event&);
  
  double getSingleMuEff(const reco::CandidateBaseRef&);
  bool triggerPassed(const reco::CandidateBaseRef&, const reco::CandidateBaseRef&);
  bool triggerPassedSS(const reco::CandidateBaseRef&, const reco::CandidateBaseRef&);

  edm::InputTag dilepton_src;
  edm::InputTag beamspot_src;
  edm::InputTag vertex_src;
  const bool use_bs_and_pv;
  const reco::BeamSpot* beamspot;
  const reco::Vertex*   vertex;
  int                   nVtx;

  std::string effType;

  bool   _usePrescaleWeight;
  int    _prescaleWeight;
  double _eventWeight;
  bool   _useMadgraphWeight;
  double _madgraphWeight;

  HardInteraction* hardInteraction;

  double pt0;
  double pt1;
  double pt2;
  double pt3;
  double pt4;
  double pt5;

  double eta0;
  double eta1;
  double eta2;
  double eta3;
  double eta4;

  std::vector< std::vector<double> > muTrigEffDATA2017;
  std::vector< std::vector<double> > muTrigEffMC2017;

  std::vector< std::vector<double> > muTrigEffDATA2018;
  std::vector< std::vector<double> > muTrigEffMC2018;

  std::vector< std::vector<double> > muTrigEffTEST99;
  std::vector< std::vector<double> > muTrigEffTEST80;
  std::vector< std::vector<double> > muTrigEffTEST50;

  TH1F *NBeamSpot;
  TH1F *NVertices;
  TH1F *NDileptons;
  TH1F *WeightMadGraph;

  TH1F *DimuonGenMass;
  TH1F *DimuonGenMassBB;
  TH1F *DimuonGenMassBE;

  TH1F *DimuonGenMassPass;
  TH1F *DimuonGenMassPassBB;
  TH1F *DimuonGenMassPassBE;

  TH1F *DimuonGenMassPassSS;
  TH1F *DimuonGenMassPassSSBB;
  TH1F *DimuonGenMassPassSSBE;

  const bool ShutUp;

  int nTotal;
  int nTotal_BB;
  int nTotal_BE;
  int nPass;
  int nPass_BB;
  int nPass_BE;
  int nPassSS;
  int nPassSS_BB;
  int nPassSS_BE;
};

Zprime2muHistosForDimuTrigEff::Zprime2muHistosForDimuTrigEff(const edm::ParameterSet& cfg)
  : dilepton_src(cfg.getParameter<edm::InputTag>("dilepton_src")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
    use_bs_and_pv(cfg.getParameter<bool>("use_bs_and_pv")),
    beamspot(0),
    vertex(0),
    nVtx(0),

    effType(cfg.getParameter<std::string>("EffType")),

    _usePrescaleWeight(cfg.getUntrackedParameter<bool>("usePrescaleWeight",false)),
    _prescaleWeight(1),
    _useMadgraphWeight(cfg.getParameter<bool>("useMadgraphWeight")),
    _madgraphWeight(1.),

    hardInteraction( new HardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")) ),

    ShutUp(cfg.getParameter<bool>("ShutUp")),

    nTotal(0),
    nTotal_BB(0),
    nTotal_BE(0),
    nPass(0),
    nPass_BB(0),
    nPass_BE(0),
    nPassSS(0),
    nPassSS_BB(0),
    nPassSS_BE(0)
{
  consumes<pat::CompositeCandidateCollection>(dilepton_src);
  consumes<reco::BeamSpot>(beamspot_src);
  consumes<reco::VertexCollection>(vertex_src);
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));
  consumes<std::vector<reco::GenParticle>>(hardInteraction->src);

  pt0 = 52.;
  pt1 = 55.;
  pt2 = 60.;
  pt3 = 120.;
  pt4 = 200.;
  pt5 = 300.;

  eta0 = 0.;
  eta1 = 0.9;
  eta2 = 1.2;
  eta3 = 2.1;
  eta4 = 2.4;

  muTrigEffDATA2017 = {
    {  0.921714, 0.903086, 0.841379, 0.646483  },
    {  0.926002, 0.908560, 0.871136, 0.717869  },
    {  0.925014, 0.910824, 0.895099, 0.797875  },
    {  0.923994, 0.916634, 0.913087, 0.843658  },
    {  0.913003, 0.917335, 0.908076, 0.867220  },
    {  0.907988, 0.895765, 0.914186, 0.785714  }
  };

  muTrigEffMC2017 = {
    {  0.952633, 0.972206, 0.858955, 0.743717  },
    {  0.956621, 0.953841, 0.887878, 0.830402  },
    {  0.954880, 0.964125, 0.908754, 0.861802  },
    {  0.955035, 0.968954, 0.921182, 0.850175  },
    {  0.946782, 0.967631, 0.921254, 0.847839  },
    {  0.942490, 0.954287, 0.912672, 0.868006  }
  };

  muTrigEffDATA2018 = {
    {  0.934448, 0.941939, 0.905695, 0.785846  },
    {  0.940030, 0.943628, 0.914531, 0.819836  },
    {  0.935614, 0.940465, 0.917430, 0.839046  },
    {  0.935093, 0.938657, 0.921874, 0.844658  },
    {  0.924876, 0.916621, 0.908925, 0.844193  },
    {  0.913842, 0.940379, 0.904000, 0.831776  }
  };

  muTrigEffMC2018 = {
    {  0.953250, 0.959069, 0.921012, 0.796532  },
    {  0.957665, 0.976723, 0.903906, 0.861013  },
    {  0.959063, 0.972395, 0.915370, 0.840214  },
    {  0.954010, 0.973465, 0.918293, 0.842394  },
    {  0.945895, 0.954395, 0.901103, 0.824672  },
    {  0.939026, 0.954834, 0.910572, 0.837956  }
  };

  muTrigEffTEST99 = { //abseta:[0.0,    0.9,      1.2,      2.1,      2.4]
    { 0.99, 0.99, 0.99, 0.99 }, //pt:[52.0,55.0]
    { 0.99, 0.99, 0.99, 0.99 }, //pt:[55.0,60.0] 
    { 0.99, 0.99, 0.99, 0.99 }, //pt:[60.0,120.0]
    { 0.99, 0.99, 0.99, 0.99 }, //pt:[120.0,200.0]
    { 0.99, 0.99, 0.99, 0.99 }, //pt:[200.0,300.0]
    { 0.99, 0.99, 0.99, 0.99 }  //pt:[300.0,]
  };

  muTrigEffTEST80 = { //abseta:[0.0,    0.9,      1.2,      2.1,      2.4]
    { 0.80, 0.80, 0.80, 0.80 }, //pt:[52.0,55.0]
    { 0.80, 0.80, 0.80, 0.80 }, //pt:[55.0,60.0]
    { 0.80, 0.80, 0.80, 0.80 }, //pt:[60.0,120.0]
    { 0.80, 0.80, 0.80, 0.80 }, //pt:[120.0,200.0]
    { 0.80, 0.80, 0.80, 0.80 }, //pt:[200.0,300.0]
    { 0.80, 0.80, 0.80, 0.80 }  //pt:[300.0,]
  };

  muTrigEffTEST50 = { //abseta:[0.0,    0.9,      1.2,      2.1,      2.4]
    { 0.50, 0.50, 0.50, 0.50 }, //pt:[52.0,55.0]
    { 0.50, 0.50, 0.50, 0.50 }, //pt:[55.0,60.0]
    { 0.50, 0.50, 0.50, 0.50 }, //pt:[60.0,120.0]
    { 0.50, 0.50, 0.50, 0.50 }, //pt:[120.0,200.0]
    { 0.50, 0.50, 0.50, 0.50 }, //pt:[200.0,300.0]
    { 0.50, 0.50, 0.50, 0.50 }  //pt:[300.0,]
  };

  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);

  NBeamSpot = fs->make<TH1F>("NBeamSpot", "# beamspots/event",  2, 0,  2);
  NVertices = fs->make<TH1F>("NVertices", "# vertices/event",  200, 0, 200);

  // Dilepton multiplicity.
  NDileptons = fs->make<TH1F>("NDileptons", "# dileptons/event", 10, 0, 10);

  // Generator weights
  WeightMadGraph = fs->make<TH1F>("WeightMadGraph", "weight per event", 4, -2,2);

  DimuonGenMass = fs->make<TH1F>("DimuonGenMass", "", 10000, 0, 10000);
  DimuonGenMassBB = fs->make<TH1F>("DimuonGenMassBB", "", 10000, 0, 10000);
  DimuonGenMassBE = fs->make<TH1F>("DimuonGenMassBE", "", 10000, 0, 10000);

  DimuonGenMassPass = fs->make<TH1F>("DimuonGenMassPass", "", 10000, 0, 10000);
  DimuonGenMassPassBB = fs->make<TH1F>("DimuonGenMassPassBB", "", 10000, 0, 10000);
  DimuonGenMassPassBE = fs->make<TH1F>("DimuonGenMassPassBE", "", 10000, 0, 10000);

  DimuonGenMassPassSS = fs->make<TH1F>("DimuonGenMassPassSS", "", 10000, 0, 10000);
  DimuonGenMassPassSSBB = fs->make<TH1F>("DimuonGenMassPassSSBB", "", 10000, 0, 10000);
  DimuonGenMassPassSSBE = fs->make<TH1F>("DimuonGenMassPassSSBE", "", 10000, 0, 10000);
}

Zprime2muHistosForDimuTrigEff::~Zprime2muHistosForDimuTrigEff()
{
  std::cout << "\n\n";

  std::cout << " nTotal : " << nTotal << std::endl;
  std::cout << "  nPass : " << nPass << std::endl;
  std::cout << "      --> " << double(nPass)/double(nTotal) << std::endl;
  std::cout << "nPassSS : " << nPassSS << std::endl;
  std::cout << "      --> " << double(nPassSS)/double(nTotal) << std::endl;

  std::cout << " nTotal_BB : " << nTotal_BB << std::endl;
  std::cout << "  nPass_BB : " << nPass_BB << std::endl;
  std::cout << "      --> " << double(nPass_BB)/double(nTotal_BB) << std::endl;
  std::cout << "nPassSS_BB : " << nPassSS_BB << std::endl;
  std::cout << "      --> " << double(nPassSS_BB)/double(nTotal_BB) << std::endl;

  std::cout << " nTotal_BE : " << nTotal_BE << std::endl;
  std::cout << "  nPass_BE : " << nPass_BE << std::endl;
  std::cout << "      --> " << double(nPass_BE)/double(nTotal_BE) << std::endl;
  std::cout << "nPassSS_BE : " << nPassSS_BE << std::endl;
  std::cout << "      --> " << double(nPassSS_BE)/double(nTotal_BE) << std::endl;

  std::cout << "\n\n";
}

void Zprime2muHistosForDimuTrigEff::getBSandPV(const edm::Event& event) {
  // We store these as bare pointers. Should find better way, but
  // don't want to pass them around everywhere...
  edm::Handle<reco::BeamSpot> hbs;
  event.getByLabel(beamspot_src, hbs);
  beamspot = hbs.isValid() ? &*hbs : 0; // nice and fragile
  NBeamSpot->Fill(beamspot != 0);

  edm::Handle<reco::VertexCollection> vertices;
  event.getByLabel(vertex_src, vertices);
  vertex = 0;
  int vertex_count = 0;
  for (reco::VertexCollection::const_iterator it = vertices->begin(), ite = vertices->end(); it != ite; ++it) {
    if (it->ndof() > 4 && fabs(it->z()) <= 24 && fabs(it->position().rho()) <= 2) {
      if (vertex == 0)
        vertex = &*it;
      ++vertex_count;
    }
  }
  nVtx = vertex_count;
  NVertices->Fill(vertex_count, _madgraphWeight );
}

double Zprime2muHistosForDimuTrigEff::getSingleMuEff(const reco::CandidateBaseRef& lep) {

  double eff = -1.;

  double pt = lep->pt();
  double abseta = fabs(lep->eta());

  std::vector< std::vector<double> > muTrigEff;

  if(effType == "DATA2017")         muTrigEff = muTrigEffDATA2017;
  else if(effType == "MC2017")      muTrigEff = muTrigEffMC2017;
  else if(effType == "DATA2018")    muTrigEff = muTrigEffDATA2018;
  else if(effType == "MC2018")      muTrigEff = muTrigEffMC2018;
  else if(effType == "TEST99")      muTrigEff = muTrigEffTEST99;
  else if(effType == "TEST80")      muTrigEff = muTrigEffTEST80;
  else if(effType == "TEST50")      muTrigEff = muTrigEffTEST50;

  if( pt0 <= pt && pt < pt1 ) {
    if(eta0 <= abseta && abseta < eta1)  eff = muTrigEff[0][0];
    if(eta1 <= abseta && abseta < eta2)  eff = muTrigEff[0][1];
    if(eta2 <= abseta && abseta < eta3)  eff = muTrigEff[0][2];
    if(eta3 <= abseta && abseta < eta4)  eff = muTrigEff[0][3];
  }
  else if( pt1 <= pt && pt < pt2 ) {
    if(eta0 <= abseta && abseta < eta1)  eff = muTrigEff[1][0];
    if(eta1 <= abseta && abseta < eta2)  eff = muTrigEff[1][1];
    if(eta2 <= abseta && abseta < eta3)  eff = muTrigEff[1][2];
    if(eta3 <= abseta && abseta < eta4)  eff = muTrigEff[1][3];
  }
  else if( pt2 <= pt && pt < pt3 ) {
    if(eta0 <= abseta && abseta < eta1)  eff = muTrigEff[2][0];
    if(eta1 <= abseta && abseta < eta2)  eff = muTrigEff[2][1];
    if(eta2 <= abseta && abseta < eta3)  eff = muTrigEff[2][2];
    if(eta3 <= abseta && abseta < eta4)  eff = muTrigEff[2][3];
  }
  else if( pt3 <= pt && pt < pt4 ) {
    if(eta0 <= abseta && abseta < eta1)  eff = muTrigEff[3][0];
    if(eta1 <= abseta && abseta < eta2)  eff = muTrigEff[3][1];
    if(eta2 <= abseta && abseta < eta3)  eff = muTrigEff[3][2];
    if(eta3 <= abseta && abseta < eta4)  eff = muTrigEff[3][3];
  }
  else if( pt4 <= pt && pt < pt5 ) {
    if(eta0 <= abseta && abseta < eta1)  eff = muTrigEff[4][0];
    if(eta1 <= abseta && abseta < eta2)  eff = muTrigEff[4][1];
    if(eta2 <= abseta && abseta < eta3)  eff = muTrigEff[4][2];
    if(eta3 <= abseta && abseta < eta4)  eff = muTrigEff[4][3];
  }
  else if( pt5 <= pt ) {
    if(eta0 <= abseta && abseta < eta1)  eff = muTrigEff[5][0];
    if(eta1 <= abseta && abseta < eta2)  eff = muTrigEff[5][1];
    if(eta2 <= abseta && abseta < eta3)  eff = muTrigEff[5][2];
    if(eta3 <= abseta && abseta < eta4)  eff = muTrigEff[5][3];
  }

  return eff;
}

bool Zprime2muHistosForDimuTrigEff::triggerPassed(const reco::CandidateBaseRef& lep0, const reco::CandidateBaseRef& lep1) {

  double eff0 = getSingleMuEff(lep0);
  double eff1 = getSingleMuEff(lep1);
  double dimuEff = 1. - (1. - eff0) * (1. - eff1);

  TRandom3 *r = new TRandom3(0);
  double rand_01 = r->Rndm();
  delete r;

  if(!ShutUp) {
    std::cout << "Zprime2muHistosForDimuTrigEff::triggerPassed" << std::endl;
    std::cout << "\tmu0 : pT=" << lep0->pt() << "\t|eta|=" << fabs(lep0->eta()) << std::endl;
    std::cout << "\t\teff=" << eff0 << std::endl;
    std::cout << "\tmu1 : pT=" << lep1->pt() << "\t|eta|=" << fabs(lep1->eta()) << std::endl;
    std::cout << "\t\teff=" << eff1 << std::endl;
    std::cout << "\t==>  Dimuon eff=" << dimuEff << "\trand=" << rand_01 << "\tisPass=" << (dimuEff > rand_01) << std::endl;
  }

  return (dimuEff > rand_01);
}

bool Zprime2muHistosForDimuTrigEff::triggerPassedSS(const reco::CandidateBaseRef& lep0, const reco::CandidateBaseRef& lep1) {

  double eff0 = getSingleMuEff(lep0);
  double eff1 = getSingleMuEff(lep1);

  TRandom3 *r1 = new TRandom3(0);
  Double_t rndm[2]; r1->RndmArray(2, rndm);
  Double_t rand0_01 = rndm[0];
  Double_t rand1_01 = rndm[1];
  delete r1;

  bool isLep0Passed = eff0 > rand0_01;
  bool isLep1Passed = eff1 > rand1_01;

  if(!ShutUp) {
    std::cout << "Zprime2muHistosForDimuTrigEff::triggerPassedSS" << std::endl;
    std::cout << "\tmu0 : pT=" << lep0->pt() << "\t|eta|=" << fabs(lep0->eta()) << std::endl;
    std::cout << "\t\teff=" << eff0 << "\trand=" << rand0_01 << "\tisPass=" << isLep0Passed << std::endl;
    std::cout << "\tmu1 : pT=" << lep1->pt() << "\t|eta|=" << fabs(lep1->eta()) << std::endl;
    std::cout << "\t\teff=" << eff1 << "\trand=" << rand1_01 << "\tisPass=" << isLep1Passed << std::endl;
  }

  return (isLep0Passed || isLep1Passed);
}

void Zprime2muHistosForDimuTrigEff::analyze(const edm::Event& event, const edm::EventSetup& setup) {

  //---- Prescales : Not using for now...
    //  edm::Handle<int> hltPrescale;
    //  edm::Handle<int> l1Prescale;
    //  event.getByLabel(edm::InputTag("getPrescales","HLTPrescale","Zprime2muAnalysis"), hltPrescale);
    //  event.getByLabel(edm::InputTag("getPrescales","L1Prescale","Zprime2muAnalysis"), l1Prescale);
    if (_usePrescaleWeight) {
      edm::Handle<int> totalPrescale;
      event.getByLabel(edm::InputTag("getPrescales","TotalPrescale","Zprime2muAnalysis"), totalPrescale);
      _prescaleWeight = *totalPrescale;
    }
    //  std::cout<<*hltPrescale<<std::endl;
    //  std::cout<<l1Prescale<<std::endl;
    //  std::cout<<totalPrescale<<std::endl;

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

  if (use_bs_and_pv)
    getBSandPV(event);

  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  event.getByLabel(dilepton_src, dileptons);

  if( !dileptons.isValid() ) {
    std::cout << "Zprime2muHistosForDimuTrigEff::analyze : !dileptons.isValid() ---> return" << std::endl;
    return;
  }

  //-- Gen-level mass
  double gM = -1.;
  hardInteraction->Fill(event);
  if( hardInteraction->IsValid() ) {
    reco::Particle::LorentzVector gen_dil = hardInteraction->dilepton();
    gM = gen_dil.mass();
  }

  NDileptons->Fill( dileptons->size(), _madgraphWeight );

  pat::CompositeCandidateCollection::const_iterator dil = dileptons->begin(), dile = dileptons->end();
  for ( ; dil != dile; ++dil) {

    const reco::CandidateBaseRef& lep0 = dileptonDaughter(*dil, 0);
    const reco::CandidateBaseRef& lep1 = dileptonDaughter(*dil, 1);
    if(lep0.isNonnull() && lep1.isNonnull()) {

      bool isBB = ( fabs(lep0->eta()) < 1.2 && fabs(lep1->eta()) < 1.2 );

      DimuonGenMass->Fill( gM, _madgraphWeight );
      nTotal += 1;
      if(isBB) {
        DimuonGenMassBB->Fill( gM, _madgraphWeight );
        nTotal_BB += 1;
      }
      else {
        DimuonGenMassBE->Fill( gM, _madgraphWeight );
        nTotal_BE += 1;
      }

      if( triggerPassed(lep0, lep1) ) {
        DimuonGenMassPass->Fill( gM, _madgraphWeight );
        nPass += 1;
        if(isBB) {
          DimuonGenMassPassBB->Fill( gM, _madgraphWeight );
          nPass_BB += 1;
        }
        else {
          DimuonGenMassPassBE->Fill( gM, _madgraphWeight );
          nPass_BE += 1;
        }
      }

      if( triggerPassedSS(lep0, lep1) ) {
        DimuonGenMassPassSS->Fill( gM, _madgraphWeight );
        nPassSS += 1;
        if(isBB) {
          DimuonGenMassPassSSBB->Fill( gM, _madgraphWeight );
          nPassSS_BB += 1;
        }
        else {
          DimuonGenMassPassSSBE->Fill( gM, _madgraphWeight );
          nPassSS_BE += 1;
        }
      }

    }


  } // for ( ; dil != dile; ++dil)

}

DEFINE_FWK_MODULE(Zprime2muHistosForDimuTrigEff);
