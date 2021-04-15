#include <boost/foreach.hpp>
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

class NtuplerLLJets : public edm::EDAnalyzer {
public:
  explicit NtuplerLLJets(const edm::ParameterSet&);
  ~NtuplerLLJets() { delete hardInteraction; }
  void analyze(const edm::Event&, const edm::EventSetup&);
  TString replace_all(const TString& a, const TString& b, const TString& c);

private:
  struct tree_t {
    unsigned run;
    unsigned lumi;
    unsigned long  event;
    float genWeight;

    float beamspot_x;
    float beamspot_x_err;
    float beamspot_y;
    float beamspot_y_err;
    float beamspot_z;
    float beamspot_z_err;

    int nvertices;

    float dil_mass;
    float dil_pt;
    float dil_rap;
    float dil_eta;
    float dil_phi;
    float dil_dR;
    float dil_dPhi;
    float dil_lep_pt[2];
    float cos_angle;
    float vertex_chi2;
    float cos_cs;
    float chi_dilepton;
    float phi_cs;

    float vertex_m;
    float vertex_m_err;
    float vertex_x;
    float vertex_x_err;
    float vertex_y;
    float vertex_y_err;
    float vertex_z;
    float vertex_z_err;

    int lep_id[2];
    float lep_p[2];
    float lep_pt[2];
    float lep_et[2];
    float lep_pt_err[2];
    float lep_px[2];
    float lep_py[2];
    float lep_pz[2];
    float lep_E[2];
    float lep_eta[2];
    float lep_phi[2];
    float lep_qOverPt[2];
    float lep_tk_p[2];
    float lep_tk_pt[2];
    float lep_tk_pt_err[2];
    float lep_tk_px[2];
    float lep_tk_py[2];
    float lep_tk_pz[2];
    float lep_tk_eta[2];
    float lep_tk_phi[2];
    float lep_tk_dz[2];
    float lep_tk_vz[2];
    float lep_tk_chi2[2];
    float lep_tk_ndf[2];
    float lep_tk_qOverPt[2];
    float lep_glb_p[2];
    float lep_glb_pt[2];
    float lep_glb_pt_err[2];
    float lep_glb_px[2];
    float lep_glb_py[2];
    float lep_glb_pz[2];
    float lep_glb_eta[2];
    float lep_glb_phi[2];
    float lep_glb_chi2[2];
    float lep_glb_ndf[2];
    float lep_glb_qOverPt[2];
    float lep_tpfms_p[2];
    float lep_tpfms_pt[2];
    float lep_tpfms_pt_err[2];
    float lep_tpfms_px[2];
    float lep_tpfms_py[2];
    float lep_tpfms_pz[2];
    float lep_tpfms_eta[2];
    float lep_tpfms_phi[2];
    float lep_tpfms_chi2[2];
    float lep_tpfms_ndf[2];
    float lep_tpfms_qOverPt[2];
    float lep_picky_p[2];
    float lep_picky_pt[2];
    float lep_picky_pt_err[2];
    float lep_picky_px[2];
    float lep_picky_py[2];
    float lep_picky_pz[2];
    float lep_picky_eta[2];
    float lep_picky_phi[2];
    float lep_picky_chi2[2];
    float lep_picky_ndf[2];
    float lep_picky_qOverPt[2];
    float lep_tuneP_p[2];
    float lep_tuneP_pt[2];
    float lep_tuneP_pt_err[2];
    float lep_tuneP_px[2];
    float lep_tuneP_py[2];
    float lep_tuneP_pz[2];
    float lep_tuneP_eta[2];
    float lep_tuneP_phi[2];
    float lep_tuneP_dz[2];
    float lep_tuneP_vz[2];
    float lep_tuneP_chi2[2];
    float lep_tuneP_ndf[2];
    float lep_tuneP_qOverPt[2];
    float lep_dyt_p[2];
    float lep_dyt_pt[2];
    float lep_dyt_pt_err[2];
    float lep_dyt_px[2];
    float lep_dyt_py[2];
    float lep_dyt_pz[2];
    float lep_dyt_eta[2];
    float lep_dyt_phi[2];
    float lep_dyt_dz[2];
    float lep_dyt_vz[2];
    float lep_dyt_chi2[2];
    float lep_dyt_ndf[2];
    float lep_dyt_qOverPt[2];

    float lep_Mu50_triggerMatchPt[2];
    float lep_Mu50_triggerMatchEta[2];
    float lep_Mu50_triggerMatchPhi[2];
    float lep_OldMu100_triggerMatchPt[2];
    float lep_OldMu100_triggerMatchEta[2];
    float lep_OldMu100_triggerMatchPhi[2];
    float lep_TkMu100_triggerMatchPt[2];
    float lep_TkMu100_triggerMatchEta[2];
    float lep_TkMu100_triggerMatchPhi[2];

    float lep_chi2dof[2];
    float lep_dB[2];
    float lep_sumPt[2];
    float lep_emEt[2];
    float lep_hadEt[2];
    float lep_hoEt[2];
    float lep_pfIso[2];
    float lep_pfIsoDB[2];
    int lep_timeNdof[2];
    float lep_timeInOut[2];
    float lep_timeOutIn[2];
    float lep_timeInOutErr[2];
    float lep_timeOutInErr[2];
    int lep_heep_id[2];
    float lep_gen_match[2];
    float lep_min_muon_dR[2];
    short lep_tk_numberOfValidTrackerHits[2];
    short lep_tk_numberOfValidTrackerLayers[2];
    short lep_tk_numberOfValidPixelHits[2];
    short lep_glb_numberOfValidTrackerHits[2]; 
    short lep_glb_numberOfValidTrackerLayers[2]; 
    short lep_glb_numberOfValidPixelHits[2];
    short lep_glb_numberOfValidMuonHits[2];
    short lep_glb_numberOfValidMuonDTHits[2];
    short lep_glb_numberOfValidMuonCSCHits[2];
    short lep_glb_numberOfValidMuonRPCHits[2];
    short lep_glb_muonStationsWithValidHits[2];
    short lep_glb_dtStationsWithValidHits[2];
    short lep_glb_cscStationsWithValidHits[2];
    short lep_glb_rpcStationsWithValidHits[2];
    short lep_glb_innermostMuonStationWithValidHits[2];
    short lep_glb_outermostMuonStationWithValidHits[2];
    short lep_tuneP_numberOfValidMuonHits[2];
    short lep_tuneP_numberOfValidMuonDTHits[2];
    short lep_tuneP_numberOfValidMuonCSCHits[2];
    short lep_tuneP_numberOfValidMuonRPCHits[2];
    short lep_tuneP_muonStationsWithValidHits[2];
    short lep_tuneP_dtStationsWithValidHits[2];
    short lep_tuneP_cscStationsWithValidHits[2];
    short lep_tuneP_rpcStationsWithValidHits[2];
    short lep_tuneP_innermostMuonStationWithValidHits[2];
    short lep_tuneP_outermostMuonStationWithValidHits[2];
    short lep_numberOfMatches[2];
    short lep_numberOfMatchedStations[2];
    short lep_numberOfMatchedRPCLayers[2];
    unsigned int lep_stationMask[2];
    int lep_numberOfChambers[2];
    int lep_numberOfChambersNoRPC[2];
    unsigned int lep_stationGapMaskDistance[2];
    unsigned int lep_stationGapMaskPull[2];
    bool lep_isGlobalMuon[2];
    bool lep_isTrackerMuon[2];

    bool GoodDataRan;
    bool GoodVtx;
    bool METFilter;

    float gen_res_mass;
    float gen_res_pt;
    float gen_res_rap;
    float gen_res_eta;
    float gen_res_phi;
    float gen_dil_mass;
    float gen_dil_pt;
    float gen_dil_rap;
    float gen_dil_eta;
    float gen_dil_phi;
    float gen_dil_dR;
    float gen_dil_dPhi;
    float gen_dil_noib_mass;
    float gen_dil_noib_pt;
    float gen_dil_noib_rap;
    float gen_dil_noib_eta;
    float gen_dil_noib_phi;
    float gen_dil_noib_dR;
    float gen_dil_noib_dPhi;
    float gen_lep_p[2];
    float gen_lep_pt[2];
    float gen_lep_px[2];
    float gen_lep_py[2];
    float gen_lep_pz[2];
    float gen_lep_E[2];
    float gen_lep_eta[2];
    float gen_lep_phi[2];
    float gen_lep_qOverPt[2];
    float gen_lep_noib_p[2];
    float gen_lep_noib_pt[2];
    float gen_lep_noib_px[2];
    float gen_lep_noib_py[2];
    float gen_lep_noib_pz[2];
    float gen_lep_noib_E[2];
    float gen_lep_noib_eta[2];
    float gen_lep_noib_phi[2];
    float gen_lep_noib_qOverPt[2];

    float met_pt;
    float met_phi;

    int nJets;
    std::vector<float> jet_pt;
    std::vector<float> jet_pt_Uncorrected;
    std::vector<float> jet_pt_L1FastJet;
    std::vector<float> jet_pt_L2Relative;
    std::vector<float> jet_pt_L3Absolute;
    std::vector<float> jet_pt_L2L3Residual;
    std::vector<float> jet_eta;
    std::vector<float> jet_phi;
    std::vector<float> jet_partonFlavour;
    std::vector<float> jet_hadronFlavour;
    std::vector<float> jet_NHF;
    std::vector<float> jet_NEMF;
    std::vector<float> jet_CHF;
    std::vector<float> jet_MUF;
    std::vector<float> jet_CEMF;
    std::vector<int> jet_NumConst;
    std::vector<int> jet_NumNeutralParticles;
    std::vector<int> jet_CHM;
    std::vector<float> jet_pfDeepCSVJetTags_probb;
    std::vector<float> jet_pfDeepCSVJetTags_probc;
    std::vector<float> jet_pfDeepCSVJetTags_probudsg;
    std::vector<float> jet_pfDeepCSVJetTags_probbb;
    std::vector<float> jet_pfDeepCSVJetTags_probcc;
  };

  tree_t t;
  TTree* tree;

  const edm::InputTag dimu_src;
  const edm::InputTag beamspot_src;
  const edm::InputTag met_src;
  const edm::InputTag jet_src;
  const edm::InputTag vertices_src;
  const bool fill_gen_info;
  const bool do_electrons;
  const edm::InputTag TriggerResults_src;
  const edm::InputTag genEventInfo_;
  std::vector<edm::InputTag> filterTags;
  HardInteraction* hardInteraction;
};

TString NtuplerLLJets::replace_all(const TString& a, const TString& b, const TString& c) {
  TString ret = a;
  ret.ReplaceAll(b, c);
  return ret;
}

NtuplerLLJets::NtuplerLLJets(const edm::ParameterSet& cfg)
  : dimu_src(cfg.getParameter<edm::InputTag>("dimu_src")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    met_src(cfg.getParameter<edm::InputTag>("met_src")),
    jet_src(cfg.getParameter<edm::InputTag>("jet_src")),
    vertices_src(cfg.getParameter<edm::InputTag>("vertices_src")),
    fill_gen_info(cfg.existsAs<edm::ParameterSet>("hardInteraction")),
    do_electrons(cfg.getParameter<bool>("doElectrons")),  
    TriggerResults_src(cfg.getParameter<edm::InputTag>("TriggerResults_src")),
    genEventInfo_(cfg.getUntrackedParameter<edm::InputTag>("genEventInfo")),
    filterTags(cfg.getParameter<std::vector<edm::InputTag> > ("metFilter")),  
    hardInteraction(fill_gen_info ? new HardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")) : 0)
{
    consumes<pat::CompositeCandidateCollection>(dimu_src);
    consumes<std::vector<pat::MET>>(met_src);
    consumes<std::vector<pat::Jet>>(jet_src);
    consumes<reco::BeamSpot>(beamspot_src);
    consumes<reco::VertexCollection>(vertices_src);
    consumes<edm::TriggerResults>(TriggerResults_src);
    consumes<GenEventInfoProduct>(genEventInfo_);
    if (fill_gen_info) consumes<std::vector<reco::GenParticle>>(hardInteraction->src);

    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("t", "");
    tree->Branch("run", &t.run, "run/i");
    tree->Branch("lumi", &t.lumi, "lumi/i");
    tree->Branch("event", &t.event, "event/l");
    tree->Branch("beamspot_x", &t.beamspot_x, "beamspot_x/F");
    tree->Branch("beamspot_x_err", &t.beamspot_x_err, "beamspot_x_err/F");
    tree->Branch("beamspot_y", &t.beamspot_y, "beamspot_y/F");
    tree->Branch("beamspot_y_err", &t.beamspot_y_err, "beamspot_y_err/F");
    tree->Branch("beamspot_z", &t.beamspot_z, "beamspot_z/F");
    tree->Branch("beamspot_z_err", &t.beamspot_z_err, "beamspot_z_err/F");
    tree->Branch("nvertices", &t.nvertices, "nvertices/I");
    tree->Branch("dil_mass", &t.dil_mass, "dil_mass/F");
    tree->Branch("dil_pt", &t.dil_pt, "dil_pt/F");
    tree->Branch("dil_rap", &t.dil_rap, "dil_rap/F");
    tree->Branch("dil_eta", &t.dil_eta, "dil_eta/F");
    tree->Branch("dil_phi", &t.dil_phi, "dil_phi/F");
    tree->Branch("dil_dR", &t.dil_dR, "dil_dR/F");
    tree->Branch("dil_dPhi", &t.dil_dPhi, "dil_dPhi/F");
    tree->Branch("dil_lep_pt", t.dil_lep_pt, "dil_lep_pt[2]/F");
    tree->Branch("cos_angle", &t.cos_angle, "cos_angle/F");
    tree->Branch("vertex_chi2", &t.vertex_chi2, "vertex_chi2/F");
    tree->Branch("cos_cs", &t.cos_cs, "cos_cs/F");
    tree->Branch("chi_dilepton", &t.chi_dilepton, "chi_dilepton/F");
    tree->Branch("phi_cs", &t.phi_cs, "phi_cs/F");
    tree->Branch("vertex_m", &t.vertex_m, "vertex_m/F");
    tree->Branch("vertex_m_err", &t.vertex_m_err, "vertex_m_err/F");
    tree->Branch("vertex_x", &t.vertex_x, "vertex_x/F");
    tree->Branch("vertex_x_err", &t.vertex_x_err, "vertex_x_err/F");
    tree->Branch("vertex_y", &t.vertex_y, "vertex_y/F");
    tree->Branch("vertex_y_err", &t.vertex_y_err, "vertex_y_err/F");
    tree->Branch("vertex_z", &t.vertex_z, "vertex_z/F");
    tree->Branch("vertex_z_err", &t.vertex_z_err, "vertex_z_err/F");
    tree->Branch("lep_id", t.lep_id, "lep_id[2]/I");
    tree->Branch("lep_heep_id", t.lep_heep_id, "lep_heep_id[2]/I");
    tree->Branch("lep_gen_match", t.lep_gen_match, "lep_gen_match[2]/F");
    tree->Branch("lep_p", t.lep_p, "lep_p[2]/F");
    tree->Branch("lep_pt", t.lep_pt, "lep_pt[2]/F");
    tree->Branch("lep_pt_err", t.lep_pt_err, "lep_pt_err[2]/F");
    tree->Branch("lep_px", t.lep_px, "lep_px[2]/F");
    tree->Branch("lep_py", t.lep_py, "lep_py[2]/F");
    tree->Branch("lep_pz", t.lep_pz, "lep_pz[2]/F");
    tree->Branch("lep_E", t.lep_E, "lep_E[2]/F");
    tree->Branch("lep_eta", t.lep_eta, "lep_eta[2]/F");
    tree->Branch("lep_et", t.lep_et, "lep_et[2]/F");
    tree->Branch("lep_phi", t.lep_phi, "lep_phi[2]/F");
    tree->Branch("lep_qOverPt", t.lep_qOverPt, "lep_qOverPt[2]/F");
    tree->Branch("lep_tk_p", t.lep_tk_p, "lep_tk_p[2]/F");
    tree->Branch("lep_tk_pt", t.lep_tk_pt, "lep_tk_pt[2]/F");
    tree->Branch("lep_tk_pt_err", t.lep_tk_pt_err, "lep_tk_pt_err[2]/F");
    tree->Branch("lep_tk_px", t.lep_tk_px, "lep_tk_px[2]/F");
    tree->Branch("lep_tk_py", t.lep_tk_py, "lep_tk_py[2]/F");
    tree->Branch("lep_tk_pz", t.lep_tk_pz, "lep_tk_pz[2]/F");
    tree->Branch("lep_tk_eta", t.lep_tk_eta, "lep_tk_eta[2]/F");
    tree->Branch("lep_tk_phi", t.lep_tk_phi, "lep_tk_phi[2]/F");
    tree->Branch("lep_tk_dz", t.lep_tk_dz, "lep_tk_dz[2]/F");
    tree->Branch("lep_tk_vz", t.lep_tk_vz, "lep_tk_vz[2]/F");
    tree->Branch("lep_tk_chi2", t.lep_tk_chi2, "lep_tk_chi2[2]/F");
    tree->Branch("lep_tk_ndf", t.lep_tk_ndf, "lep_tk_ndf[2]/F");
    tree->Branch("lep_tk_qOverPt", t.lep_tk_qOverPt, "lep_tk_qOverPt[2]/F");
    tree->Branch("lep_glb_p", t.lep_glb_p, "lep_glb_p[2]/F");
    tree->Branch("lep_glb_pt", t.lep_glb_pt, "lep_glb_pt[2]/F");
    tree->Branch("lep_glb_pt_err", t.lep_glb_pt_err, "lep_glb_pt_err[2]/F");
    tree->Branch("lep_glb_px", t.lep_glb_px, "lep_glb_px[2]/F");
    tree->Branch("lep_glb_py", t.lep_glb_py, "lep_glb_py[2]/F");
    tree->Branch("lep_glb_pz", t.lep_glb_pz, "lep_glb_pz[2]/F");
    tree->Branch("lep_glb_eta", t.lep_glb_eta, "lep_glb_eta[2]/F");
    tree->Branch("lep_glb_phi", t.lep_glb_phi, "lep_glb_phi[2]/F");
    tree->Branch("lep_glb_chi2", t.lep_glb_chi2, "lep_glb_chi2[2]/F");
    tree->Branch("lep_glb_ndf", t.lep_glb_ndf, "lep_glb_ndf[2]/F");
    tree->Branch("lep_glb_qOverPt", t.lep_glb_qOverPt, "lep_glb_qOverPt[2]/F");
    tree->Branch("lep_tpfms_p", t.lep_tpfms_p, "lep_tpfms_p[2]/F");
    tree->Branch("lep_tpfms_pt", t.lep_tpfms_pt, "lep_tpfms_pt[2]/F");
    tree->Branch("lep_tpfms_pt_err", t.lep_tpfms_pt_err, "lep_tpfms_pt_err[2]/F");
    tree->Branch("lep_tpfms_px", t.lep_tpfms_px, "lep_tpfms_px[2]/F");
    tree->Branch("lep_tpfms_py", t.lep_tpfms_py, "lep_tpfms_py[2]/F");
    tree->Branch("lep_tpfms_pz", t.lep_tpfms_pz, "lep_tpfms_pz[2]/F");
    tree->Branch("lep_tpfms_eta", t.lep_tpfms_eta, "lep_tpfms_eta[2]/F");
    tree->Branch("lep_tpfms_phi", t.lep_tpfms_phi, "lep_tpfms_phi[2]/F");
    tree->Branch("lep_tpfms_chi2", t.lep_tpfms_chi2, "lep_tpfms_chi2[2]/F");
    tree->Branch("lep_tpfms_ndf", t.lep_tpfms_ndf, "lep_tpfms_ndf[2]/F");
    tree->Branch("lep_tpfms_qOverPt", t.lep_tpfms_qOverPt, "lep_tpfms_qOverPt[2]/F");
    tree->Branch("lep_picky_p", t.lep_picky_p, "lep_picky_p[2]/F");
    tree->Branch("lep_picky_pt", t.lep_picky_pt, "lep_picky_pt[2]/F");
    tree->Branch("lep_picky_pt_err", t.lep_picky_pt_err, "lep_picky_pt_err[2]/F");
    tree->Branch("lep_picky_px", t.lep_picky_px, "lep_picky_px[2]/F");
    tree->Branch("lep_picky_py", t.lep_picky_py, "lep_picky_py[2]/F");
    tree->Branch("lep_picky_pz", t.lep_picky_pz, "lep_picky_pz[2]/F");
    tree->Branch("lep_picky_eta", t.lep_picky_eta, "lep_picky_eta[2]/F");
    tree->Branch("lep_picky_phi", t.lep_picky_phi, "lep_picky_phi[2]/F");
    tree->Branch("lep_picky_chi2", t.lep_picky_chi2, "lep_picky_chi2[2]/F");
    tree->Branch("lep_picky_ndf", t.lep_picky_ndf, "lep_picky_ndf[2]/F");
    tree->Branch("lep_picky_qOverPt", t.lep_picky_qOverPt, "lep_picky_qOverPt[2]/F");
    tree->Branch("lep_tuneP_p", t.lep_tuneP_p, "lep_tuneP_p[2]/F");
    tree->Branch("lep_tuneP_pt", t.lep_tuneP_pt, "lep_tuneP_pt[2]/F");
    tree->Branch("lep_tuneP_pt_err", t.lep_tuneP_pt_err, "lep_tuneP_pt_err[2]/F");
    tree->Branch("lep_tuneP_px", t.lep_tuneP_px, "lep_tuneP_px[2]/F");
    tree->Branch("lep_tuneP_py", t.lep_tuneP_py, "lep_tuneP_py[2]/F");
    tree->Branch("lep_tuneP_pz", t.lep_tuneP_pz, "lep_tuneP_pz[2]/F");
    tree->Branch("lep_tuneP_eta", t.lep_tuneP_eta, "lep_tuneP_eta[2]/F");
    tree->Branch("lep_tuneP_phi", t.lep_tuneP_phi, "lep_tuneP_phi[2]/F");
    tree->Branch("lep_tuneP_chi2", t.lep_tuneP_chi2, "lep_tuneP_chi2[2]/F");
    tree->Branch("lep_tuneP_ndf", t.lep_tuneP_ndf, "lep_tuneP_ndf[2]/F");
    tree->Branch("lep_tuneP_qOverPt", t.lep_tuneP_qOverPt, "lep_tuneP_qOverPt[2]/F");
    tree->Branch("lep_dyt_p", t.lep_dyt_p, "lep_dyt_p[2]/F");
    tree->Branch("lep_dyt_pt", t.lep_dyt_pt, "lep_dyt_pt[2]/F");
    tree->Branch("lep_dyt_pt_err", t.lep_dyt_pt_err, "lep_dyt_pt_err[2]/F");
    tree->Branch("lep_dyt_px", t.lep_dyt_px, "lep_dyt_px[2]/F");
    tree->Branch("lep_dyt_py", t.lep_dyt_py, "lep_dyt_py[2]/F");
    tree->Branch("lep_dyt_pz", t.lep_dyt_pz, "lep_dyt_pz[2]/F");
    tree->Branch("lep_dyt_eta", t.lep_dyt_eta, "lep_dyt_eta[2]/F");
    tree->Branch("lep_dyt_phi", t.lep_dyt_phi, "lep_dyt_phi[2]/F");
    tree->Branch("lep_dyt_chi2", t.lep_dyt_chi2, "lep_dyt_chi2[2]/F");
    tree->Branch("lep_dyt_ndf", t.lep_dyt_ndf, "lep_dyt_ndf[2]/F");
    tree->Branch("lep_dyt_qOverPt", t.lep_dyt_qOverPt, "lep_dyt_qOverPt[2]/F");
    tree->Branch("lep_Mu50_triggerMatchPt", t.lep_Mu50_triggerMatchPt, "lep_Mu50_triggerMatchPt[2]/F");
    tree->Branch("lep_Mu50_triggerMatchEta", t.lep_Mu50_triggerMatchEta, "lep_Mu50_triggerMatchEta[2]/F");
    tree->Branch("lep_Mu50_triggerMatchPhi", t.lep_Mu50_triggerMatchPhi, "lep_Mu50_triggerMatchPhi[2]/F");
    tree->Branch("lep_OldMu100_triggerMatchPt", t.lep_OldMu100_triggerMatchPt, "lep_OldMu100_triggerMatchPt[2]/F");
    tree->Branch("lep_OldMu100_triggerMatchEta", t.lep_OldMu100_triggerMatchEta, "lep_OldMu100_triggerMatchEta[2]/F");
    tree->Branch("lep_OldMu100_triggerMatchPhi", t.lep_OldMu100_triggerMatchPhi, "lep_OldMu100_triggerMatchPhi[2]/F");
    tree->Branch("lep_TkMu100_triggerMatchPt", t.lep_TkMu100_triggerMatchPt, "lep_TkMu100_triggerMatchPt[2]/F");
    tree->Branch("lep_TkMu100_triggerMatchEta", t.lep_TkMu100_triggerMatchEta, "lep_TkMu100_triggerMatchEta[2]/F");
    tree->Branch("lep_TkMu100_triggerMatchPhi", t.lep_TkMu100_triggerMatchPhi, "lep_TkMu100_triggerMatchPhi[2]/F");
    tree->Branch("lep_chi2dof", t.lep_chi2dof, "lep_chi2dof[2]/F");
    tree->Branch("lep_dB", t.lep_dB, "lep_dB[2]/F");
    tree->Branch("lep_sumPt", t.lep_sumPt, "lep_sumPt[2]/F");
    tree->Branch("lep_emEt", t.lep_emEt, "lep_emEt[2]/F");
    tree->Branch("lep_hadEt", t.lep_hadEt, "lep_hadEt[2]/F");
    tree->Branch("lep_hoEt", t.lep_hoEt, "lep_hoEt[2]/F");
    tree->Branch("lep_pfIso", t.lep_pfIso, "lep_pfIso[2]/F");
    tree->Branch("lep_pfIsoDB", t.lep_pfIsoDB, "lep_pfIsoDB[2]/F");
    tree->Branch("lep_timeNdof", t.lep_timeNdof, "lep_timeNdof[2]/I");
    tree->Branch("lep_timeInOut", t.lep_timeInOut, "lep_timeInOut[2]/F");
    tree->Branch("lep_timeOutIn", t.lep_timeOutIn, "lep_timeOutIn[2]/F");
    tree->Branch("lep_timeInOutErr", t.lep_timeInOutErr, "lep_timeInOutErr[2]/F");
    tree->Branch("lep_timeOutInErr", t.lep_timeOutInErr, "lep_timeOutInErr[2]/F");
    tree->Branch("lep_min_muon_dR", t.lep_min_muon_dR, "lep_min_muon_dR[2]/F");
    tree->Branch("lep_tk_numberOfValidTrackerHits", t.lep_tk_numberOfValidTrackerHits, "lep_tk_numberOfValidTrackerHits[2]/S");
    tree->Branch("lep_tk_numberOfValidTrackerLayers", t.lep_tk_numberOfValidTrackerLayers, "lep_tk_numberOfValidTrackerLayers[2]/S");
    tree->Branch("lep_tk_numberOfValidPixelHits", t.lep_tk_numberOfValidPixelHits, "lep_tk_numberOfValidPixelHits[2]/S");
    tree->Branch("lep_glb_numberOfValidTrackerHits", t.lep_glb_numberOfValidTrackerHits, "lep_glb_numberOfValidTrackerHits[2]/S");
    tree->Branch("lep_glb_numberOfValidTrackerLayers", t.lep_glb_numberOfValidTrackerLayers, "lep_glb_numberOfValidTrackerLayers[2]/S");
    tree->Branch("lep_glb_numberOfValidPixelHits", t.lep_glb_numberOfValidPixelHits, "lep_glb_numberOfValidPixelHits[2]/S");
    tree->Branch("lep_glb_numberOfValidMuonHits", t.lep_glb_numberOfValidMuonHits, "lep_glb_numberOfValidMuonHits[2]/S");
    tree->Branch("lep_glb_numberOfValidMuonDTHits", t.lep_glb_numberOfValidMuonDTHits, "lep_glb_numberOfValidMuonDTHits[2]/S");
    tree->Branch("lep_glb_numberOfValidMuonCSCHits", t.lep_glb_numberOfValidMuonCSCHits, "lep_glb_numberOfValidMuonCSCHits[2]/S");
    tree->Branch("lep_glb_numberOfValidMuonRPCHits", t.lep_glb_numberOfValidMuonRPCHits, "lep_glb_numberOfValidMuonRPCHits[2]/S");
    tree->Branch("lep_glb_muonStationsWithValidHits", t.lep_glb_muonStationsWithValidHits, "lep_glb_muonStationsWithValidHits[2]/S");
    tree->Branch("lep_glb_dtStationsWithValidHits", t.lep_glb_dtStationsWithValidHits, "lep_glb_dtStationsWithValidHits[2]/S");
    tree->Branch("lep_glb_cscStationsWithValidHits", t.lep_glb_cscStationsWithValidHits, "lep_glb_cscStationsWithValidHits[2]/S");
    tree->Branch("lep_glb_rpcStationsWithValidHits", t.lep_glb_rpcStationsWithValidHits, "lep_glb_rpcStationsWithValidHits[2]/S");
    tree->Branch("lep_glb_innermostMuonStationWithValidHits", t.lep_glb_innermostMuonStationWithValidHits, "lep_glb_innermostMuonStationWithValidHits[2]/S");
    tree->Branch("lep_glb_outermostMuonStationWithValidHits", t.lep_glb_outermostMuonStationWithValidHits, "lep_glb_outermostMuonStationWithValidHits[2]/S");
    tree->Branch("lep_tuneP_numberOfValidMuonHits", t.lep_tuneP_numberOfValidMuonHits, "lep_tuneP_numberOfValidMuonHits[2]/S");
    tree->Branch("lep_tuneP_numberOfValidMuonDTHits", t.lep_tuneP_numberOfValidMuonDTHits, "lep_tuneP_numberOfValidMuonDTHits[2]/S");
    tree->Branch("lep_tuneP_numberOfValidMuonCSCHits", t.lep_tuneP_numberOfValidMuonCSCHits, "lep_tuneP_numberOfValidMuonCSCHits[2]/S");
    tree->Branch("lep_tuneP_numberOfValidMuonRPCHits", t.lep_tuneP_numberOfValidMuonRPCHits, "lep_tuneP_numberOfValidMuonRPCHits[2]/S");
    tree->Branch("lep_tuneP_muonStationsWithValidHits", t.lep_tuneP_muonStationsWithValidHits, "lep_tuneP_muonStationsWithValidHits[2]/S");
    tree->Branch("lep_tuneP_dtStationsWithValidHits", t.lep_tuneP_dtStationsWithValidHits, "lep_tuneP_dtStationsWithValidHits[2]/S");
    tree->Branch("lep_tuneP_cscStationsWithValidHits", t.lep_tuneP_cscStationsWithValidHits, "lep_tuneP_cscStationsWithValidHits[2]/S");
    tree->Branch("lep_tuneP_rpcStationsWithValidHits", t.lep_tuneP_rpcStationsWithValidHits, "lep_tuneP_rpcStationsWithValidHits[2]/S");
    tree->Branch("lep_tuneP_innermostMuonStationWithValidHits", t.lep_tuneP_innermostMuonStationWithValidHits, "lep_tuneP_innermostMuonStationWithValidHits[2]/S");
    tree->Branch("lep_tuneP_outermostMuonStationWithValidHits", t.lep_tuneP_outermostMuonStationWithValidHits, "lep_tuneP_outermostMuonStationWithValidHits[2]/S");
    tree->Branch("lep_numberOfMatches", t.lep_numberOfMatches, "lep_numberOfMatches[2]/S");
    tree->Branch("lep_numberOfMatchedStations", t.lep_numberOfMatchedStations, "lep_numberOfMatchedStations[2]/S");
    tree->Branch("lep_numberOfMatchedRPCLayers",t.lep_numberOfMatchedRPCLayers, "lep_numberOfMatchedRPCLayers[2]/S");
    tree->Branch("lep_stationMask", t.lep_stationMask, "lep_stationMask[2]/I");
    tree->Branch("lep_numberOfChambers", t.lep_numberOfChambers, "lep_numberOfChambers[2]/I");
    tree->Branch("lep_numberOfChambersNoRPC", t.lep_numberOfChambersNoRPC, "lep_numberOfChambersNoRPC[2]/I");
    tree->Branch("lep_stationGapMaskDistance", t.lep_stationGapMaskDistance, "lep_stationGapMaskDistance[2]/I");
    tree->Branch("lep_stationGapMaskPull", t.lep_stationGapMaskPull, "lep_stationGapMaskPull[2]/I");
    tree->Branch("lep_isGlobalMuon", t.lep_isGlobalMuon, "lep_isGlobalMuon[2]/O");
    tree->Branch("lep_isTrackerMuon", t.lep_isTrackerMuon, "lep_isTrackerMuon[2]/O");
    tree->Branch("GoodDataRan", &t.GoodDataRan, "GoodDataRan/O");
    tree->Branch("GoodVtx", &t.GoodVtx, "GoodVtx/O");
    tree->Branch("METFilter", &t.METFilter, "METFilter/O");
    tree->Branch("met_pt", &t.met_pt, "met_pt/F");
    tree->Branch("met_phi", &t.met_phi, "met_phi/F");
    tree->Branch("nJets", &t.nJets, "nJets/I");
    tree->Branch("jet_pt", &t.jet_pt);
    tree->Branch("jet_pt_Uncorrected", &t.jet_pt_Uncorrected);
    tree->Branch("jet_pt_L1FastJet", &t.jet_pt_L1FastJet);
    tree->Branch("jet_pt_L2Relative", &t.jet_pt_L2Relative);
    tree->Branch("jet_pt_L3Absolute", &t.jet_pt_L3Absolute);
    tree->Branch("jet_pt_L2L3Residual", &t.jet_pt_L2L3Residual);
    tree->Branch("jet_eta", &t.jet_eta);
    tree->Branch("jet_phi", &t.jet_phi);
    tree->Branch("jet_partonFlavour", &t.jet_partonFlavour);
    tree->Branch("jet_hadronFlavour", &t.jet_hadronFlavour);
    tree->Branch("jet_NHF", &t.jet_NHF);
    tree->Branch("jet_NEMF", &t.jet_NEMF);
    tree->Branch("jet_CHF", &t.jet_CHF);
    tree->Branch("jet_MUF", &t.jet_MUF);
    tree->Branch("jet_CEMF", &t.jet_CEMF);
    tree->Branch("jet_NumConst", &t.jet_NumConst);
    tree->Branch("jet_NumNeutralParticles", &t.jet_NumNeutralParticles);
    tree->Branch("jet_CHM", &t.jet_CHM);
    tree->Branch("jet_pfDeepCSVJetTags_probb", &t.jet_pfDeepCSVJetTags_probb);
    tree->Branch("jet_pfDeepCSVJetTags_probc", &t.jet_pfDeepCSVJetTags_probc);
    tree->Branch("jet_pfDeepCSVJetTags_probudsg", &t.jet_pfDeepCSVJetTags_probudsg);
    tree->Branch("jet_pfDeepCSVJetTags_probbb", &t.jet_pfDeepCSVJetTags_probbb);
    tree->Branch("jet_pfDeepCSVJetTags_probcc", &t.jet_pfDeepCSVJetTags_probcc);

    if (fill_gen_info) {
        tree->Branch("genWeight", &t.genWeight, "genWeight/F");
        tree->Branch("gen_res_mass", &t.gen_res_mass, "gen_res_mass/F");
        tree->Branch("gen_res_pt", &t.gen_res_pt, "gen_res_pt/F");
        tree->Branch("gen_res_rap", &t.gen_res_rap, "gen_res_rap/F");
        tree->Branch("gen_res_eta", &t.gen_res_eta, "gen_res_eta/F");
        tree->Branch("gen_res_phi", &t.gen_res_phi, "gen_res_phi/F");
        tree->Branch("gen_dil_mass", &t.gen_dil_mass, "gen_dil_mass/F");
        tree->Branch("gen_dil_pt", &t.gen_dil_pt, "gen_dil_pt/F");
        tree->Branch("gen_dil_rap", &t.gen_dil_rap, "gen_dil_rap/F");
        tree->Branch("gen_dil_eta", &t.gen_dil_eta, "gen_dil_eta/F");
        tree->Branch("gen_dil_phi", &t.gen_dil_phi, "gen_dil_phi/F");
        tree->Branch("gen_dil_dR", &t.gen_dil_dR, "gen_dil_dR/F");
        tree->Branch("gen_dil_dPhi", &t.gen_dil_dPhi, "gen_dil_dPhi/F");
        tree->Branch("gen_dil_noib_mass", &t.gen_dil_noib_mass, "gen_dil_noib_mass/F");
        tree->Branch("gen_dil_noib_pt", &t.gen_dil_noib_pt, "gen_dil_noib_pt/F");
        tree->Branch("gen_dil_noib_rap", &t.gen_dil_noib_rap, "gen_dil_noib_rap/F");
        tree->Branch("gen_dil_noib_eta", &t.gen_dil_noib_eta, "gen_dil_noib_eta/F");
        tree->Branch("gen_dil_noib_phi", &t.gen_dil_noib_phi, "gen_dil_noib_phi/F");
        tree->Branch("gen_dil_noib_dR", &t.gen_dil_noib_dR, "gen_dil_noib_dR/F");
        tree->Branch("gen_dil_noib_dPhi", &t.gen_dil_noib_dPhi, "gen_dil_noib_dPhi/F");
        tree->Branch("gen_lep_p", t.gen_lep_p, "gen_lep_p[2]/F");
        tree->Branch("gen_lep_pt", t.gen_lep_pt, "gen_lep_pt[2]/F");
        tree->Branch("gen_lep_px", t.gen_lep_px, "gen_lep_px[2]/F");
        tree->Branch("gen_lep_py", t.gen_lep_py, "gen_lep_py[2]/F");
        tree->Branch("gen_lep_pz", t.gen_lep_pz, "gen_lep_pz[2]/F");
        tree->Branch("gen_lep_E", t.gen_lep_E, "gen_lep_E[2]/F");
        tree->Branch("gen_lep_eta", t.gen_lep_eta, "gen_lep_eta[2]/F");
        tree->Branch("gen_lep_phi", t.gen_lep_phi, "gen_lep_phi[2]/F");
        tree->Branch("gen_lep_qOverPt", t.gen_lep_qOverPt, "gen_lep_qOverPt[2]/F");
        tree->Branch("gen_lep_noib_pt", t.gen_lep_noib_pt, "gen_lep_noib_pt[2]/F");
        tree->Branch("gen_lep_noib_px", t.gen_lep_noib_px, "gen_lep_noib_px[2]/F");
        tree->Branch("gen_lep_noib_py", t.gen_lep_noib_py, "gen_lep_noib_py[2]/F");
        tree->Branch("gen_lep_noib_pz", t.gen_lep_noib_pz, "gen_lep_noib_pz[2]/F");
        tree->Branch("gen_lep_noib_E", t.gen_lep_noib_E, "gen_lep_noib_E[2]/F");
        tree->Branch("gen_lep_noib_eta", t.gen_lep_noib_eta, "gen_lep_noib_eta[2]/F");
        tree->Branch("gen_lep_noib_phi", t.gen_lep_noib_phi, "gen_lep_noib_phi[2]/F");
        tree->Branch("gen_lep_noib_qOverPt", t.gen_lep_noib_qOverPt, "gen_lep_noib_qOverPt[2]/F");
    }

  tree->SetAlias("OppSign",  "lep_id[0]*lep_id[1] < 0");
}

template <typename T>
float userFloat(const T& patobj, const char* name, float def=-999.) {
  return patobj.hasUserFloat(name) ? patobj.userFloat(name) : def;
}

template <typename T>
int userInt(const T& patobj, const char* name, int def=-999) {
  return patobj.hasUserInt(name) ? patobj.userInt(name) : def;
}

void NtuplerLLJets::analyze(const edm::Event& event, const edm::EventSetup&) {
    memset(&t, 0, sizeof(tree_t));

    //
    // Event Information
    //
    t.run = event.id().run();
    t.lumi = event.luminosityBlock();
    t.event = event.id().event();
 
    // Get Trigger information
    edm::Handle<edm::TriggerResults> respat;
    event.getByLabel(TriggerResults_src, respat);
    const edm::TriggerNames& namespat = event.triggerNames(*respat);
    if (namespat.triggerIndex("Flag_goodVertices") < respat->size()) {
        t.GoodDataRan = 1;
        t.GoodVtx = respat->accept(namespat.triggerIndex("Flag_goodVertices"));
        bool metFilterAccept = true;
        for ( std::vector<edm::InputTag>::iterator filterTag_i = filterTags.begin(); filterTag_i != filterTags.end(); ++filterTag_i ) {
            std::string filterTag = (*filterTag_i).label();
            metFilterAccept  *= respat->accept(namespat.triggerIndex(filterTag));
        }
        t.METFilter = metFilterAccept;
    }

    // Get Beamspot information
    edm::Handle<reco::BeamSpot> bs;
    event.getByLabel(beamspot_src, bs);
    t.beamspot_x     = bs->x0();
    t.beamspot_x_err = bs->x0Error();
    t.beamspot_y     = bs->y0();
    t.beamspot_y_err = bs->y0Error();
    t.beamspot_z     = bs->z0();
    t.beamspot_z_err = bs->z0Error();

    // Get Vertex information
    edm::Handle<reco::VertexCollection> pvs;
    event.getByLabel(vertices_src, pvs);
    t.nvertices = 0;
    BOOST_FOREACH(const reco::Vertex& vtx, *pvs) {
        if (vtx.ndof() > 4 && fabs(vtx.z()) <= 24 && fabs(vtx.position().rho()) <= 2)
            t.nvertices += 1;
    }

    if (fill_gen_info) {
        hardInteraction->Fill(event);
        double EventWeight = 1.;
        edm::Handle<GenEventInfoProduct> gen_ev_info;

        event.getByLabel(genEventInfo_, gen_ev_info);
        if (gen_ev_info.isValid()) {
            EventWeight = gen_ev_info->weight();
            t.genWeight = ( EventWeight > 0.0 ) ? 1.0 : -1.0;
        }
        else {
            EventWeight = 1.0;
            t.genWeight = 1.0;
        }
        //
        // Store Generator Level information
        //
        if (hardInteraction->IsValidForRes()) {
            t.gen_res_mass = hardInteraction->resonance->mass();
            t.gen_res_pt   = hardInteraction->resonance->pt();
            t.gen_res_rap  = hardInteraction->resonance->rapidity();
            t.gen_res_eta  = hardInteraction->resonance->eta();
            t.gen_res_phi  = hardInteraction->resonance->phi();

            t.gen_dil_noib_mass = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).mass();
            t.gen_dil_noib_pt   = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).pt();
            t.gen_dil_noib_rap  = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).Rapidity();
            t.gen_dil_noib_eta  = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).eta();
            t.gen_dil_noib_phi  = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).phi();
            t.gen_dil_noib_dR   = deltaR(*hardInteraction->lepMinusNoIB, *hardInteraction->lepPlusNoIB);
            t.gen_dil_noib_dPhi = deltaPhi(*hardInteraction->lepMinusNoIB, *hardInteraction->lepPlusNoIB);

            t.gen_lep_noib_pt[0]  = hardInteraction->lepMinusNoIB->pt();
            t.gen_lep_noib_px[0]  = hardInteraction->lepMinusNoIB->px();
            t.gen_lep_noib_py[0]  = hardInteraction->lepMinusNoIB->py();
            t.gen_lep_noib_pz[0]  = hardInteraction->lepMinusNoIB->pz();
            t.gen_lep_noib_E[0]  = hardInteraction->lepMinusNoIB->energy();
            t.gen_lep_noib_eta[0] = hardInteraction->lepMinusNoIB->eta();
            t.gen_lep_noib_phi[0] = hardInteraction->lepMinusNoIB->phi();

            t.gen_lep_noib_pt[1]  = hardInteraction->lepPlusNoIB->pt();
            t.gen_lep_noib_px[1]  = hardInteraction->lepPlusNoIB->px();
            t.gen_lep_noib_py[1]  = hardInteraction->lepPlusNoIB->py();
            t.gen_lep_noib_pz[1]  = hardInteraction->lepPlusNoIB->pz();
            t.gen_lep_noib_E[1]  = hardInteraction->lepPlusNoIB->energy();
            t.gen_lep_noib_eta[1] = hardInteraction->lepPlusNoIB->eta();
            t.gen_lep_noib_phi[1] = hardInteraction->lepPlusNoIB->phi();
        } // end if hardInteraction->IsValidForRes()
        else {
            t.gen_res_mass = -999;
            t.gen_res_pt   = -999;
            t.gen_res_rap  = -999;
            t.gen_res_eta  = -999;
            t.gen_res_phi  = -999;

            t.gen_dil_noib_mass = -999;
            t.gen_dil_noib_pt   = -999;
            t.gen_dil_noib_rap  = -999;
            t.gen_dil_noib_eta  = -999;
            t.gen_dil_noib_phi  = -999;
            t.gen_dil_noib_dR   = -999;
            t.gen_dil_noib_dPhi = -999;

            t.gen_lep_noib_pt[0]  = -999;
            t.gen_lep_noib_px[0]  = -999;
            t.gen_lep_noib_py[0]  = -999;
            t.gen_lep_noib_pz[0]  = -999;
            t.gen_lep_noib_E[0]  = -999;
            t.gen_lep_noib_eta[0] = -999;
            t.gen_lep_noib_phi[0] = -999;

            t.gen_lep_noib_pt[1]  = -999;
            t.gen_lep_noib_px[1]  = -999;
            t.gen_lep_noib_py[1]  = -999;
            t.gen_lep_noib_pz[1]  = -999;
            t.gen_lep_noib_E[1]  = -999;
            t.gen_lep_noib_eta[1] = -999;
            t.gen_lep_noib_phi[1] = -999;
        }

        if (hardInteraction->IsValid()) {
            t.gen_dil_mass = (hardInteraction->lepPlus->p4() + hardInteraction->lepMinus->p4()).mass();
            t.gen_dil_pt   = (hardInteraction->lepPlus->p4() + hardInteraction->lepMinus->p4()).pt();
            t.gen_dil_rap  = (hardInteraction->lepPlus->p4() + hardInteraction->lepMinus->p4()).Rapidity();
            t.gen_dil_eta  = (hardInteraction->lepPlus->p4() + hardInteraction->lepMinus->p4()).eta();
            t.gen_dil_phi  = (hardInteraction->lepPlus->p4() + hardInteraction->lepMinus->p4()).phi();
            t.gen_dil_dR   = deltaR(*hardInteraction->lepMinus, *hardInteraction->lepPlus);
            t.gen_dil_dPhi = deltaPhi(*hardInteraction->lepMinus, *hardInteraction->lepPlus);

            t.gen_lep_p[0]  = hardInteraction->lepMinus->p();
            t.gen_lep_pt[0]  = hardInteraction->lepMinus->pt();
            t.gen_lep_px[0]  = hardInteraction->lepMinus->px();
            t.gen_lep_py[0]  = hardInteraction->lepMinus->py();
            t.gen_lep_pz[0]  = hardInteraction->lepMinus->pz();
            t.gen_lep_E[0]  = hardInteraction->lepMinus->energy();
            t.gen_lep_eta[0] = hardInteraction->lepMinus->eta();
            t.gen_lep_phi[0] = hardInteraction->lepMinus->phi();
            t.gen_lep_qOverPt[0] = hardInteraction->lepMinus->charge() / hardInteraction->lepMinus->pt();

            t.gen_lep_p[1]  = hardInteraction->lepPlus->p();
            t.gen_lep_pt[1]  = hardInteraction->lepPlus->pt();
            t.gen_lep_px[1]  = hardInteraction->lepPlus->px();
            t.gen_lep_py[1]  = hardInteraction->lepPlus->py();
            t.gen_lep_pz[1]  = hardInteraction->lepPlus->pz();
            t.gen_lep_E[1]  = hardInteraction->lepPlus->energy();
            t.gen_lep_eta[1] = hardInteraction->lepPlus->eta();
            t.gen_lep_phi[1] = hardInteraction->lepPlus->phi();
            t.gen_lep_qOverPt[1] = hardInteraction->lepPlus->charge() / hardInteraction->lepPlus->pt();
        } // end if hardInteraction->IsValid()
        else {
            t.gen_dil_mass = -999;
            t.gen_dil_pt   = -999;
            t.gen_dil_rap  = -999;
            t.gen_dil_eta  = -999;
            t.gen_dil_phi  = -999;
            t.gen_dil_dR   = -999;
            t.gen_dil_dPhi = -999;

            t.gen_lep_p[0]  = -999;
            t.gen_lep_pt[0]  = -999;
            t.gen_lep_px[0]  = -999;
            t.gen_lep_py[0]  = -999;
            t.gen_lep_pz[0]  = -999;
            t.gen_lep_E[0]  = -999;
            t.gen_lep_eta[0] = -999;
            t.gen_lep_phi[0] = -999;
            t.gen_lep_qOverPt[0] = -999;

            t.gen_lep_p[1]  = -999;
            t.gen_lep_pt[1]  = -999;
            t.gen_lep_px[1]  = -999;
            t.gen_lep_py[1]  = -999;
            t.gen_lep_pz[1]  = -999;
            t.gen_lep_E[1]  = -999;
            t.gen_lep_eta[1] = -999;
            t.gen_lep_phi[1] = -999;
            t.gen_lep_qOverPt[1] = -999;
        }
    } // end if fill_gen_info

    //
    // Get dilepton collection
    //
    edm::Handle<pat::CompositeCandidateCollection> dils;
    event.getByLabel(dimu_src, dils);

    //
    // Loop over dil candidates in dils
    //
    BOOST_FOREACH(const pat::CompositeCandidate& dil, *dils) {
        t.dil_mass = dil.mass();
        t.dil_pt = dil.pt();
        t.dil_rap = dil.rapidity();
        t.dil_eta = dil.eta();
        t.dil_phi = dil.phi();
        t.dil_dR = deltaR(*dil.daughter(0), *dil.daughter(1));
        t.dil_dPhi = deltaPhi(*dil.daughter(0), *dil.daughter(1));
        t.dil_lep_pt[0] = dil.daughter(0)->pt();
        t.dil_lep_pt[1] = dil.daughter(1)->pt();

        // Only deal with dileptons composed of e,mu for now.
        assert(dil.numberOfDaughters() == 2);
        assert(abs(dil.daughter(0)->pdgId()) == 11 || abs(dil.daughter(0)->pdgId()) == 13);
        assert(abs(dil.daughter(1)->pdgId()) == 11 || abs(dil.daughter(1)->pdgId()) == 13);

        // set opp_sign and diff_flavor
        const bool opp_sign = dil.daughter(0)->charge() + dil.daughter(1)->charge() == 0;
        const bool diff_flavor = abs(dil.daughter(0)->pdgId()) != abs(dil.daughter(1)->pdgId());
        //const bool dimuon = abs(dil.daughter(0)->pdgId()) == 13 && abs(dil.daughter(1)->pdgId()) == 13;
        
        //
        // Loop over dil.daughters
        //
        for (size_t i = 0; i < 2; ++i) {

            // For e-mu dileptons, put the muon first. For opposite-sign
            // dileptons, always put the negative lepton first. Otherwise
            // don't mess with the order.
            size_t w = i;
            if (diff_flavor) w = abs(dil.daughter(i)->pdgId()) == 13 ? 0 : 1;
            else if (opp_sign) w = dil.daughter(i)->charge() < 0 ? 0 : 1;

            // Set lepton information
            t.lep_id[w] = dil.daughter(i)->pdgId();
            if (do_electrons){
                const reco::CandidateBaseRef& lep = dileptonDaughter(dil, i);
                const pat::Electron* ele = toConcretePtr<pat::Electron>(lep);
                t.lep_eta[w] = ele->superCluster()->eta();
            }
            else t.lep_eta[w] = dil.daughter(i)->eta();
            t.lep_phi[w] = dil.daughter(i)->phi();

            //
            // Non-muon information
            //
            if (abs(t.lep_id[w]) != 13) {
                t.lep_pt_err[w] = -999;
                t.lep_tk_p[w] = -999;
                t.lep_tk_pt[w] = -999;
                t.lep_tk_pt_err[w] = -999;
                t.lep_tk_px[w] = -999;
                t.lep_tk_py[w] = -999;
                t.lep_tk_pz[w] = -999;
                t.lep_tk_eta[w] = -999;
                t.lep_tk_phi[w] = -999;
                t.lep_tk_vz[w] = -999;
                t.lep_tk_dz[w] = -999;
                t.lep_tk_chi2[w] = -999;
                t.lep_tk_ndf[w] = -999;
                t.lep_tk_qOverPt[w] = -999;
                t.lep_glb_p[w] = -999;
                t.lep_glb_pt[w] = -999;
                t.lep_glb_pt_err[w] = -999;
                t.lep_glb_px[w] = -999;
                t.lep_glb_py[w] = -999;
                t.lep_glb_pz[w] = -999;
                t.lep_glb_eta[w] = -999;
                t.lep_glb_phi[w] = -999;
                t.lep_glb_chi2[w] = -999;
                t.lep_glb_ndf[w] = -999;
                t.lep_glb_qOverPt[w] = -999;
                t.lep_tpfms_p[w] = -999;
                t.lep_tpfms_pt[w] = -999;
                t.lep_tpfms_pt_err[w] = -999;
                t.lep_tpfms_px[w] = -999;
                t.lep_tpfms_py[w] = -999;
                t.lep_tpfms_pz[w] = -999;
                t.lep_tpfms_eta[w] = -999;
                t.lep_tpfms_phi[w] = -999;
                t.lep_tpfms_chi2[w] = -999;
                t.lep_tpfms_ndf[w] = -999;
                t.lep_tpfms_qOverPt[w] = -999;
                t.lep_picky_p[w] = -999;
                t.lep_picky_pt[w] = -999;
                t.lep_picky_pt_err[w] = -999;
                t.lep_picky_px[w] = -999;
                t.lep_picky_py[w] = -999;
                t.lep_picky_pz[w] = -999;
                t.lep_picky_eta[w] = -999;
                t.lep_picky_phi[w] = -999;
                t.lep_picky_chi2[w] = -999;
                t.lep_picky_ndf[w] = -999;
                t.lep_picky_qOverPt[w] = -999;
                t.lep_tuneP_p[w] = -999;
                t.lep_tuneP_pt[w] = -999;
                t.lep_tuneP_pt_err[w] = -999;
                t.lep_tuneP_px[w] = -999;
                t.lep_tuneP_py[w] = -999;
                t.lep_tuneP_pz[w] = -999;
                t.lep_tuneP_eta[w] = -999;
                t.lep_tuneP_phi[w] = -999;
                t.lep_tuneP_chi2[w] = -999;
                t.lep_tuneP_ndf[w] = -999;
                t.lep_tuneP_qOverPt[w] = -999;
                t.lep_dyt_p[w] = -999;
                t.lep_dyt_pt[w] = -999;
                t.lep_dyt_pt_err[w] = -999;
                t.lep_dyt_px[w] = -999;
                t.lep_dyt_py[w] = -999;
                t.lep_dyt_pz[w] = -999;
                t.lep_dyt_eta[w] = -999;
                t.lep_dyt_phi[w] = -999;
                t.lep_dyt_chi2[w] = -999;
                t.lep_dyt_ndf[w] = -999;
                t.lep_dyt_qOverPt[w] = -999;
                t.lep_Mu50_triggerMatchPt[w] = -999;
                t.lep_Mu50_triggerMatchEta[w] = -999;
                t.lep_Mu50_triggerMatchPhi[w] = -999;
                t.lep_OldMu100_triggerMatchPt[w] = -999;
                t.lep_OldMu100_triggerMatchEta[w] = -999;
                t.lep_OldMu100_triggerMatchPhi[w] = -999;
                t.lep_TkMu100_triggerMatchPt[w] = -999;
                t.lep_TkMu100_triggerMatchEta[w] = -999;
                t.lep_TkMu100_triggerMatchPhi[w] = -999;
                t.lep_chi2dof[w] = -999;
                t.lep_dB[w] = -999;
                t.lep_sumPt[w] = -999;
                t.lep_emEt[w] = -999;
                t.lep_hadEt[w] = -999;
                t.lep_hoEt[w] = -999;
                t.lep_pfIso[w] = -999;
                t.lep_pfIsoDB[w] = -999;
                t.lep_timeNdof[w] = -999;
                t.lep_timeInOut[w] = -999;
                t.lep_timeOutIn[w] = -999;
                t.lep_timeInOutErr[w] = -999;
                t.lep_timeOutInErr[w] = -999;
                t.lep_tk_numberOfValidTrackerHits[w] = -999; 
                t.lep_tk_numberOfValidTrackerLayers[w] = -999; 
                t.lep_tk_numberOfValidPixelHits[w] = -999;
                t.lep_glb_numberOfValidTrackerHits[w] = -999; 
                t.lep_glb_numberOfValidTrackerLayers[w] = -999; 
                t.lep_glb_numberOfValidPixelHits[w] = -999;
                t.lep_glb_numberOfValidMuonHits[w] = -999;
                t.lep_glb_numberOfValidMuonDTHits[w] = -999;
                t.lep_glb_numberOfValidMuonCSCHits[w] = -999;
                t.lep_glb_numberOfValidMuonRPCHits[w] = -999;
                t.lep_glb_muonStationsWithValidHits[w] = -999;
                t.lep_glb_dtStationsWithValidHits[w] = -999;
                t.lep_glb_cscStationsWithValidHits[w] = -999;
                t.lep_glb_rpcStationsWithValidHits[w] = -999;
                t.lep_glb_innermostMuonStationWithValidHits[w] = -999;
                t.lep_glb_outermostMuonStationWithValidHits[w] = -999;
                t.lep_tuneP_numberOfValidMuonHits[w] = -999;
                t.lep_tuneP_numberOfValidMuonDTHits[w] = -999;
                t.lep_tuneP_numberOfValidMuonCSCHits[w] = -999;
                t.lep_tuneP_numberOfValidMuonRPCHits[w] = -999;
                t.lep_tuneP_muonStationsWithValidHits[w] = -999;
                t.lep_tuneP_dtStationsWithValidHits[w] = -999;
                t.lep_tuneP_cscStationsWithValidHits[w] = -999;
                t.lep_tuneP_rpcStationsWithValidHits[w] = -999;
                t.lep_tuneP_innermostMuonStationWithValidHits[w] = -999;
                t.lep_tuneP_outermostMuonStationWithValidHits[w] = -999;
                t.lep_numberOfMatches[w] = -999;
                t.lep_numberOfMatchedStations[w] = -999;
                t.lep_numberOfMatchedRPCLayers[w] = -999;
                t.lep_stationMask[w] = 999;
                t.lep_isGlobalMuon[w] = false;
                t.lep_isTrackerMuon[w] = false;

                //
                // Electron Information
                //
                if (abs(t.lep_id[w]) == 11) {
                    const pat::Electron* el = toConcretePtr<pat::Electron>(dileptonDaughter(dil, i));
                    assert(el);
                    t.lep_p[w] = dil.daughter(i)->p();
                    t.lep_pt[w] = dil.daughter(i)->pt();
                    t.lep_et[w] = dil.daughter(i)->et();
                    t.lep_px[w] = dil.daughter(i)->px();
                    t.lep_py[w] = dil.daughter(i)->py();
                    t.lep_pz[w] = dil.daughter(i)->pz();
                    t.lep_E[w] = dil.daughter(i)->energy();
                    t.lep_gen_match[w] = userFloat(*el, "genMatch", 0);
                    t.lep_heep_id[w] = userInt(*el, "HEEPId", 999);
                    t.lep_min_muon_dR[w] = userFloat(*el, "min_muon_dR", 999);
                } // end if electron
            } // end if !muon

            //
            // Muon Information
            //
            else { // else of if (abs(t.lep_id[w]) != 13) 
                t.lep_heep_id[w] = 999;
                t.lep_gen_match[w] = 1;
                t.lep_min_muon_dR[w] = 999;

                //
                // Muon is always from dilepton object
                //
                const pat::Muon* mu = toConcretePtr<pat::Muon>(dileptonDaughter(dil, i));
                assert(mu);
                //
                // Default Muon info (tuneP)
                //
                const reco::Track* tk = patmuon::getPickedTrack(*mu).get();
                t.lep_p[w]     = tk->p();
                t.lep_pt[w]     = tk->pt();
                t.lep_px[w]     = tk->px();
                t.lep_py[w]     = tk->py();
                t.lep_pz[w]     = tk->pz();
                t.lep_qOverPt[w] = tk->charge() / tk->pt();
                t.lep_pt_err[w] = tk->ptError();

                //
                // Tracker Track Muon Information
                //
                if (mu->innerTrack().isNull()){
                    t.lep_tk_p[w] = -999;
                    t.lep_tk_pt[w] = -999;
                    t.lep_tk_pt_err[w] = -999;
                    t.lep_tk_px[w] = -999;
                    t.lep_tk_py[w] = -999;
                    t.lep_tk_pz[w] = -999;
                    t.lep_tk_eta[w] = -999;
                    t.lep_tk_phi[w] = -999;
                    t.lep_tk_vz[w] = -999;
                    t.lep_tk_dz[w] = -999;
                    t.lep_tk_chi2[w] = -999;
                    t.lep_tk_ndf[w] = -999;
                    t.lep_tk_qOverPt[w] = -999;
                }
                else {
                    t.lep_tk_p[w] = mu->innerTrack()->p();
                    t.lep_tk_pt[w] = mu->innerTrack()->pt();
                    t.lep_tk_pt_err[w] = mu->innerTrack()->ptError();
                    t.lep_tk_px[w] = mu->innerTrack()->px();
                    t.lep_tk_py[w] = mu->innerTrack()->py();
                    t.lep_tk_pz[w] = mu->innerTrack()->pz();
                    t.lep_tk_eta[w] = mu->innerTrack()->eta();
                    t.lep_tk_phi[w] = mu->innerTrack()->phi();
                    t.lep_tk_vz[w] = mu->innerTrack()->vz();
                    t.lep_tk_dz[w] = mu->innerTrack()->dz();
                    t.lep_tk_chi2[w] = mu->innerTrack()->chi2();
                    t.lep_tk_ndf[w] = mu->innerTrack()->ndof();
                    t.lep_tk_qOverPt[w] = (mu->charge())/(mu->innerTrack()->pt());
                }
                //
                // Global Muon Information
                //
                if (mu->globalTrack().isNull()){
                    t.lep_glb_p[w] = -999;
                    t.lep_glb_pt[w] = -999;
                    t.lep_glb_pt_err[w] = -999;
                    t.lep_glb_px[w] = -999;
                    t.lep_glb_py[w] = -999;
                    t.lep_glb_pz[w] = -999;
                    t.lep_glb_eta[w] = -999;
                    t.lep_glb_phi[w] = -999;
                    t.lep_glb_chi2[w] = -999;
                    t.lep_glb_ndf[w] = -999;
                    t.lep_glb_qOverPt[w] = -999;
                }
                else {
                    t.lep_glb_p[w] = mu->globalTrack()->p();
                    t.lep_glb_pt[w] = mu->globalTrack()->pt();
                    t.lep_glb_pt_err[w] = mu->globalTrack()->ptError();
                    t.lep_glb_px[w] = mu->globalTrack()->px();
                    t.lep_glb_py[w] = mu->globalTrack()->py();
                    t.lep_glb_pz[w] = mu->globalTrack()->pz();
                    t.lep_glb_eta[w] = mu->globalTrack()->eta();
                    t.lep_glb_phi[w] = mu->globalTrack()->phi();
                    t.lep_glb_chi2[w] = mu->globalTrack()->chi2();
                    t.lep_glb_ndf[w] = mu->globalTrack()->ndof();
                    t.lep_glb_qOverPt[w] = (mu->charge())/(mu->globalTrack()->pt());
                }
                //
                // Tracker Plus First Muon Station Muon Information
                //
                if (!(mu->tpfmsTrack().refCore().isAvailable())) {
                    t.lep_tpfms_p[w] = -999;
                    t.lep_tpfms_pt[w] = -999;
                    t.lep_tpfms_pt_err[w] = -999;
                    t.lep_tpfms_px[w] = -999;
                    t.lep_tpfms_py[w] = -999;
                    t.lep_tpfms_pz[w] = -999;
                    t.lep_tpfms_eta[w] = -999;
                    t.lep_tpfms_phi[w] = -999;
                    t.lep_tpfms_chi2[w] = -999;
                    t.lep_tpfms_ndf[w] = -999;
                    t.lep_tpfms_qOverPt[w] = -999;
                }
                else {
                    t.lep_tpfms_p[w] = mu->tpfmsTrack()->p();
                    t.lep_tpfms_pt[w] = mu->tpfmsTrack()->pt();
                    t.lep_tpfms_pt_err[w] = mu->tpfmsTrack()->ptError();
                    t.lep_tpfms_px[w] = mu->tpfmsTrack()->px();
                    t.lep_tpfms_py[w] = mu->tpfmsTrack()->py();
                    t.lep_tpfms_pz[w] = mu->tpfmsTrack()->pz();
                    t.lep_tpfms_eta[w] = mu->tpfmsTrack()->eta();
                    t.lep_tpfms_phi[w] = mu->tpfmsTrack()->phi();
                    t.lep_tpfms_chi2[w] = mu->tpfmsTrack()->chi2();
                    t.lep_tpfms_ndf[w] = mu->tpfmsTrack()->ndof();
                    t.lep_tpfms_qOverPt[w] = (mu->charge())/(mu->tpfmsTrack()->pt());
                }
                //
                // Picky Muon Information
                //
                if (!(mu->pickyTrack().refCore().isAvailable())) {
                    t.lep_picky_p[w] = -999;
                    t.lep_picky_pt[w] = -999;
                    t.lep_picky_pt_err[w] = -999;
                    t.lep_picky_px[w] = -999;
                    t.lep_picky_py[w] = -999;
                    t.lep_picky_pz[w] = -999;
                    t.lep_picky_eta[w] = -999;
                    t.lep_picky_phi[w] = -999;
                    t.lep_picky_chi2[w] = -999;
                    t.lep_picky_ndf[w] = -999;
                    t.lep_picky_qOverPt[w] = -999;
                }
                else {
                    t.lep_picky_p[w] = mu->pickyTrack()->p();
                    t.lep_picky_pt[w] = mu->pickyTrack()->pt();
                    t.lep_picky_pt_err[w] = mu->pickyTrack()->ptError();
                    t.lep_picky_px[w] = mu->pickyTrack()->px();
                    t.lep_picky_py[w] = mu->pickyTrack()->py();
                    t.lep_picky_pz[w] = mu->pickyTrack()->pz();
                    t.lep_picky_eta[w] = mu->pickyTrack()->eta();
                    t.lep_picky_phi[w] = mu->pickyTrack()->phi();
                    t.lep_picky_chi2[w] = mu->pickyTrack()->chi2();
                    t.lep_picky_ndf[w] = mu->pickyTrack()->ndof();
                    t.lep_picky_qOverPt[w] = (mu->charge())/(mu->pickyTrack()->pt());
                }
                //
                // dyt track info
                //
                if (!(mu->dytTrack().refCore().isAvailable())) {
                    t.lep_dyt_p[w] = -999;
                    t.lep_dyt_pt[w] = -999;
                    t.lep_dyt_pt_err[w] = -999;
                    t.lep_dyt_px[w] = -999;
                    t.lep_dyt_py[w] = -999;
                    t.lep_dyt_pz[w] = -999;
                    t.lep_dyt_eta[w] = -999;
                    t.lep_dyt_phi[w] = -999;
                    t.lep_dyt_chi2[w] = -999;
                    t.lep_dyt_ndf[w] = -999;
                    t.lep_dyt_qOverPt[w] = -999;
                }
                else {
                    t.lep_dyt_p[w] = mu->dytTrack()->p();
                    t.lep_dyt_pt[w] = mu->dytTrack()->pt();
                    t.lep_dyt_pt_err[w] = mu->dytTrack()->ptError();
                    t.lep_dyt_px[w] = mu->dytTrack()->px();
                    t.lep_dyt_py[w] = mu->dytTrack()->py();
                    t.lep_dyt_pz[w] = mu->dytTrack()->pz();
                    t.lep_dyt_eta[w] = mu->dytTrack()->eta();
                    t.lep_dyt_phi[w] = mu->dytTrack()->phi();
                    t.lep_dyt_chi2[w] = mu->dytTrack()->chi2();
                    t.lep_dyt_ndf[w] = mu->dytTrack()->ndof();
                    t.lep_dyt_qOverPt[w] = (mu->charge())/(mu->dytTrack()->pt());
                }
                //
                // tuneP track info
                //
                if (!(mu->tunePMuonBestTrack().refCore().isAvailable())) {
                    t.lep_tuneP_p[w] = -999;
                    t.lep_tuneP_pt[w] = -999;
                    t.lep_tuneP_pt_err[w] = -999;
                    t.lep_tuneP_px[w] = -999;
                    t.lep_tuneP_py[w] = -999;
                    t.lep_tuneP_pz[w] = -999;
                    t.lep_tuneP_eta[w] = -999;
                    t.lep_tuneP_phi[w] = -999;
                    t.lep_tuneP_chi2[w] = -999;
                    t.lep_tuneP_ndf[w] = -999;
                    t.lep_tuneP_qOverPt[w] = -999;
                    t.lep_tuneP_numberOfValidMuonHits[w] = -999;
                    t.lep_tuneP_numberOfValidMuonDTHits[w] = -999;
                    t.lep_tuneP_numberOfValidMuonCSCHits[w] = -999;
                    t.lep_tuneP_numberOfValidMuonRPCHits[w] = -999;
                    t.lep_tuneP_muonStationsWithValidHits[w] = -999;
                    t.lep_tuneP_dtStationsWithValidHits[w] = -999;
                    t.lep_tuneP_cscStationsWithValidHits[w] = -999;
                    t.lep_tuneP_rpcStationsWithValidHits[w] = -999;
                    t.lep_tuneP_innermostMuonStationWithValidHits[w] = -999;
                    t.lep_tuneP_outermostMuonStationWithValidHits[w] = -999;
                }
                else {
                    t.lep_tuneP_p[w] = mu->tunePMuonBestTrack()->p();
                    t.lep_tuneP_pt[w] = mu->tunePMuonBestTrack()->pt();
                    t.lep_tuneP_pt_err[w] = mu->tunePMuonBestTrack()->ptError();
                    t.lep_tuneP_px[w] = mu->tunePMuonBestTrack()->px();
                    t.lep_tuneP_py[w] = mu->tunePMuonBestTrack()->py();
                    t.lep_tuneP_pz[w] = mu->tunePMuonBestTrack()->pz();
                    t.lep_tuneP_eta[w] = mu->tunePMuonBestTrack()->eta();
                    t.lep_tuneP_phi[w] = mu->tunePMuonBestTrack()->phi();
                    t.lep_tuneP_chi2[w] = mu->tunePMuonBestTrack()->chi2();
                    t.lep_tuneP_ndf[w] = mu->tunePMuonBestTrack()->ndof();
                    t.lep_tuneP_qOverPt[w] = (mu->charge())/(mu->tunePMuonBestTrack()->pt());
                    // Valid Muon (all), DT, CSC, RPC hits
                    t.lep_tuneP_numberOfValidMuonHits[w] = mu->tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits();
                    t.lep_tuneP_numberOfValidMuonDTHits[w] = mu->tunePMuonBestTrack()->hitPattern().numberOfValidMuonDTHits();
                    t.lep_tuneP_numberOfValidMuonCSCHits[w] = mu->tunePMuonBestTrack()->hitPattern().numberOfValidMuonCSCHits();
                    t.lep_tuneP_numberOfValidMuonRPCHits[w] = mu->tunePMuonBestTrack()->hitPattern().numberOfValidMuonRPCHits();
                    // Valid Muon, DT, CSC, RPC, innermost, outermost Station Hits
                    t.lep_tuneP_muonStationsWithValidHits[w] = mu->tunePMuonBestTrack()->hitPattern().muonStationsWithValidHits();
                    t.lep_tuneP_dtStationsWithValidHits[w] = mu->tunePMuonBestTrack()->hitPattern().dtStationsWithValidHits();
                    t.lep_tuneP_cscStationsWithValidHits[w] = mu->tunePMuonBestTrack()->hitPattern().cscStationsWithValidHits();
                    t.lep_tuneP_rpcStationsWithValidHits[w] = mu->tunePMuonBestTrack()->hitPattern().rpcStationsWithValidHits();
                    t.lep_tuneP_innermostMuonStationWithValidHits[w] = mu->tunePMuonBestTrack()->hitPattern().innermostMuonStationWithValidHits();
                    t.lep_tuneP_outermostMuonStationWithValidHits[w] = mu->tunePMuonBestTrack()->hitPattern().outermostMuonStationWithValidHits();
                }
                //
                // Trigger Match Information
                //
                t.lep_Mu50_triggerMatchPt[w]  = userFloat(*mu, "Mu50_TriggerMatchPt",  -999);
                t.lep_Mu50_triggerMatchEta[w] = userFloat(*mu, "Mu50_TriggerMatchEta", -999);
                t.lep_Mu50_triggerMatchPhi[w] = userFloat(*mu, "Mu50_TriggerMatchPhi", -999);
                t.lep_OldMu100_triggerMatchPt[w]  = userFloat(*mu, "OldMu100_TriggerMatchPt",  -999);
                t.lep_OldMu100_triggerMatchEta[w] = userFloat(*mu, "OldMu100_TriggerMatchEta", -999);
                t.lep_OldMu100_triggerMatchPhi[w] = userFloat(*mu, "OldMu100_TriggerMatchPhi", -999);
                t.lep_TkMu100_triggerMatchPt[w]  = userFloat(*mu, "TkMu100_TriggerMatchPt",  -999);
                t.lep_TkMu100_triggerMatchEta[w] = userFloat(*mu, "TkMu100_TriggerMatchEta", -999);
                t.lep_TkMu100_triggerMatchPhi[w] = userFloat(*mu, "TkMu100_TriggerMatchPhi", -999);
                //
                // Misc. event quantities
                //
                t.lep_dB[w] = mu->dB();
                t.lep_sumPt[w] = mu->isolationR03().sumPt;
                t.lep_emEt[w] = mu->isolationR03().emEt;
                t.lep_hadEt[w] = mu->isolationR03().hadEt;
                t.lep_hoEt[w] = mu->isolationR03().hoEt;
                t.lep_pfIso[w] = mu->pfIsolationR04().sumChargedHadronPt + 
                                 mu->pfIsolationR04().sumNeutralHadronEt + 
                                 mu->pfIsolationR04().sumPhotonEt;
                t.lep_pfIsoDB[w] = mu->pfIsolationR04().sumChargedHadronPt + 
                                   std::max(mu->pfIsolationR04().sumNeutralHadronEt + 
                                            mu->pfIsolationR04().sumPhotonEt - 
                                            0.5*mu->pfIsolationR04().sumPUPt
                                        ,0.0);
                if (mu->isTimeValid()) {
                    t.lep_timeNdof[w] = mu->time().nDof;
                    t.lep_timeInOut[w] = mu->time().timeAtIpInOut;
                    t.lep_timeOutIn[w] = mu->time().timeAtIpOutIn;
                    t.lep_timeInOutErr[w] = mu->time().timeAtIpInOutErr;
                    t.lep_timeOutInErr[w] = mu->time().timeAtIpOutInErr;
                }
                else {
                    t.lep_timeNdof[w] = -999;
                    t.lep_timeInOut[w] = -999;
                    t.lep_timeOutIn[w] = -999;
                    t.lep_timeInOutErr[w] = -999;
                    t.lep_timeOutInErr[w] = -999;
                }

                //
                // Tracker track
                //
                t.lep_isTrackerMuon[w] = mu->isTrackerMuon();
                // Tracker Track tracker information
                if (mu->isTrackerMuon()){
                    t.lep_tk_numberOfValidTrackerHits[w] = mu->innerTrack()->hitPattern().numberOfValidTrackerHits();
                    t.lep_tk_numberOfValidTrackerLayers[w] = mu->innerTrack()->hitPattern().trackerLayersWithMeasurement();
                    t.lep_tk_numberOfValidPixelHits[w] = mu->innerTrack()->hitPattern().numberOfValidPixelHits();
                }
                //
                // Global track
                //
                t.lep_isGlobalMuon[w] = mu->isGlobalMuon();
                // Global track tracker information
                if (mu->isGlobalMuon()){
                    t.lep_glb_numberOfValidTrackerHits[w] = mu->globalTrack()->hitPattern().numberOfValidTrackerHits();
                    t.lep_chi2dof[w] = mu->globalTrack()->normalizedChi2();
                    t.lep_glb_numberOfValidTrackerLayers[w] = mu->globalTrack()->hitPattern().trackerLayersWithMeasurement();
                    t.lep_glb_numberOfValidPixelHits[w] = mu->globalTrack()->hitPattern().numberOfValidPixelHits();
                    // Valid Muon (all), DT, CSC, RPC hits
                    t.lep_glb_numberOfValidMuonHits[w] = mu->globalTrack()->hitPattern().numberOfValidMuonHits();
                    t.lep_glb_numberOfValidMuonDTHits[w] = mu->globalTrack()->hitPattern().numberOfValidMuonDTHits();
                    t.lep_glb_numberOfValidMuonCSCHits[w] = mu->globalTrack()->hitPattern().numberOfValidMuonCSCHits();
                    t.lep_glb_numberOfValidMuonRPCHits[w] = mu->globalTrack()->hitPattern().numberOfValidMuonRPCHits();
                    // Valid Muon, DT, CSC, RPC, innermost, outermost Station Hits
                    t.lep_glb_muonStationsWithValidHits[w] = mu->globalTrack()->hitPattern().muonStationsWithValidHits();
                    t.lep_glb_dtStationsWithValidHits[w] = mu->globalTrack()->hitPattern().dtStationsWithValidHits();
                    t.lep_glb_cscStationsWithValidHits[w] = mu->globalTrack()->hitPattern().cscStationsWithValidHits();
                    t.lep_glb_rpcStationsWithValidHits[w] = mu->globalTrack()->hitPattern().rpcStationsWithValidHits();
                    t.lep_glb_innermostMuonStationWithValidHits[w] = mu->globalTrack()->hitPattern().innermostMuonStationWithValidHits();
                    t.lep_glb_outermostMuonStationWithValidHits[w] = mu->globalTrack()->hitPattern().outermostMuonStationWithValidHits();
                }
                // number of chambers with matched segments
                t.lep_numberOfMatches[w] = mu->numberOfMatches();
                // number of stations with matched segments
                t.lep_numberOfMatchedStations[w] = mu->numberOfMatchedStations();
                // number of layers with matched rpc hits
                t.lep_numberOfMatchedRPCLayers[w] = mu->numberOfMatchedRPCLayers();
                // get bit map of stations with matched segments
                // bits 0-1-2-3 = DT stations 1-2-3-4
                // bits 4-5-6-7 = CSC stations 1-2-3-4
                t.lep_stationMask[w] = mu->stationMask();
                // number of chambers
                t.lep_numberOfChambers[w] = mu->numberOfChambers();
                // number of chambers not including RPC matches
                t.lep_numberOfChambersNoRPC[w] = mu->numberOfChambersCSCorDT();
                // distanceCut = 10cm by default (distance in cm)
                t.lep_stationGapMaskDistance[w] = mu->stationGapMaskDistance();
                // sigmaCut = 3 by default (in # sigmas)
                t.lep_stationGapMaskPull[w] = mu->stationGapMaskPull();
            } // end else of if (abs(t.lep_id[w]) != 13) 

        } // end for (int i = 0; i < 2; i++) // Loop over dilepton leptons

        //
        // more event quantites
        //
        t.cos_angle    = userFloat(dil, "cos_angle", 999);
        t.vertex_chi2  = userFloat(dil, "vertex_chi2");
        t.vertex_m     = userFloat(dil, "vertexM");
        t.vertex_m_err = userFloat(dil, "vertexMError");
        t.vertex_x     = userFloat(dil, "vertexX");
        t.vertex_x_err = userFloat(dil, "vertexXError");
        t.vertex_y     = userFloat(dil, "vertexY");
        t.vertex_y_err = userFloat(dil, "vertexYError");
        t.vertex_z     = userFloat(dil, "vertexZ");
        t.vertex_z_err = userFloat(dil, "vertexZError");

        if (opp_sign) {
            const reco::CandidateBaseRef mum = dileptonDaughterByCharge(dil, -1);
            const reco::CandidateBaseRef mup = dileptonDaughterByCharge(dil, +1);

            t.cos_cs = calcCosThetaCSAnal(mum->pz(), mum->energy(), mup->pz(), mup->energy(), dil.pt(), dil.pz(), dil.mass());
            t.chi_dilepton = exp(std::abs(mum->p4().Rapidity()-mup->p4().Rapidity()));
            t.phi_cs = calcPhiCSAnal(mum->px(), mum->py(), mup->px(), mup->py(), dil.pt(), dil.eta(), dil.phi(), dil.mass(), true);
        } // end if opp_sign
        else {
            t.cos_cs = -999;
            t.chi_dilepton = -999;
            t.phi_cs = -999;
        } // end if !opp_sign

        //
        // Fill tree
        //
        edm::Handle< std::vector< pat::MET > > mets;
        event.getByLabel(met_src, mets);
        t.met_pt = mets->front().pt();
        t.met_phi = mets->front().phi();

        edm::Handle< std::vector< pat::Jet > > jets;
        event.getByLabel(jet_src, jets);

        int nJets = 0;
        for (std::vector<pat::Jet>::const_iterator itJet = jets->begin(); itJet != jets->end(); itJet++) {
            if( itJet->pt() < 20.0 || abs(itJet->eta()) > 2.5 )
                continue;

            t.jet_pt.push_back( itJet->pt() );
            t.jet_eta.push_back( itJet->eta() );
            t.jet_phi.push_back( itJet->phi() );
            t.jet_partonFlavour.push_back( itJet->partonFlavour() );
            t.jet_hadronFlavour.push_back( itJet->hadronFlavour() );
            t.jet_NHF.push_back( itJet->neutralHadronEnergyFraction() );
            t.jet_NEMF.push_back( itJet->neutralEmEnergyFraction() );
            t.jet_CHF.push_back( itJet->chargedHadronEnergyFraction() );
            t.jet_MUF.push_back( itJet->muonEnergyFraction() );
            t.jet_CEMF.push_back( itJet->chargedEmEnergyFraction() );
            t.jet_NumConst.push_back( itJet->chargedMultiplicity()+itJet->neutralMultiplicity() );
            t.jet_NumNeutralParticles.push_back( itJet->neutralMultiplicity() );
            t.jet_CHM.push_back( itJet->chargedMultiplicity() );
            t.jet_pfDeepCSVJetTags_probb.push_back( itJet->bDiscriminator("pfDeepCSVJetTags:probb") );
            t.jet_pfDeepCSVJetTags_probc.push_back( itJet->bDiscriminator("pfDeepCSVJetTags:probc") );
            t.jet_pfDeepCSVJetTags_probudsg.push_back( itJet->bDiscriminator("pfDeepCSVJetTags:probudsg") );
            t.jet_pfDeepCSVJetTags_probbb.push_back( itJet->bDiscriminator("pfDeepCSVJetTags:probbb") );
            t.jet_pfDeepCSVJetTags_probcc.push_back( itJet->bDiscriminator("pfDeepCSVJetTags:probcc") );

            // JEC test
            // TString str_0 = TString::Format("%5d: (%.2f, %.2f, %.2f, %d, %d)",
            //                                 nJets,
            //                                 itJet->pt(),
            //                                 itJet->eta(),
            //                                 itJet->phi(),
            //                                 itJet->partonFlavour(),
            //                                 itJet->hadronFlavour());
            // std::cout << str_0 << std::endl;
            // for(auto set : itJet->availableJECSets()) {
            //     std::cout << "\tJET set: " << set << std::endl;
            //     for(auto level : itJet->availableJECLevels(set)) {
            //         auto corrJet = itJet->correctedJet(level, "none", set);
            //         TString str_1 = TString::Format("\t\t%s: (%.2f, %.2f, %.2f, %d, %d)",
            //                                         level.c_str(),
            //                                         corrJet.pt(),
            //                                         corrJet.eta(),
            //                                         corrJet.phi(),
            //                                         corrJet.partonFlavour(),
            //                                         corrJet.hadronFlavour());
            //         std::cout << str_1 << std::endl;
            //     }
            // }

            // JEC
            auto set = itJet->currentJECSet();
            auto levels = itJet->availableJECLevels(set);
            if (std::find(levels.begin(), levels.end(), "Uncorrected") != levels.end()) {
                t.jet_pt_Uncorrected.push_back(itJet->correctedJet("Uncorrected", "none", set));
            }
            else {
                t.jet_pt_Uncorrected.push_back(-999.);
            }
            if (std::find(levels.begin(), levels.end(), "L1FastJet") != levels.end()) {
                t.jet_pt_L1FastJet.push_back(itJet->correctedJet("L1FastJet", "none", set));
            }
            else {
                t.jet_pt_L1FastJet.push_back(-999.);
            }
            if (std::find(levels.begin(), levels.end(), "L2Relative") != levels.end()) {
                t.jet_pt_L2Relative.push_back(itJet->correctedJet("L2Relative", "none", set));
            }
            else {
                t.jet_pt_L2Relative.push_back(-999.);
            }
            if (std::find(levels.begin(), levels.end(), "L3Absolute") != levels.end()) {
                t.jet_pt_L3Absolute.push_back(itJet->correctedJet("L3Absolute", "none", set));
            }
            else {
                t.jet_pt_L3Absolute.push_back(-999.);
            }
            if (std::find(levels.begin(), levels.end(), "L2L3Residual") != levels.end()) {
                t.jet_pt_L2L3Residual.push_back(itJet->correctedJet("L2L3Residual", "none", set));
            }
            else {
                t.jet_pt_L2L3Residual.push_back(-999.);
            }

            nJets++;
        }
        t.nJets = nJets;

        tree->Fill();
    }
} // end NtuplerLLJets::analyze

DEFINE_FWK_MODULE(NtuplerLLJets);
