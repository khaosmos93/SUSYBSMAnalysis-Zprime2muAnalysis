
import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, offline_pt_threshold, trigger_match_2018, prescaled_trigger_match_2018
from SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2018_cff import loose_cut, trigger_match, tight_cut, allDimuons

# Below cuts should be part of loose_cut in  python/OurSelection2018_cff.py
nm1_cuts = [
    ('Pt',      'pt > 53'),
    ('DB',      'abs(dB) < 0.2'),
    ('Iso',     'isolationR03.sumPt / innerTrack.pt < 0.10'),
    ('TkLayers','globalTrack.hitPattern.trackerLayersWithMeasurement > 5'),
    ('PxHits',  'globalTrack.hitPattern.numberOfValidPixelHits >= 1'),
    ('MuHits',  'globalTrack.hitPattern.numberOfValidMuonHits > 0'),
    ('MuHits', '( (globalTrack.hitPattern.numberOfValidMuonHits > 0) || (tunePMuonBestTrack.hitPattern.numberOfValidMuonHits > 0) )'),
    ('MuMatch', ('(( numberOfMatchedStations>1 ) || ( numberOfMatchedStations==1 && ( expectedNnumberOfMatchedStations<2 || !(stationMask==1 || stationMask==16) || numberOfMatchedRPCLayers>2)))', 'isTrackerMuon')),
]


def addNminus1Histos(process, leptons_name, name, HistosFromPAT):

  path_list = []

  allDimuonsTemp = getattr(process, 'all'+name)
  dimuonsTemp    = getattr(process, name)


  for nm1_cut_name, nm1_cut in nm1_cuts:
    if type(nm1_cut) != tuple:
      nm1_cut = (nm1_cut,)

    lc = loose_cut
    for c in nm1_cut:
      if c not in lc:
        raise ValueError('nm1_cut "%s" not in nm1_cut string "%s"' % (c, lc))
      lc = lc.replace(' && ' + c, '') # Relies on none of the cuts above being first in the list.

    obj_no = allDimuonsTemp.clone(loose_cut = lc)
    setattr(process, 'allDimuonsNo' + nm1_cut_name, obj_no)

    # obj_ti = obj_no.clone(tight_cut = tight_cut + ' && ' + ' && '.join(nm1_cut))
    # setattr(process, 'allDimuonsTi' + nm1_cut_name, obj_ti)

  process.allDimuonsNoNo       = allDimuonsTemp.clone()
  process.allDimuonsNoTrgMatch = allDimuonsTemp.clone(tight_cut = tight_cut.replace(trigger_match_2018, ''))

  alldimus = [x for x in dir(process) if 'allDimuonsNo' in x or 'allDimuonsTi' in x]

  # Sanity check that the replaces above did something.
  for x in alldimus:
    if 'NoNo' in x:
      continue
    o = getattr(process, x)
    assert o.loose_cut.value() != loose_cut or o.tight_cut.value() != tight_cut

  # For all the allDimuons producers, make dimuons producers, and
  # analyzers to make the histograms.
  for alld in alldimus:
    dimu = dimuonsTemp.clone(src = alld)
    dimu_name = alld.replace('allD', 'd')
    setattr(process, dimu_name, dimu)
    hists = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(dimu_name))
    setattr(process, dimu_name.replace('dimuons', ''), hists)
    path_list.append(getattr(process, alld) * dimu * hists)

  # Handle the cuts that have to be applied at the
  # Zprime2muCompositeCandidatePicker level.
  process.allDimuonsBASE = allDimuonsTemp.clone()
  path_list.append(process.allDimuonsBASE)
  process.dimuonsNoB2B     = dimuonsTemp.clone(src = 'allDimuonsBASE')
  process.dimuonsNoVtxProb = dimuonsTemp.clone(src = 'allDimuonsBASE')
  process.dimuonsNoDptPt   = dimuonsTemp.clone(src = 'allDimuonsBASE')
  delattr(process.dimuonsNoB2B,     'back_to_back_cos_angle_min')
  delattr(process.dimuonsNoVtxProb, 'vertex_chi2_max')
  delattr(process.dimuonsNoDptPt,   'dpt_over_pt_max')

  for dimu in ['dimuonsNoB2B', 'dimuonsNoVtxProb', 'dimuonsNoDptPt']:
    hists = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(dimu))
    setattr(process, dimu.replace('dimuons', ''), hists)
    path_list.append(getattr(process, dimu) * hists)

  # Special case to remove |dB| and B2B cuts simultaneously, as they can
  # be correlated (anti-cosmics).
  process.allDimuonsNoCosm = allDimuonsTemp.clone(loose_cut = loose_cut.replace(' && abs(dB) < 0.2', ''))
  process.dimuonsNoCosm = dimuonsTemp.clone(src = 'allDimuonsNoCosm')
  delattr(process.dimuonsNoCosm, 'back_to_back_cos_angle_min')
  process.NoCosm = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = 'dimuonsNoCosm')
  path_list.append(process.allDimuonsNoCosm * process.dimuonsNoCosm * process.NoCosm)

  return path_list

