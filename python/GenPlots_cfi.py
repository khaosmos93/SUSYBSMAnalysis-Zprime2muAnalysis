import FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction,hardInteraction_MiniAOD

GenPlots_MiniAOD = cms.EDAnalyzer('GenPlots',
                    isMC = cms.bool(False),
                    hardInteraction = hardInteraction_MiniAOD,
                    useMadgraphWeight = cms.bool(True),
                    both_in_acc = cms.bool(False),
)
