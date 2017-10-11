#!/usr/bin/env python

miniAOD = True

# intime_bin numbering: bin 0 = 0-5, bin 1 = 6-11, bin 2 = 12-26
# late_bin numbering: bin 0 = 0-9, bin 2 = 10-26
intime_bin, late_bin = -1, -1
check_prescaled_path = False


################################################################################

import sys, os
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import cms, process

process.maxEvents.input = -1
process.source.fileNames = [#'/store/mc/RunIISpring16MiniAODv1/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/30000/02244373-7E03-E611-B581-003048F5B2B4.root',
            #'file:/u/user/msoh/Samples/MC/DY/M6000/0A222CC8-16C7-E611-BFC2-FA163E39B5B6.root', #M4500-6000
            #'file:/u/user/msoh/Samples/MC/DY/M6000/52FEE967-17C7-E611-9A69-FA163E2EED86.root',
            #'file:/u/user/msoh/Samples/MC/DY/M6000/A4999F0B-17C7-E611-A246-1866DAEA6D0C.root',
            #'file:/u/user/msoh/Samples/MC/DY/M6000/A8CCA36B-17C7-E611-8E4A-0090FAA57D64.root',
            #'file:/u/user/msoh/Samples/MC/DY/M6000/B86337F0-16C7-E611-BA6A-24BE05C44B91.root',
            #'file:/u/user/msoh/Samples/MC/DY/M6000/EA325C83-17C7-E611-AFA4-001E67E5E8B6.root',
            #'file:/u/user/msoh/Samples/MC/DY/M6000/F214122A-17C7-E611-9134-001E674FC800.root',

                            ]
process.options.wantSummary = True
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # default 1000


from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import  leptons,leptonsMini, muonPhotonMatchMiniAOD, muonPhotonMatch, allDimuons, dimuons, rec_level_module

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.HardInteractionFilter_cfi')
process.HardInteractionFilterRes = process.HardInteractionFilter.clone(use_resonance_mass=True)

process.load('SUSYBSMAnalysis.Zprime2muAnalysis.EfficiencyFromMC_cfi')

import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2016_cff as OurSelection
process.allDimuonsOur = OurSelection.allDimuons.clone()
process.dimuonsOur = OurSelection.dimuons.clone(src = 'allDimuonsOur')
process.OurEfficiencyFromMCMini = process.EfficiencyFromMCMini.clone(dimuon_src = 'dimuonsOur', acceptance_max_eta_2 = 2.4)
process.OurEfficiencyFromMC = process.EfficiencyFromMC.clone(dimuon_src = 'dimuonsOur', acceptance_max_eta_2 = 2.4)
process.OurEfficiencyFromMCnoTrigger = process.EfficiencyFromMCnoTrigger.clone(dimuon_src = 'dimuonsOur', acceptance_max_eta_2 = 2.4) ### NO TRIGGER PROCESS

# Temporarily disable explicit checks on Level-1 decision until we
# figure out which branch to use.
#for eff in [process.EfficiencyFromMC, process.OurEfficiencyFromMC]:
### NO TRIGGER PROCESSES
for eff in [process.EfficiencyFromMCnoTrigger, process.OurEfficiencyFromMCnoTrigger, process.OurEfficiencyFromMCMini, process.EfficiencyFromMCMini]:
    eff.check_l1 = False
    eff.useWeight = True # Apply SF weight
    eff.useRandom = False #True # weight with random gaus


# this will get all the Histostunep, Histospicky, Histosglobal, etc. below.
if miniAOD:
    process.HardInteractionFilterRes.hardInteraction.src = cms.InputTag('prunedGenParticles')
    process.HardInteractionFilter.hardInteraction.src = cms.InputTag('prunedGenParticles')
    from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
    electrons_miniAOD(process)


    process.leptons = process.leptonsMini.clone()
    process.Zprime2muAnalysisSequencePlain = cms.Sequence(process.HardInteractionFilterRes *process.muonPhotonMatchMiniAOD * process.egmGsfElectronIDSequence * process.leptons * process.allDimuons * process.dimuons)

else:
    process.Zprime2muAnalysisSequencePlain = cms.Sequence(process.HardInteractionFilterRes *process.muonPhotonMatch * process.leptons * process.allDimuons * process.dimuons)


p2 = process.Zprime2muAnalysisSequencePlain

if miniAOD:
    p2 = p2 * process.EfficiencyFromMCMini * process.allDimuonsOur * process.dimuonsOur * process.OurEfficiencyFromMCMini

else:
    p2 = p2 * process.EfficiencyFromMC * process.allDimuonsOur * process.dimuonsOur * process.OurEfficiencyFromMC

process.p2 = cms.Path(p2)

f = file('outfile', 'w')
f.write(process.dumpPython())
f.close()

###############################################################################


