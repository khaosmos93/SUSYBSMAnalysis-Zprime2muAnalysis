#!/usr/bin/env python


miniAOD = True
Electrons = False
ex = '20190301'

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable
# Set temporary global tags here, may be changed later
MCGT = '102X_upgrade2018_realistic_v12'

import sys, os, FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import goodDataFiltersMiniAOD
from SUSYBSMAnalysis.Zprime2muAnalysis.NtupleFromPAT_cfi import NtupleFromPAT_MiniAOD,NtupleFromPAT
from SUSYBSMAnalysis.Zprime2muAnalysis.MCSamples import samples

process.source.fileNames = [
        #'/store/data/Run2018A/SingleMuon/MINIAOD/06Jun2018-v1/410000/CCA4DBD1-FF83-E811-988F-FA163E5991FE.root'
        #'/store/data/Run2018D/SingleMuon/MINIAOD/PromptReco-v2/000/322/068/00000/F8DCA3B9-41B0-E811-8B23-FA163E279E4C.root'
        '/store/mc/RunIIAutumn18MiniAOD/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/120000/078DB2B1-40DD-634D-A3CF-D2E377CAFA48.root'
           ]

process.maxEvents.input = -1
#process.options.wantSummary = cms.untracked.bool(True)# false di default
process.MessageLogger.cerr.FwkReport.reportEvery = 10000 # default 1000

from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, prescaled_trigger_match, trigger_paths, prescaled_trigger_paths, overall_prescale, offline_pt_threshold, prescaled_offline_pt_threshold, trigger_filters, trigger_path_names, prescaled_trigger_filters, prescaled_trigger_path_names, prescaled_trigger_match_2018, trigger_match_2018

# The histogramming module that will be cloned multiple times below
# for making histograms with different cut/dilepton combinations.

if miniAOD:
    from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
    electrons_miniAOD(process)

    from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT_MiniAOD as HistosFromPAT
    ####################################
    HistosFromPAT.leptonsFromDileptons = True
    HistosFromPAT.usekFactor = False #### Set TRUE to use K Factor #####
    ####################################
    ZSkim = False #### Set TRUE to skim dy50to120 with a Z pt < 100 GeV #####
    ####################################
    doNminus1 = True  #### Set True to produce N-1 histograms

else:
    from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT
    HistosFromPAT.leptonsFromDileptons = True


# Since the prescaled trigger comes with different prescales in
# different runs/lumis, this filter prescales it to a common factor to
# make things simpler.
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrescaleToCommon_cff')
process.PrescaleToCommon.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommon.overall_prescale = overall_prescale
process.PrescaleToCommonMiniAOD.trigger_paths = prescaled_trigger_paths
process.PrescaleToCommonMiniAOD.overall_prescale = overall_prescale # 500 for 2018

# These modules define the basic selection cuts. For the monitoring
# sets below, we don't need to define a whole new module, since they
# just change one or two cuts -- see below.
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelectionDec2012_cff as OurSelectionDec2012
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2016_cff as OurSelection2016
import SUSYBSMAnalysis.Zprime2muAnalysis.OurSelection2018_cff as OurSelection2018

# CandCombiner includes charge-conjugate decays with no way to turn it
# off. To get e.g. mu+mu+ separate from mu-mu-, cut on the sum of the
# pdgIds (= -26 for mu+mu+).
dils = [
    ('MuonsPlusMuonsMinus', '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-','daughter(0).pdgId() + daughter(1).pdgId() == 0'),
    ('MuonsPlusMuonsPlus',  '%(leptons_name)s:muons@+ %(leptons_name)s:muons@+','daughter(0).pdgId() + daughter(1).pdgId() == -26'),
    ('MuonsMinusMuonsMinus','%(leptons_name)s:muons@- %(leptons_name)s:muons@-','daughter(0).pdgId() + daughter(1).pdgId() == 26'),
    ('MuonsSameSign',       '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',''),
    ('MuonsAllSigns',       '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',''),
    ]

# Define sets of cuts for which to make plots. If using a selection
# that doesn't have a trigger match, need to re-add a hltHighLevel
# filter somewhere below.
cuts = {
    #'Our2012'  : OurSelectionDec2012,
    #'Our2016'  : OurSelection2016,
    'Our2018'  : OurSelection2018,
    'Our2018MuPrescaled' : OurSelection2018,
    'Simple'   : OurSelection2018, # The selection cuts in the module listed here are ignored below.
    }

if miniAOD and Electrons:
    dils = [\
        ('MuonsPlusMuonsMinus',     '%(leptons_name)s:muons@+ %(leptons_name)s:muons@-',    'daughter(0).pdgId() + daughter(1).pdgId() == 0'),
        ('MuonsPlusMuonsPlus',      '%(leptons_name)s:muons@+ %(leptons_name)s:muons@+',    'daughter(0).pdgId() + daughter(1).pdgId() == -26'),
        ('MuonsMinusMuonsMinus',    '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',    'daughter(0).pdgId() + daughter(1).pdgId() == 26'),
        ('MuonsSameSign',           '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',    ''),
        ('MuonsAllSigns',           '%(leptons_name)s:muons@- %(leptons_name)s:muons@-',    ''),
        ('MuonsPlusElectronsMinus', '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@-','daughter(0).pdgId() + daughter(1).pdgId() == -2'),
        ('MuonsMinusElectronsPlus', '%(leptons_name)s:muons@- %(leptons_name)s:electrons@+','daughter(0).pdgId() + daughter(1).pdgId() == 2'),
        ('MuonsPlusElectronsPlus',  '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+','daughter(0).pdgId() + daughter(1).pdgId() == -24'),
        ('MuonsMinusElectronsMinus','%(leptons_name)s:muons@- %(leptons_name)s:electrons@-','daughter(0).pdgId() + daughter(1).pdgId() == 24'),
        ('MuonsElectronsOppSign',   '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@-',''),
        ('MuonsElectronsSameSign',  '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',''),
        ('MuonsElectronsAllSigns',  '%(leptons_name)s:muons@+ %(leptons_name)s:electrons@+',''),
        ]
    
    cuts = {
        'Our2012'  : OurSelectionDec2012,
        'Our2016'  : OurSelection2016,
        'Our2018'  : OurSelection2018,
        'EmuVeto'  : OurSelectionDec2012, # this switches on the dRMuEl veto
        'Simple'   : OurSelectionDec2012, # The selection cuts in the module listed here are ignored below.
        }
    

# Loop over all the cut sets defined and make the lepton, allDilepton
# (combinatorics only), and dilepton (apply cuts) modules for them.
for cut_name, Selection in cuts.iteritems():
    # Keep track of modules to put in the path for this set of cuts.
    path_list = []

    # Clone the LeptonProducer to make leptons with the set of cuts
    # we're doing here flagged.  I.e., muon_cuts in LeptonProducer
    # just marks each muon with a userInt "cutFor" that is 0 if it
    # passes the cuts, and non-0 otherwise; it does not actually drop
    # any of the muons. The cutFor flag actually gets ignored by the
    # LooseTightPairSelector in use for all the cuts above, at
    # present.
    if miniAOD: path_list.append(process.egmGsfElectronIDSequence)
        
    leptons_name = cut_name + 'Leptons'
    if cut_name == 'Simple':
        muon_cuts = ''
    elif 'MuPrescaled' in cut_name:
        muon_cuts = Selection.loose_cut.replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
    else:
        muon_cuts = Selection.loose_cut

    if miniAOD:
        leptons = process.leptonsMini.clone(muon_cuts = muon_cuts)
        if (len(trigger_filters)>0 or len(prescaled_trigger_filters)>0) and ('Our2018' in cut_name or cut_name=='Simple'):
            leptons.trigger_filters = trigger_filters
            leptons.trigger_path_names = trigger_path_names
            leptons.prescaled_trigger_filters = prescaled_trigger_filters
            leptons.prescaled_trigger_path_names = prescaled_trigger_path_names
    else:
        leptons = process.leptons.clone(muon_cuts = muon_cuts)

    if  miniAOD and Electrons:
        if cut_name == 'EmuVeto':
            leptons.electron_muon_veto_dR = 0.1

    # Keep using old TuneP for past selections
    #if 'Dec2012' not in Selection.__file__:
    #    leptons.muon_track_for_momentum = cms.string('TunePNew')
    setattr(process, leptons_name, leptons)
    path_list.append(leptons)

    # Make all the combinations of dileptons we defined above.
    for dil_name, dil_decay, dil_cut in dils:
        # For the EmuVeto path, we only care about e-mu events.
        if cut_name == 'EmuVeto' and 'Electron' not in dil_name:
            continue

        # For the MuPrescaled paths, we don't care about e-mu events.
        if 'MuPrescaled' in cut_name and 'Electron' in dil_name:
            continue

        # Unique names for the modules: allname for the allDileptons,
        # and name for dileptons.
        name = cut_name + dil_name
        allname = 'all' + name

        alldil = Selection.allDimuons.clone(decay = dil_decay % locals(), cut = dil_cut)
        if 'AllSigns' in dil_name:
            alldil.checkCharge = cms.bool(False)

        dil = Selection.dimuons.clone(src = cms.InputTag(allname))

        # Implement the differences to the selections; currently, as
        # in Zprime2muCombiner, the cuts in loose_cut and
        # tight_cut are the ones actually used to drop leptons, and
        # not the ones passed into the LeptonProducer to set cutFor above.
        if cut_name == 'Simple':
            alldil.electron_cut_mask = cms.uint32(0)
            alldil.loose_cut = 'isGlobalMuon && pt > 20.'
            alldil.tight_cut = ''
            dil.max_candidates = 100
            dil.sort_by_pt = True
            dil.do_remove_overlap = False
            dil.prefer_Z=False
            if hasattr(dil, 'back_to_back_cos_angle_min'):
                delattr(dil, 'back_to_back_cos_angle_min')
            if hasattr(dil, 'vertex_chi2_max'):
                delattr(dil, 'vertex_chi2_max')
            if hasattr(dil, 'dpt_over_pt_max'):
                delattr(dil, 'dpt_over_pt_max')
        elif cut_name == 'OurNoIso':
            alldil.loose_cut = alldil.loose_cut.value().replace(' && isolationR03.sumPt / innerTrack.pt < 0.10', '')
        elif 'MuPrescaled' in cut_name:
            alldil.loose_cut = alldil.loose_cut.value().replace('pt > %s' % offline_pt_threshold, 'pt > %s' % prescaled_offline_pt_threshold)
            assert alldil.tight_cut == trigger_match_2018
            if len(prescaled_trigger_filters)>0:
                alldil.tight_cut = prescaled_trigger_match_2018
            else:
                alldil.tight_cut = prescaled_trigger_match

        # Histos now just needs to know which leptons and dileptons to use.
        histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'muons'), dilepton_src = cms.InputTag(name))

        # Add all these modules to the process and the path list.
        setattr(process, allname, alldil)
        setattr(process, name, dil)
        setattr(process, name + 'Histos', histos)
        path_list.append(alldil * dil * histos)

        # -- Add N-1 histograms -- #
        if doNminus1 and cut_name == 'Our2018' and dil_name == 'MuonsPlusMuonsMinus':
            from addNminus1 import *
            nm1_path_list = addNminus1Histos(process, leptons_name, name, HistosFromPAT)
            path_list = path_list + nm1_path_list

    # Finally, make the path for this set of cuts.
    pathname = 'path' + cut_name
    if miniAOD:
        process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DileptonPreselector_cfi')
        process.load("SUSYBSMAnalysis.Zprime2muAnalysis.EventCounter_cfi")
        pobj = process.EventCounter * process.dileptonPreselector *  process.muonPhotonMatchMiniAOD * reduce(lambda x,y: x*y, path_list)
    else:
        pobj = process.muonPhotonMatch * reduce(lambda x,y: x*y, path_list)


    if 'MuPrescaled' in cut_name:
        if miniAOD:
            pobj = process.PrescaleToCommonMiniAOD * pobj 
        else:
            pobj = process.PrescaleToCommon * pobj 

    path = cms.Path(pobj)
    setattr(process, pathname, path)


def apply_gen_filters(process,sampleName):
    from SUSYBSMAnalysis.Zprime2muAnalysis.MCFilters_cfi import DYPtZskim, TTbarGenMassFilter, DibosonGenMassFilter, TauTauFilter
    addFilter = False
    if miniAOD:
        process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
    if ('dy50to120' in sampleName or 'dyInclusive' in sampleName) and ZSkim:
        mcFilter = DYPtZskim.clone()
        addFilter = True
    elif 'ttbar_lep50to500' in sampleName:
        mcFilter = TTbarGenMassFilter.clone()
        addFilter = True
    elif 'WWinclusive' in sampleName:
        mcFilter = DibosonGenMassFilter.clone()
        addFilter = True
    elif 'dyInclusive50' in sampleName:
        mcFilter = TauTauFilter.clone()
        addFilter = True
    if addFilter:
        if not miniAOD:
            mcFilter.src = cms.InputTag('prunedMCLeptons')
        setattr(process,sampleName+'Filter',mcFilter)
        mcFilterPath = getattr(process,sampleName+'Filter')
        for path_name, path in process.paths.iteritems():
            getattr(process,path_name).insert(0,mcFilterPath)

def ntuplify(process, cut='Simple', dil_name='MuonsAllSigns', fill_gen_info=False):
    dimu_src_tag = cut+dil_name
    if miniAOD:
        ntuple = NtupleFromPAT_MiniAOD.clone(dimu_src=cms.InputTag(dimu_src_tag))
    else:
        ntuple = NtupleFromPAT.clone(dimu_src=cms.InputTag(dimu_src_tag))

    if fill_gen_info:
        from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction
        ntuple.hardInteraction = hardInteraction
        if miniAOD:
            ntuple.TriggerResults_src = cms.InputTag('TriggerResults','','PAT')
            process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
            obj = process.prunedMCLeptons
            obj.src = cms.InputTag('prunedGenParticles')

    ntuple_name = cut+dil_name+'Ntuple'
    setattr(process,ntuple_name,ntuple)
    ntuplepath = getattr(process,ntuple_name)

    if hasattr(process,'path'+cut): 
        path = getattr(process,'path'+cut)
        if fill_gen_info:
            path *= obj * ntuplepath
        else:
            path *= ntuplepath

def printify(process):
    process.MessageLogger.categories.append('PrintEvent')

    process.load('HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi')
    process.triggerSummaryAnalyzerAOD.inputTag = cms.InputTag('hltTriggerSummaryAOD','','HLT')
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.triggerSummaryAnalyzerAOD

    process.PrintOriginalMuons = cms.EDAnalyzer('PrintEvent', muon_src = cms.InputTag('cleanPatMuonsTriggerMatch'), trigger_results_src = cms.InputTag('TriggerResults','','HLT'))
    process.pathSimple *= process.PrintOriginalMuons

    pe = process.PrintEventSimple = cms.EDAnalyzer('PrintEvent', dilepton_src = cms.InputTag('SimpleMuonsPlusMuonsMinus'))
    if hasattr(process, 'pathSimple'):
        process.pathSimple *= process.PrintEventSimple

    #- 2011-2012 selection (Nlayers > 8)
    #process.PrintEventOurNew = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsPlusMuonsMinus'))
    #process.PrintEventOurNewSS = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsSameSign'))
    #process.PrintEventOurNewEmu = pe.clone(dilepton_src = cms.InputTag('OurNewMuonsElectronsOppSign'))
    #process.pathOurNew *= process.PrintEventOurNew * process.PrintEventOurNewSS * process.PrintEventOurNewEmu

    #- December 2012 selection (Nlayers > 5, re-tuned TuneP, dpT/pT < 0.3)
    if hasattr(process, 'pathOur2012'):
        process.PrintEventOur2012    = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsPlusMuonsMinus'))
        process.PrintEventOur2012SS  = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsSameSign'))
        process.PrintEventOur2012Emu = pe.clone(dilepton_src = cms.InputTag('Our2012MuonsElectronsOppSign'))
        process.pathOur2012 *= process.PrintEventOur2012 * process.PrintEventOur2012SS * process.PrintEventOur2012Emu

def check_prescale(process, trigger_paths, hlt_process_name='HLT'):
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.CheckPrescale_cfi')
    process.CheckPrescale.trigger_paths = cms.vstring(*trigger_paths)
    process.CheckPrescale.dump_prescales = cms.untracked.bool(False)
    process.pCheckPrescale = cms.Path(process.CheckPrescale)

def add_filters(process,is_mc=True):
    for cut_name, Selection in cuts.iteritems():
        path_name = 'path'+cut_name
        if hasattr(process,path_name) and cut_name != 'Simple':
            process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
            for dataFilter in goodDataFiltersMiniAOD:
                if is_mc:
                    dataFilter.src = cms.InputTag('TriggerResults','','PAT') # to submit MC
                getattr(process,path_name).insert(0,dataFilter)

def for_data(process):
    # Add filters
    add_filters(process,is_mc=False)
    # make a SimpleMuonsAllSignsNtuple
    ntuplify(process) 
    # make a Our2018MuonsPlusMuonsMinusNtuple
    ntuplify(process,cut='Our2018',dil_name='MuonsPlusMuonsMinus') 
    # make a Our2018MuonsPlusMuonsMinusNtuple
    ntuplify(process,cut='Our2018MuPrescaled',dil_name='MuonsPlusMuonsMinus') 
    if Electrons:
        # make a Our2018MuonsElectronsAllSigns
        ntuplify(process,cut='Simple',dil_name='MuonsElectronsAllSigns') 
    #check_prescale(process, prescaled_trigger_paths)

def for_mc(process, hlt_process_name, fill_gen_info):
    process.GlobalTag.globaltag = MCGT
    # Add filters
    add_filters(process)
    # make a SimpleMuonsAllSignsNtuple
    ntuplify(process,fill_gen_info=fill_gen_info) 
    # make a Our2018MuonsPlusMuonsMinusNtuple
    ntuplify(process,cut='Our2018',dil_name='MuonsPlusMuonsMinus',fill_gen_info=fill_gen_info)
    if Electrons:
        # make a Our2018MuonsElectronsAllSigns
        ntuplify(process,cut='Simple',dil_name='MuonsElectronsAllSigns',fill_gen_info=fill_gen_info) 
    # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)
    switch_hlt_process_name(process, hlt_process_name)

def get_dataset(run):
    #JMTBAD common with dataset_details in submit below, make a DataSamples.py?
    run = int(run)
    if 190450 <= run <= 191284:
        return '/SingleMu/tucker-datamc_SingleMuRun2012A_Prompt_190450_191284_20120418134612-57b19813ab8f2ab142c4566dc6738156/USER'
    else:
        raise ValueError('dunno how to do run %i' % run)

if 'int_data' in sys.argv:
    for_data(process)
    #printify(process)
    
if 'int_mc' in sys.argv:
    for_mc(process, 'HLT', True)
    #printify(process)
    
if 'gogo' in sys.argv:
    for_data(process)
    printify(process)
    
    n = sys.argv.index('gogo')
    run, lumi, event = sys.argv[n+1], sys.argv[n+2], sys.argv[n+3]
    print run, lumi, event
    run = int(run)
    lumi = int(lumi)
    event = int(event)
    filename = [x for x in sys.argv if x.endswith('.root')]
    if filename:
        filename = filename[0]
    else:
        dataset = get_dataset(run)
        print dataset
        output = os.popen('dbs search --url https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet --query="find file where dataset=%s and run=%s and lumi=%s"' % (dataset, run, lumi)).read()
        print repr(output)
        filename = [x for x in output.split('\n') if x.endswith('.root')][0]
    print filename
    process.source.fileNames = [filename]
    from SUSYBSMAnalysis.Zprime2muAnalysis.cmsswtools import set_events_to_process
    set_events_to_process(process, [(run, event)])

#f = file('outfile_histos1', 'w')
#f.write(process.dumpPython())
#f.close()

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config,getUsernameFromSiteDB
config = config()
config.General.requestName = 'ana_datamc_%(name)s%(extra)s'
config.General.workArea = 'crab'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_crab.py'   
config.Data.inputDataset =  '%(ana_dataset)s'
config.Data.inputDBS = 'global'
config.Data.publication = False
job_control
config.Data.outputDatasetTag = 'ana_datamc_%(name)s'
config.Data.outLFNDirBase = '/store/user/'+getUsernameFromSiteDB()
config.Site.storageSite = 'T2_CH_CERN'
'''
#config.Data.outLFNDirBase = '/store/group/phys_exotica/dimuon/2018/datamc'
    
    just_testing = 'testing' in sys.argv
    extra = '_'+ex if ex!='' else ''
    # Run on data.
    if 'no_data' not in sys.argv:
        from SUSYBSMAnalysis.Zprime2muAnalysis.goodlumis import *

        dataset_details = [
            # what is this dataset? Use only run 315267?
            #('SingleMuonRun2018A-22May2018-v1', '/SingleMuon/Run2018A-22May2018-v1/MINIAOD'), 
            # To be replaced by 06June2018 (these datasets had data deletion problems)
            #('SingleMuonRun2018A-PromptReco-v1', '/SingleMuon/Run2018A-PromptReco-v1/MINIAOD'), 
            #('SingleMuonRun2018A-PromptReco-v2', '/SingleMuon/Run2018A-PromptReco-v2/MINIAOD'),

            # PromptReco A-C
            # PPD recommendation for 2018A PromptReco 
            # 06Jun2018-v1 + PromptReco-v3
            #('SingleMuonRun2018A-06June2018-v1', '/SingleMuon/Run2018A-06Jun2018-v1/MINIAOD'), 
            #('SingleMuonRun2018A-PromptReco-v3', '/SingleMuon/Run2018A-PromptReco-v3/MINIAOD'),
            #('SingleMuonRun2018B-PromptReco-v1', '/SingleMuon/Run2018B-PromptReco-v1/MINIAOD'),
            #('SingleMuonRun2018B-PromptReco-v2', '/SingleMuon/Run2018B-PromptReco-v2/MINIAOD'),
            #('SingleMuonRun2018C-PromptReco-v1', '/SingleMuon/Run2018C-PromptReco-v1/MINIAOD'),
            #('SingleMuonRun2018C-PromptReco-v2', '/SingleMuon/Run2018C-PromptReco-v2/MINIAOD'),
            #('SingleMuonRun2018C-PromptReco-v3', '/SingleMuon/Run2018C-PromptReco-v3/MINIAOD'),
            # Good to use
            ('SingleMuonRun2018A-17Sep2018-v2',  '/SingleMuon/Run2018A-17Sep2018-v2/MINIAOD'),
            ('SingleMuonRun2018B-17Sep2018-v1',  '/SingleMuon/Run2018B-17Sep2018-v1/MINIAOD'),
            ('SingleMuonRun2018C-17Sep2018-v1',  '/SingleMuon/Run2018C-17Sep2018-v1/MINIAOD'),
            ('SingleMuonRun2018D-PromptReco-v2', '/SingleMuon/Run2018D-PromptReco-v2/MINIAOD'),
        ]

        lumi_lists = ['Run2018MuonsOnly']

        jobs = []
        for lumi_name in lumi_lists:
            ll = eval(lumi_name + '_ll') if lumi_name != 'NoLumiMask' else None
            for dd in dataset_details:
                jobs.append(dd + (lumi_name, ll))
                
        for dataset_name, ana_dataset, lumi_name, lumi_list in jobs:

            json_fn = 'tmp.json'
            lumi_list.writeJSON(json_fn)
            lumi_mask = json_fn

            name = '%s_%s' % (lumi_name, dataset_name)
            print name

            new_py = open('histos.py').read()
            new_py += "\nfor_data(process)\n"
            if '17Sep2018' in dataset_name:
                new_py += "\nprocess.GlobalTag.globaltag = '102X_dataRun2_Sep2018Rereco_v1'\n"
            else:
                new_py += "\nprocess.GlobalTag.globaltag = '102X_dataRun2_Prompt_v11'\n"
            open('histos_crab.py', 'wt').write(new_py)

            new_crab_cfg = crab_cfg % locals()

            job_control = '''
config.Data.splitting = 'Automatic'
config.Data.lumiMask = '%(lumi_mask)s'
''' % locals()

            new_crab_cfg = new_crab_cfg.replace('job_control', job_control)
            open('crabConfig.py', 'wt').write(new_crab_cfg)

            if not just_testing:
                os.system('crab submit -c crabConfig.py')
            else:
                cmd = 'diff histos.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crab.py'
                print cmd
                os.system(cmd)

        if not just_testing:
            os.system('rm crabConfig.py crabConfig.pyc histos_crab.py histos_crab.pyc tmp.json')

    if 'no_mc' not in sys.argv:
        # Set crab_cfg for MC.
        crab_cfg = crab_cfg.replace('job_control','''
config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 10000
    ''')

       
        for sample in samples:
            name = sample.name
            ana_dataset = sample.dataset
            print name, ana_dataset

            new_py = open('histos.py').read()
            new_py += "\nfor_mc(process,'HLT',True)\n"
            new_py += "\napply_gen_filters(process,\"%(name)s\")\n"%locals()
            open('histos_crab.py', 'wt').write(new_py)

            open('crabConfig.py', 'wt').write(crab_cfg % locals())
            if not just_testing:
                os.system('crab submit -c crabConfig.py')
            else:
                cmd = 'diff histos.py histos_crab.py | less'
                print cmd
                os.system(cmd)
                cmd = 'less crabConfig.py'
                print cmd
                os.system(cmd)

        if not just_testing:
            os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')

