#!/usr/bin/env python
import FWCore.ParameterSet.Config as cms
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_hlt_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import switch_reco_process_name
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cfg import process
from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import goodDataFiltersMiniAOD

process.source.fileNames =['/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/EB675885-31D1-9448-A802-4D8AE9B37566.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/51B1931F-8D12-8C4A-87BB-3B1FE0502C6C.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/CD05A289-9842-754A-98E6-148388C4EF41.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/9703A373-2475-4F45-BA2B-30CBDDA74298.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/8E033556-7497-5341-AC31-02CCB9CDCADF.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/543B865A-C341-EF4F-8DD6-D56D6DFEBCFA.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/89E5100A-81A9-E64F-A2DB-C378D11A6C64.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/9898EABB-AFFE-AA41-9355-26DD08016FDC.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/99BD5B79-F881-6D40-8970-01B0088F74D1.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/4F3043C6-B0F4-C046-A075-9EB1174C26AF.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/1A885AD3-3BDE-AA43-88FB-CA7AE2F88BF2.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/3E77D9B6-3963-754A-9110-EEC969468B5E.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/3A00FE0E-09F5-6148-8A43-53BEBD2F39B8.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/802AE29A-757E-1B41-8182-D074726C27FB.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/90842FB5-209B-9B4E-8833-1EE978B9D4DB.root','/store/mc/RunIIAutumn18MiniAOD/ZToEE_NNPDF31_TuneCP5_13TeV-powheg_M_1400_2300/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/000C927C-4DC5-E94E-9FFD-DC522900611D.root',]

process.maxEvents.input = -1
isMC = True
addNTuples = False
year = 2018
sampleName = 'dy1400to2300'
process.GlobalTag.globaltag = '94X_mc2017_realistic_v17'
process.options.wantSummary = cms.untracked.bool(True)# false di default
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # default 1000
#import for high pT muon triggers
from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import trigger_match, trigger_paths, overall_prescale, offline_pt_threshold, trigger_filters, trigger_filters2016, trigger_path_names, trigger_path_names2016, trigger_match_2018, trigger_match_2016,trigger_path_full_names, trigger_path_full_names2016
#import for prescaled low pT muon triggers
from SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi import prescaled_trigger_pt_threshold, overall_prescale_2016, overall_prescale_2017, overall_prescale_2018, prescaled_trigger_filters_16, prescaled_trigger_path_names_16, prescaled_trigger_path_full_names_16, prescaled_trigger_match_2016, prescaled_trigger_filters_18, prescaled_trigger_path_names_18, prescaled_trigger_path_full_names_18, prescaled_trigger_match_2018, prescaled_trigger_path_name_list_16, prescaled_trigger_path_name_list_17, prescaled_trigger_path_name_list_18, prescaled_offline_pt_threshold

# The histogramming module that will be cloned multiple times below
# for making histograms with different cut/dilepton combinations.

from SUSYBSMAnalysis.Zprime2muAnalysis.Zprime2muAnalysis_cff import electrons_miniAOD
electrons_miniAOD(process)
#!/usr/bin/env python
from SUSYBSMAnalysis.Zprime2muAnalysis.HistosFromPAT_cfi import HistosFromPAT_MiniAOD as HistosFromPAT
HistosFromPAT.leptonsFromDileptons = True
HistosFromPAT.usekFactor = True #### Set TRUE to use K Factor on DY. If used, the k factor will be applied to ALL samples submitted. #####
HistosFromPAT.useTTBarWeight = False #### Set TRUE to use NNPDF Weights for ttbar #####
	
# These modules define the basic selection cuts. For the monitoring
# sets below, we don't need to define a whole new module, since they
# just change one or two cuts -- see below.
import SUSYBSMAnalysis.Zprime2muAnalysis.ElectronSelection_cff as ElectronSelection
import SUSYBSMAnalysis.Zprime2muAnalysis.ElectronSelection2016_cff as ElectronSelection2016
import SUSYBSMAnalysis.Zprime2muAnalysis.ElectronSelection2018_cff as ElectronSelection2018

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
if year == 2016:
    setupEgammaPostRecoSeq(process, era='2016-Legacy')
if year == 2017:
    setupEgammaPostRecoSeq(process, era='2017-Nov17ReReco')
if year == 2018:
	print "setting up"
	setupEgammaPostRecoSeq(process, era='2018-Prompt')

from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
if year == 2016 or year == 2017:
	process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
    		DataEra = cms.string("2017BtoF"), #Use 2016BtoH for 2016
    		UseJetEMPt = cms.bool(False),
    		PrefiringRateSystematicUncty = cms.double(0.2),
    		SkipWarnings = False)

	if year==2016:
    		process.prefiringweight.DataEra = cms.string("2016BtoH")

dils = [
#	('ElectronsOppSign',        '%(leptons_name)s:electrons@+ %(leptons_name)s:electrons@-',     ''),
#	('ElectronsSameSign',       '%(leptons_name)s:electrons@+ %(leptons_name)s:electrons@+',     ''),
	('ElectronsAllSigns',       '%(leptons_name)s:electrons@+ %(leptons_name)s:electrons@+',     ''),
	]

cuts = {
	'ElectronSelection'  : ElectronSelection,
	}
	
if year == 2016:
	cuts = {
	"ElectronSelection" : ElectronSelection2016,
	}
if year == 2018:
	cuts = {
	"ElectronSelection" : ElectronSelection2018,
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
    # present
    path_list.append(process.egmGsfElectronIDSequence)
	    
    leptons_name = cut_name + 'Leptons'
    leptons = process.leptonsMini.clone()
    if year == 2016 and ("03Feb" in sampleName or "23Sep" in sampleName or "Prompt" in sampleName or ("dy" in sampleName  and not "Inclusive" in sampleName) or "CI" in sampleName or "ADD" in sampleName) and not sampleName == "dyMCAtNLO":
    #if year == 2016:
	leptons.trigger_summary = cms.InputTag('selectedPatTrigger')
    if year == 2018:
	leptons.hlt_filter_ele = cms.vstring('hltDiEle25CaloIdLMWPMS2UnseededFilter')
    if year == 2016:
	leptons.hlt_filter_ele = cms.vstring('hltDiEle33CaloIdLMWPMS2UnseededFilter','hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter')


    leptons.electron_src = cms.InputTag('slimmedElectrons',"","Zprime2muAnalysis")
    setattr(process, leptons_name, leptons)
    path_list.append(leptons)
	
    # Make all the combinations of dileptons we defined above.
    for dil_name, dil_decay, dil_cut in dils:

        # Unique names for the modules: allname for the allDileptons,
        # and name for dileptons.
        name = cut_name + dil_name
        allname = 'all' + name

        alldil = Selection.allDielectrons.clone(decay = dil_decay % locals(), cut = dil_cut)
        if 'AllSigns' in dil_name:
            alldil.checkCharge = cms.bool(False)
        dil = Selection.dielectrons.clone(src = cms.InputTag(allname))
	
	# Implement the differences to the selections; currently, as
        # in Zprime2muCombiner, the cuts in loose_cut and
        # tight_cut are the ones actually used to drop leptons, and
        # not the ones passed into the LeptonProducer to set cutFor above.
        # Histos now just needs to know which leptons and dileptons to use.
	if isMC:
		if year == 2018:
			histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'electrons'), dilepton_src = cms.InputTag(name),doElectrons = cms.bool(True),pu_weights = cms.vstring("dy1400to2300","data_2018"), year = cms.int32(year))
		elif year == 2017:	
			histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'electrons'), dilepton_src = cms.InputTag(name),doElectrons = cms.bool(True),pu_weights = cms.vstring("mc_2017","data_2017"), year = cms.int32(year))
		else:	
			histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'electrons'), dilepton_src = cms.InputTag(name),doElectrons = cms.bool(True),pu_weights = cms.vstring("mc_2016","data_2016"), year = cms.int32(year))
	else:	
		histos = HistosFromPAT.clone(lepton_src = cms.InputTag(leptons_name, 'electrons'), dilepton_src = cms.InputTag(name),doElectrons = cms.bool(True))
	#if not isMC:
	#	delattr(histos,'hardInteraction')

	histos.hardInteraction.doingElectrons = True

	if 'ConLR' in sampleName or 'DesLR' in sampleName or 'ConRL' in sampleName or 'DesRL' in sampleName:
		L = 10000	
		if '16TeV' in sampleName:
			L = 16000
		if '100kTeV' in sampleName:
			L = 100000000
		if '1TeV' in sampleName:
			L = 1000
		if '22TeV' in sampleName:
			L = 22000
		if '24TeV' in sampleName:
			L = 24000
		if '28TeV' in sampleName:
			L = 28000
		if '32TeV' in sampleName:
			L = 32000
		if '34TeV' in sampleName:
			L = 34000
		if '40TeV' in sampleName:
			L = 40000
		if "Con" in sampleName:
			histos.lrWeightProducer.interference = -1
		else:	
			histos.lrWeightProducer.interference = 1
		histos.lrWeightProducer.Lambda = L	
		histos.lrWeightProducer.calculate = True
		histos.lrWeightProducer.doingElectrons = True
		if "RL" in sampleName:
			histos.lrWeightProducer.doingLR = False
	
        # Add all these modules to the process and the path list.
        setattr(process, allname, alldil)
        setattr(process, name, dil)
        setattr(process, name + 'Histos', histos)
	if not isMC:
		#del histos.hardInteraction
		#histos.useMadgraphWeight  = False
		
		trig = Selection.dielectronHLT
		if "DoubleEG2016H" in sampleName:
			trig.triggerConditions = cms.vstring("HLT_DoubleEle33_CaloIdL_MW_v*")
		trigName = cut_name + "HLTFilter"
        	setattr(process, trigName, trig)
		delattr(getattr(process,name + 'Histos'),'hardInteraction')	
        	path_list.append(trig * alldil * dil * histos)
	else:
		if year == 2018:
			alldil.loose_cut_ele = cms.string('et > 35 && abs(userFloat("etaSC")) < 2.5 && !(abs(userFloat("etaSC")) > 1.4442 && abs(userFloat("etaSC")) < 1.566) && userInt("cutFor2018") == 1')	
		else:
			alldil.loose_cut_ele = cms.string('et > 35 && abs(userFloat("etaSC")) < 2.5 && !(abs(userFloat("etaSC")) > 1.4442 && abs(userFloat("etaSC")) < 1.566) && userInt("cutFor") == 1')	
		alldil.tight_cut_ele = cms.string("")
		#if not year == 2017:	
		alldil.ele_match_l1 = cms.bool(False)
		if "CI" in sampleName or "ADD" in sampleName:
			getattr(process, name + 'Histos', histos).hardInteraction.matchTaus = cms.bool(False)
        	path_list.append(alldil * dil * histos)
	
    # Finally, make the path for this set of cuts.
    pathname = 'path' + cut_name
    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.DielectronPreselector_cfi')
    process.load("SUSYBSMAnalysis.Zprime2muAnalysis.EventCounter_cfi")
    if year == 2016 or year == 2017:
    	pobj = process.EventCounter * process.egammaPostRecoSeq * process.prefiringweight * process.dielectronPreseletor *  process.muonPhotonMatchMiniAOD * reduce(lambda x,y: x*y, path_list)
    else:	
	pobj = process.EventCounter * process.egammaPostRecoSeq * process.dielectronPreseletor *  process.muonPhotonMatchMiniAOD * reduce(lambda x,y: x*y, path_list)

    path = cms.Path(pobj)
    setattr(process, pathname, path)

    process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
    for dataFilter in goodDataFiltersMiniAOD:
	#setattr(process,dataFilter 
	getattr(process,pathname).insert(1,dataFilter)



if addNTuples:
	  
	process.SimpleNtupler = cms.EDAnalyzer('SimpleNtupler_miniAOD',
					   dimu_src = cms.InputTag('ElectronSelectionElectronsAllSigns'),
					   met_src = cms.InputTag("slimmedMETs"),
					   jet_src = cms.InputTag("slimmedJets"),
					   beamspot_src = cms.InputTag('offlineBeamSpot'),
					   vertices_src = cms.InputTag('offlineSlimmedPrimaryVertices'),
					   #TriggerResults_src = cms.InputTag('TriggerResults', '', 'PAT'),	#mc
					   TriggerResults_src = cms.InputTag('TriggerResults', '', 'RECO'),	#data
					   genEventInfo = cms.untracked.InputTag('generator'),
					   metFilter = cms.VInputTag( cms.InputTag("Flag_HBHENoiseFilter"), cms.InputTag("Flag_HBHENoiseIsoFilter"), cms.InputTag("Flag_EcalDeadCellTriggerPrimitiveFilter"), cms.InputTag("Flag_eeBadScFilter"), cms.InputTag("Flag_globalTightHalo2016Filter")),
					   doElectrons = cms.bool(True),
					   )


	if isMC:
		#process.load('SUSYBSMAnalysis.Zprime2muAnalysis.PrunedMCLeptons_cfi')
		#obj = process.prunedMCLeptons
		#obj.src = cms.InputTag('prunedGenParticles')
	 
		from SUSYBSMAnalysis.Zprime2muAnalysis.HardInteraction_cff import hardInteraction_MiniAOD as hardInteraction
		hardInteraction.doingElectrons = True
		process.SimpleNtupler.hardInteraction = hardInteraction
		if hasattr(process, 'pathElectronSelection'):
			#process.pathElectronSelection *=obj * process.SimpleNtupler 
			process.pathElectronSelection *= process.SimpleNtupler 

	else:
		if hasattr(process, 'pathElectronSelection'):
			process.pathElectronSelection *= process.SimpleNtupler 


if isMC:
	switch_reco_process_name(process, "PAT") # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)
    #switch_hlt_process_name(process, hlt_process_name) # this must be done last (i.e. after anything that might have an InputTag for something HLT-related)