
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'AccXEff_noWei_noRand_20171010_dyInf'
config.General.workArea = 'crab_AccXEff_20171010'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_crab.py'
#config.JobType.priority = 1

config.Data.outLFNDirBase = '/store/user/moh/AccXEff_R/20171010/'
config.Data.inputDataset =  '/ZToMuMu_NNPDF30_13TeV-powheg_M_6000_Inf/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 20000
#config.Data.unitsPerJob  = 100
config.Data.publication = False
config.Data.outputDatasetTag = 'AccXEff_noWei_noRand_20171010_dyInf'
#config.Data.outLFNDirBase = '/store/user/moh'
config.Data.ignoreLocality = True

config.Site.storageSite = 'T2_KR_KNU'

