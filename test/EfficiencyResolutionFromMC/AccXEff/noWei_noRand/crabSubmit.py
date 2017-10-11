from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys, os
import datetime

now = datetime.datetime.now()
date = now.strftime('%Y%m%d')

submitVersion = 'AccXEff_R/%s' % (date)### modified
mainOutputDir = '/store/user/%s/%s' % (getUsernameFromSiteDB(),submitVersion) ### modified

if __name__ == '__main__' and 'submit' in sys.argv:
    crab_cfg = '''
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'AccXEff_noWei_noRand_%(date)s_%(name)s'
config.General.workArea = 'crab_AccXEff_%(date)s'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'histos_crab.py'
#config.JobType.priority = 1

config.Data.outLFNDirBase = '%(mainOutputDir)s/'
config.Data.inputDataset =  '%(dataset)s'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'FileBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 20000
#config.Data.unitsPerJob  = 100
config.Data.publication = False
config.Data.outputDatasetTag = 'AccXEff_noWei_noRand_%(date)s_%(name)s'
#config.Data.outLFNDirBase = '/store/user/moh'
config.Data.ignoreLocality = True

config.Site.storageSite = 'T2_KR_KNU'

'''

    samples = [
             #  ('dy120','/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',50, 120),
             #  ('dy200','/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',120, 200),
             #  ('dy400','/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',200, 400),
             #  ('dy800','/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',400, 800),
             #  ('dy1400','/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',800, 1400),
             #  ('dy2300','/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',1400, 2300),
             #  ('dy3500','/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',2300,3500),
             #  ('dy4500','/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',3500,4500),
             #  ('dy6000','/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',4500,6000),
             #  ('dyInf','/ZToMuMu_NNPDF30_13TeV-powheg_M_6000_Inf/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM',6000,-1),

                ('dy120','/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',50,120),
                ('dy200','/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',120,200),
                ('dy400','/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',200,400),
                ('dy800','/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',400,800),
                ('dy1400','/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',800,1400),
                ('dy2300','/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',1400,2300),
                ('dy3500','/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',2300,3500),
                ('dy4500','/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',3500,4500),
                ('dy6000','/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',4500,6000),
                ('dyInf','/ZToMuMu_NNPDF30_13TeV-powheg_M_6000_Inf/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',6000,1e99),

       ]

    just_testing = 'testing' in sys.argv

    for name, dataset, lo, hi in samples:
        open('crabConfig.py', 'wt').write(crab_cfg % locals())

        new_py = open('histos_eff.py').read()
        new_py += '\nprocess.HardInteractionFilter.min_mass = "%i"\n' % lo
        new_py += '\nprocess.HardInteractionFilter.max_mass = "%i"\n' % hi
        new_py += '\nprocess.HardInteractionFilterRes.min_mass = "%i"\n' % lo
        new_py += '\nprocess.HardInteractionFilterRes.max_mass = "%i"\n' % hi
        open('histos_crab.py', 'wt').write(new_py)

        if not just_testing:
            os.system('crab submit  -c crabConfig.py')
            #os.system('rm crabConfig.py histos_crab.py histos_crab.pyc')
