import sys, os, shutil
if __name__ == '__main__':

  dirs = [
                ('noWei_noRand_dy120','noWei_noRand/crab_AccXEff_20171010/crab_AccXEff_noWei_noRand_20171010_dy120'),
                ('noWei_noRand_dy200','noWei_noRand/crab_AccXEff_20171010/crab_AccXEff_noWei_noRand_20171010_dy200'),
                ('noWei_noRand_dy400','noWei_noRand/crab_AccXEff_20171010/crab_AccXEff_noWei_noRand_20171010_dy400'),
                ('noWei_noRand_dy800','noWei_noRand/crab_AccXEff_20171010/crab_AccXEff_noWei_noRand_20171010_dy800'),
                ('noWei_noRand_dy1400','noWei_noRand/crab_AccXEff_20171010/crab_AccXEff_noWei_noRand_20171010_dy1400'),
                ('noWei_noRand_dy2300','noWei_noRand/crab_AccXEff_20171010/crab_AccXEff_noWei_noRand_20171010_dy2300'),
                ('noWei_noRand_dy3500','noWei_noRand/crab_AccXEff_20171010/crab_AccXEff_noWei_noRand_20171010_dy3500'),
                ('noWei_noRand_dy4500','noWei_noRand/crab_AccXEff_20171010/crab_AccXEff_noWei_noRand_20171010_dy4500'),
                ('noWei_noRand_dy6000','noWei_noRand/crab_AccXEff_20171010/crab_AccXEff_noWei_noRand_20171010_dy6000'),
                ('noWei_noRand_dyInf','noWei_noRand/crab_AccXEff_20171010/crab_AccXEff_noWei_noRand_20171010_dyInf'),

                ('yesWei_noRand_dy120','yesWei_noRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_noRand_20171010_dy120'),
                ('yesWei_noRand_dy200','yesWei_noRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_noRand_20171010_dy200'),
                ('yesWei_noRand_dy400','yesWei_noRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_noRand_20171010_dy400'),
                ('yesWei_noRand_dy800','yesWei_noRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_noRand_20171010_dy800'),
                ('yesWei_noRand_dy1400','yesWei_noRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_noRand_20171010_dy1400'),
                ('yesWei_noRand_dy2300','yesWei_noRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_noRand_20171010_dy2300'),
                ('yesWei_noRand_dy3500','yesWei_noRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_noRand_20171010_dy3500'),
                ('yesWei_noRand_dy4500','yesWei_noRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_noRand_20171010_dy4500'),
                ('yesWei_noRand_dy6000','yesWei_noRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_noRand_20171010_dy6000'),
                ('yesWei_noRand_dyInf','yesWei_noRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_noRand_20171010_dyInf'),

                ('yesWei_yesRand_dy120','yesWei_yesRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_yesRand_20171010_dy120'),
                ('yesWei_yesRand_dy200','yesWei_yesRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_yesRand_20171010_dy200'),
                ('yesWei_yesRand_dy400','yesWei_yesRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_yesRand_20171010_dy400'),
                ('yesWei_yesRand_dy800','yesWei_yesRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_yesRand_20171010_dy800'),
                ('yesWei_yesRand_dy1400','yesWei_yesRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_yesRand_20171010_dy1400'),
                ('yesWei_yesRand_dy2300','yesWei_yesRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_yesRand_20171010_dy2300'),
                ('yesWei_yesRand_dy3500','yesWei_yesRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_yesRand_20171010_dy3500'),
                ('yesWei_yesRand_dy4500','yesWei_yesRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_yesRand_20171010_dy4500'),
                ('yesWei_yesRand_dy6000','yesWei_yesRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_yesRand_20171010_dy6000'),
                ('yesWei_yesRand_dyInf','yesWei_yesRand/crab_AccXEff_20171010/crab_AccXEff_yesWei_yesRand_20171010_dyInf'),

          ]


  for name, di in dirs:

    #name = '|| -- '+name+' -- ||'
    print '\n',name
    if 'mergeoutputs' in sys.argv:
      thePath = '/u/user/msoh/ZPrime/AccXEff/CMSSW_8_0_26/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/EfficiencyResolutionFromMC/AccXEff/' + di + '/results'
      os.chdir(thePath)
      print '\n',os.getcwd()
      theFileName = 'histos_' + name + '.root'
      hadd_command = 'hadd ' + theFileName + ' zp2mu_histos*.root'
      os.system(hadd_command)

      outputPath = '/u/user/msoh/ZPrime/AccXEff/CMSSW_8_0_26/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/EfficiencyResolutionFromMC/AccXEff/outputs'
      shutil.copy2(thePath + '/' + theFileName, outputPath + '/' +theFileName)
      os.chdir('/u/user/msoh/ZPrime/AccXEff/CMSSW_8_0_26/src/SUSYBSMAnalysis/Zprime2muAnalysis/test/EfficiencyResolutionFromMC/AccXEff/')

    print '\n\n\tfinish!\n'
