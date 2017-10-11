import sys, os
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

    name = '|| -- '+name+' -- ||'
    if 'status' in sys.argv:
      command = 'crab status -d ' + di + ' --verboseErrors --long'
    elif 'resubmit' in sys.argv:
      command = 'crab resubmit -d ' + di
    elif 'getoutput' in sys.argv:
      command = 'crab getoutput -d ' + di + ' --quantity=all'
    elif 'kill' in sys.argv:
      command = 'crab kill -d ' + di

    print '\n',name
    os.system(command)
