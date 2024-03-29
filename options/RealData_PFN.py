#--   622 events - 12 files
#--   Data from first collisions 23-11-09
##--  Energy per beam:  0.45 TeV 

from Gaudi.Configuration import * 

EventSelector().Input   = [
    #### First Collisions - Does not seem to work for now but can try again - Magnet off ### 
    ##"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005616/0000/00005616_00000002_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    ##"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005616/0000/00005616_00000003_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    ##"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005616/0000/00005616_00000005_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    ##"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005616/0000/00005616_00000006_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    ##"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005616/0000/00005616_00000007_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    ##"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005616/0000/00005616_00000008_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    ##"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005616/0000/00005616_00000009_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    ##"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005616/0000/00005616_00000010_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    ##"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005616/0000/00005616_00000011_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    ##"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005616/0000/00005616_00000012_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    ##"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005616/0000/00005617_00000001_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    ##"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005616/0000/00005617_00000002_1.dst' TYP='POOL_ROOTTREE' OPT='READ'"#,

    ### More collisions, magnet on. 264837 events ###
   #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005717/0000/00005717_00000010_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005717/0000/00005717_00000008_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005717/0000/00005717_00000007_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005717/0000/00005717_00000005_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005717/0000/00005717_00000004_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005717/0000/00005717_00000001_1.dst' TYP='POOL_ROOTTREE' OPT='READ'"#,
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000037_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000035_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000018_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000013_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000012_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000011_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000010_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000009_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000008_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000007_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000006_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000005_1.dst' TYP='POOL_ROOTTREE' OPT='READ'", 
   #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000004_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000003_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000002_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005710/0000/00005710_00000001_1.dst' TYP='POOL_ROOTTREE' OPT='READ'"#,
    
    ### More collisions - run 5727 ###
    "DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005727/0000/00005727_00000031_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    "DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005727/0000/00005727_00000007_1.dst' TYP='POOL_ROOTTREE' OPT='READ'"#,
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005842/0000/00005842_00000194_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005842/0000/00005842_00000188_1.dst' TYP='POOL_ROOTTREE' OPT='READ'"
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005727/0000/00005727_00000008_1.dst' TYP='POOL_ROOTTREE' OPT='READ'"#,
    ### RecoToDST-07 ###
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005842/0000/00005842_00000197_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    

    ### Even more collisions - This time with beam energy of 1.18TeV mag down, VELO open ###
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005785/0000/00005785_00000001_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='PFN:/castor/cern.ch/grid/lhcb/data/2009/DST/00005785/0000/00005785_00000002_1.dst' TYP='POOL_ROOTTREE' OPT='READ'"#,

    ### LFNs ###
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000044_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000043_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000042_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000041_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000040_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000039_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000038_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000037_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000036_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000035_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000034_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000033_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000032_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000031_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000030_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000029_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000028_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000027_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000026_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000025_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000024_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000023_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000022_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000021_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000020_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000019_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000018_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000017_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000016_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000015_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000014_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000013_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000012_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000011_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000010_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000008_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000007_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000005_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000004_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
   #"DATAFILE='LFN:/lhcb/data/2009/DST/00005717/0000/00005717_00000001_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000037_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000035_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000018_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000013_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000012_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000011_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000010_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000009_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000008_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000007_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000006_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000005_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000004_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000003_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000002_1.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    #"DATAFILE='LFN:/lhcb/data/2009/DST/00005710/0000/00005710_00000001_1.dst' TYP='POOL_ROOTTREE' OPT='READ'"#,

    ]
