#--   ~200,000 events - 12 files
#--   Data from first collisions, starting 23-11-09

from Gaudi.Configuration import * 

EventSelector().Input   = [
    #### Full run data 2009, cleaned up to include only interesting events (reduced from ~3.5e6 events to ~200,000); Beam energy = 0.45TeV  ### 
    "DATAFILE='PFN:/castor/cern.ch/user/r/rhuston/Data/2009/5727/0/Data2009.PhysEvent.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    "DATAFILE='PFN:/castor/cern.ch/user/r/rhuston/Data/2009/5727/1/Data2009.PhysEvent.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    "DATAFILE='PFN:/castor/cern.ch/user/r/rhuston/Data/2009/5727/2/Data2009.PhysEvent.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    "DATAFILE='PFN:/castor/cern.ch/user/r/rhuston/Data/2009/5727/3/Data2009.PhysEvent.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    "DATAFILE='PFN:/castor/cern.ch/user/r/rhuston/Data/2009/5727/4/Data2009.PhysEvent.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    "DATAFILE='PFN:/castor/cern.ch/user/r/rhuston/Data/2009/5727/5/Data2009.PhysEvent.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    "DATAFILE='PFN:/castor/cern.ch/user/r/rhuston/Data/2009/5727/6/Data2009.PhysEvent.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    "DATAFILE='PFN:/castor/cern.ch/user/r/rhuston/Data/2009/5727/7/Data2009.PhysEvent.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    "DATAFILE='PFN:/castor/cern.ch/user/r/rhuston/Data/2009/5727/8/Data2009.PhysEvent.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    "DATAFILE='PFN:/castor/cern.ch/user/r/rhuston/Data/2009/5727/9/Data2009.PhysEvent.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    "DATAFILE='PFN:/castor/cern.ch/user/r/rhuston/Data/2009/5727/10/Data2009.PhysEvent.dst' TYP='POOL_ROOTTREE' OPT='READ'",
    "DATAFILE='PFN:/castor/cern.ch/user/r/rhuston/Data/2009/5727/11/Data2009.PhysEvent.dst' TYP='POOL_ROOTTREE' OPT='READ'"#,

    ]
