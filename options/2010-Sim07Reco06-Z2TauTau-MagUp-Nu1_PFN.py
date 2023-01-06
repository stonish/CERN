#-- GAUDI jobOptions generated on Fri Nov  5 15:48:13 2010
#-- Contains event types : 
#--   42100000 - 3 files - 54200 events - 13.53 GBytes

from Gaudi.Configuration import * 

EventSelector().Input   = [
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/2010/DST/00007700/0000/00007700_00000001_1.dst?svcClass=lhcbdata' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/2010/DST/00007700/0000/00007700_00000002_1.dst?svcClass=lhcbdata' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/2010/DST/00007700/0000/00007700_00000003_1.dst?svcClass=lhcbdata' TYP='POOL_ROOTTREE' OPT='READ'"]
FileCatalog().Catalogs = [ 'xmlcatalog_file:2010-Sim07Reco06-Z2TauTau-MagUp-Nu1_PFN.xml' ]
