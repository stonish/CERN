#-- GAUDI jobOptions generated on Fri Nov  5 15:44:20 2010
#-- Contains event types : 
#--   42100000 - 3 files - 53000 events - 13.27 GBytes

from Gaudi.Configuration import * 

EventSelector().Input   = [
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/2010/DST/00007724/0000/00007724_00000001_1.dst?svcClass=lhcbdata' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/2010/DST/00007724/0000/00007724_00000002_1.dst?svcClass=lhcbdata' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/2010/DST/00007724/0000/00007724_00000003_1.dst?svcClass=lhcbdata' TYP='POOL_ROOTTREE' OPT='READ'"]
FileCatalog().Catalogs = [ 'xmlcatalog_file:2010-Sim07Reco06-Z2TauTau-MagDown-Nu1_PFN.xml' ]
