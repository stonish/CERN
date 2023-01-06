from Gaudi.Configuration import * 

EventSelector().Input   = [
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/LHCb/Collision11/EW.DST/00010830/0000/00010830_00002156_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/LHCb/Collision11/EW.DST/00010830/0000/00010830_00002098_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/LHCb/Collision11/EW.DST/00010830/0000/00010830_00002045_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/LHCb/Collision11/EW.DST/00010830/0000/00010830_00002002_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
#"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/LHCb/Collision11/EW.DST/00010830/0000/00010830_00001925_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",#File can no longer be found at CERN
#"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/LHCb/Collision11/EW.DST/00010830/0000/00010830_00001879_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",#File can no longer be found at CERN
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/LHCb/Collision11/EW.DST/00010830/0000/00010830_00001826_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/LHCb/Collision11/EW.DST/00010830/0000/00010830_00001766_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/LHCb/Collision11/EW.DST/00010830/0000/00010830_00001630_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'"]
FileCatalog().Catalogs = [ 'xmlcatalog_file:2011Data_MagDown_EWKstrip13b_PFN_sample.xml' ]
