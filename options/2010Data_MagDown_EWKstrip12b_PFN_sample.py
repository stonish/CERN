from Gaudi.Configuration import * 

EventSelector().Input   = [
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00008380/0000/00008380_00000177_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00008380/0000/00008380_00000160_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00008380/0000/00008380_00000033_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00008380/0000/00008380_00000076_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00008380/0000/00008380_00000088_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00008380/0000/00008380_00000187_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00008380/0000/00008380_00000002_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00008380/0000/00008380_00000032_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'",
"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00008380/0000/00008380_00000724_1.ew.dst' TYP='POOL_ROOTTREE' OPT='READ'"]
FileCatalog().Catalogs = [ 'xmlcatalog_file:2010Data_MagDown_EWKstrip12b_PFN_sample.xml' ]
