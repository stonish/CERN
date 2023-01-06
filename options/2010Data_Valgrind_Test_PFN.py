# Just ten (one) DST(s) used to run a large (small) Valgrind job to check for memory leaks etc.
# File is from Strip12b, MagUp

from Gaudi.Configuration import *

EventSelector().Input   = [
    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00008380/0000/00008380_00000696_1.ew.dst' TYP='POOL_ROOTTREE' Opt='\
    READ'"] # A DST which will not work for some reason!
    #"   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00008394/0000/00008394_00000302_1.ew.dst' TYP='POOL_ROOTTREE' Opt='READ'", #Extra - Check "illegal leaves" in motherTuple and recAllTuple
#    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00007958/0000/00007958_00000174_1.ew.dst' TYP='POOL_ROOTTREE' Opt='READ'",
#    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00007958/0000/00007958_00000313_1.ew.dst' TYP='POOL_ROOTTREE' Opt='READ'",
#    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00007958/0000/00007958_00000214_1.ew.dst' TYP='POOL_ROOTTREE' Opt='READ'",
#    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00007958/0000/00007958_00000182_1.ew.dst' TYP='POOL_ROOTTREE' Opt='READ'",
#    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00007958/0000/00007958_00000312_1.ew.dst' TYP='POOL_ROOTTREE' Opt='READ'",
#    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00007958/0000/00007958_00000055_1.ew.dst' TYP='POOL_ROOTTREE' Opt='READ'",
#    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00007958/0000/00007958_00000282_1.ew.dst' TYP='POOL_ROOTTREE' Opt='READ'",
#    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00007958/0000/00007958_00000008_1.ew.dst' TYP='POOL_ROOTTREE' Opt='READ'"]

#   "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/data/2010/EW.DST/00008394/0000/00008394_00000925_1.ew.dst' TYP='POOL_ROOTTREE' Opt='READ'"]

