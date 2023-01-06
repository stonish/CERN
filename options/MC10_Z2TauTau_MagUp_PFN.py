# MC10 DSTs for MagUp data sample
# Z2TauTau with tau-> anything
# List produced on 8th April 2011 (Total of 44 DSTs available at this stage)

from Gaudi.Configuration import *

EventSelector().Input   = [
    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00009068/0000/00009068_00000041_1.allstreams.dst' TYP='POOL_ROOTTREE' Opt='READ'",
    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00009068/0000/00009068_00000011_1.allstreams.dst' TYP='POOL_ROOTTREE' Opt='READ'",
    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00009068/0000/00009068_00000018_1.allstreams.dst' TYP='POOL_ROOTTREE' Opt='READ'",
    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00009068/0000/00009068_00000003_1.allstreams.dst' TYP='POOL_ROOTTREE' Opt='READ'",
    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00009068/0000/00009068_00000032_1.allstreams.dst' TYP='POOL_ROOTTREE' Opt='READ'"
    ]

