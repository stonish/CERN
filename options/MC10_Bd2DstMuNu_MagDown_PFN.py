# MC10 DSTs for MagDown data sample
# Bd2DstMuNu 
# List produced on 19th April 2011 (Total of 84 DSTs (~10^9 events) available at this stage)

from Gaudi.Configuration import *

EventSelector().Input   = [
     ### Bd_Dstmunu,Kpi MC10 ### 
    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00008617/0000/00008617_00000032_1.allstreams.dst' TYP='POOL_ROOTTREE' Opt='READ'",
    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00008617/0000/00008617_00000025_1.allstreams.dst' TYP='POOL_ROOTTREE' Opt='READ'",
    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00008617/0000/00008617_00000027_1.allstreams.dst' TYP='POOL_ROOTTREE' Opt='READ'",
    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00008617/0000/00008617_00000042_1.allstreams.dst' TYP='POOL_ROOTTREE' Opt='READ'",
    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00008617/0000/00008617_00000038_1.allstreams.dst' TYP='POOL_ROOTTREE' Opt='READ'",
    "   DATAFILE='root://castorlhcb.cern.ch//castor/cern.ch/grid/lhcb/MC/MC10/ALLSTREAMS.DST/00008617/0000/00008617_00000009_1.allstreams.dst' TYP='POOL_ROOTTREE' Opt='READ'"
    ]

