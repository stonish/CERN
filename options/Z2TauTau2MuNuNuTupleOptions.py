########################################################################
from Gaudi.Configuration import *
#######################################################################
#
# 1) Let's define a sequence
#
from Configurables import GaudiSequencer
z2tautauseq = GaudiSequencer("Z2TauTauSeq")
#######################################################################
#
# 2) Add the Tutorial algorithm
#
from Configurables import Z2TauTauTupleAlgorithm_Err, PhysDesktop
z2tautau = Z2TauTauTupleAlgorithm_Err("Z2TauTau")
z2tautauseq.Members += [ z2tautau ]
z2tautau.addTool( PhysDesktop() )
z2tautau.PhysDesktop.InputLocations = [ "StdLooseMuons" ]
from GaudiKernel.SystemOfUnits import MeV
#z2tautau.MassWindow = 30*MeV
#z2tautau.MaxChi2 = 100
z2tautau.Particle = "Z0"
#######################################################################
#
# 3) Configure the application
#
from Configurables import DaVinci
nEvts = "5000"
DaVinci().HistogramFile = "Z2TauTau2MuNuNuHistos_Err_" + nEvts + "evt.root"     # Histogram file
DaVinci().TupleFile = "Z2TauTau2MuNuNu_Tree_Err_" + nEvts + "evt.root"          # Tuple file
DaVinci().EvtMax = 5000                                                         # Number of events
DaVinci().DataType = "2008"                                                     # Default is "DC06"
DaVinci().Simulation   = True                                                   # It's MC
#
# Add our own stuff
#
DaVinci().UserAlgorithms = [ z2tautauseq ]
DaVinci().MainOptions  = ""                    # None
########################################################################
#
# To run in shell :
# gaudirun.py options/Z2TauTau2MuNuNuOptions.py /options/Z2TauTau2MuNuNu-AllEvts.py
#
########################################################################
