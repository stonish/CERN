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
from Configurables import Z2TauTauTupleAlgorithm
z2tautau = Z2TauTauTupleAlgorithm("Z2TauTau")
z2tautauseq.Members += [ z2tautau ]
z2tautau.InputLocations = [ "StdTightMuons" ]
from GaudiKernel.SystemOfUnits import MeV
#z2tautau.MassWindow = 30*MeV
#z2tautau.MaxChi2 = 100
z2tautau.Particle = "Z0"
#######################################################################
#
# 3) Decay Tree Tuple
#
from Configurables import DecayTreeTuple
tuple = DecayTreeTuple()
tuple.InputLocations = [ "Z2TauTau" ] 
tuple.ToolList += [ "TupleToolTrigger", "TupleToolMCTruth" ]
tuple.Decay = "[Z0 -> ^mu+ ^mu- ]cc"
#######################################################################
#
# 4) Event Tuple
#
from Configurables import EventTuple
etuple = EventTuple()
etuple.ToolList = [ "TupleToolEventInfo", "TupleToolGeneration", "TupleToolTrigger" ]
#######################################################################
#
# 5) Needed to analyse real data?? (i.e. not MC data)
#
#from Configurables import CondDB
#CondDB(UseOracle = True)
#######################################################################
#
# 6) Configure the application
#
from Configurables import DaVinci
#nEvts                   = "5000"
nEvts                  = "Real"
temp                    = "" #"/tmp/rhuston/DaVinci/"
DaVinci().HistogramFile = temp + "Z2TauTau2MuNuNuHistos_Tight_" + nEvts + "evt.root"#_MC09_Calo.root"    # Histogram file
DaVinci().TupleFile     = temp + "Z2TauTau2MuNuNu_Tree_Tight_" + nEvts + "evt.root"#_MC09_Calo.root"     # Tuple file
DaVinci().EvtMax        = 100 #-1 #121939                                                              # Number of events
DaVinci().DataType      = "MC09"#"2009"#"MC09"                                                           # Default is "MC09" for v24 onwards (was "DC06" before)
DaVinci().Simulation    = True                                                                    # It's MC
#
# Add our own stuff
#
DaVinci().UserAlgorithms = [ z2tautauseq]#, tuple, etuple ]
#DaVinci().HltType = 'Hlt1' # Changed in DaVinci v24r4
#DaVinci().Hlt = True
DaVinci().MainOptions  = ""                                                                # None
########################################################################
#
# To run in shell :
# gaudirun.py options/Z2TauTau2MuNuNuOptions.py /options/Z2TauTau2MuNuNu-AllEvts.py
#
########################################################################
