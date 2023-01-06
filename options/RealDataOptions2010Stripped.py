########################################################################
#import ROOT
from Gaudi.Configuration import *
#######################################################################
#
# 1) Let's define a sequence
#
from Configurables import GaudiSequencer
realseq = GaudiSequencer("RealSeq")
#######################################################################
#
# 2) Add the Tutorial algorithm
#
from Configurables import Z2TauTauTupleAlgorithm
real = Z2TauTauTupleAlgorithm("Real")
#real.OutputLevel = 2
realseq.Members += [ real ]
real.InputLocations = [ "StdLooseMuons" ] #[ "StdNoPIDsMuons" ]#[ "StdLoosePions" ]#
from GaudiKernel.SystemOfUnits import MeV
#real.MassWindow = 30*MeV
#real.MaxChi2 = 100
real.Particle = "Z0"
#######################################################################
#
# 3) Decay Tree Tuple
#
from Configurables import DecayTreeTuple
tuple = DecayTreeTuple()
tuple.InputLocations = [ "Real" ] 
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
# 5) Needed to analyse real data?? -- Doesnt seem to work anynmore?
#
from Configurables import CondDB
#CondDB(UseOracle = True)
CondDB().IgnoreHeartBeat = True
#######################################################################
#
# 6) Configure the application
#
from Configurables import DaVinci
nEvts                   = "TEST"#"run6234_FakeMuon" #"run5727_M1HitDist"
type                    = "loose" #"NoPIDs"#"LoosePions"#"Tight" 
date                    = "14Sep10"
temp                    = "/tmp/rhuston/DaVinci/" #"PFN:castor:///castor/cern.ch/user/r/rhuston/624/0/" #"/tmp/rhuston/Ganga/624/2/"
DaVinci().HistogramFile = temp + "RealHistos_" + nEvts + "_" + type + "_" + date + ".root"    # Histogram file
DaVinci().TupleFile     = temp + "Real_Tree_" + nEvts + "_" + type + "_" + date + ".root"     # Tuple file
DaVinci().EvtMax        = 1000 #-1 #121939                                                    # Number of events
#DaVinci().DataType      = "2009" #"MC09"                                                       # Default is "MC09" for v24 onwards (was "DC06" before)
#DaVinci().Simulation    = True                                                                   # It's MC
#DaVinci().SkipEvents    = 4000 #20000

# TEST
DaVinci().DDDBtag       = "head-20100119"
DaVinci().CondDBtag     = "head-20100303"
DaVinci().DataType      = "2010"
DaVinci().Simulation    = False
DaVinci().Lumi          = False # Need to be able to set this to true!
importOptions("$APPCONFIGOPTS/DisableLFC.py")
#MessageSvc().OutputLevel = 5
#
# Add our own stuff
#
DaVinci().UserAlgorithms = [ realseq]#, tuple, etuple ]
#DaVinci().HltType = 'Hlt1' # Changed in DaVinci v24r4
#DaVinci().Hlt = True
DaVinci().MainOptions  = ""                                                                # None
########################################################################
#
# To run in shell :
# gaudirun.py options/Z2TauTau2MuNuNuOptions.py /options/Z2TauTau2MuNuNu-AllEvts.py
#
########################################################################
