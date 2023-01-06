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
#from Configurables import FilterDesktop
#filter = FilterDesktop("Filter", Code = "(P>30*GeV")
#filter.InputLocations = [ "StdLooseMuons" ]
#filter.Code = "(P>30*GeV)"
#filter.OutputLevel = 4#2 
#realseq.Members += [ filter ]

from Configurables import Z2TauTauTupleAlgorithm
real = Z2TauTauTupleAlgorithm("Real")
real.OutputLevel = 4#2
realseq.Members += [ real ]
real.InputLocations = [ "StdLooseMuons" ]#[ "StdNoPIDsMuons" ]#[ "StdLoosePions" ]
#real.InputLocations = [ "Filter" ]
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
tuple.OutputLevel = 2
tuple.InputLocations = [ "Test" ] 
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
# 5) Needed to analyse real data -- Only true for "recent data"
#
#from Configurables import CondDB
#CondDB(UseOracle = True) #Can comment out for interactive but works with --use-grid (Only true for some of the data) 
#CondDB().IgnoreHeartBeat = True
#######################################################################
#
# 6) Configure the application
#
from Configurables import DaVinci
nEvts                   = "run7019_MuonEffAndRes_DiMuonStrip05" #"run5727_M1HitDist"
type                    = "Loose"#_P30GeV"#"Loose"#"LoosePions"#"Tight"#"NoPIDs" 
date                    = "03Aug10"
temp                    = ""#"/afs/cern.ch/user/r/rhuston/scratch0/DaVinci/"#"/tmp/rhuston/DaVinci/"#"PFN:castor:///castor/cern.ch/user/r/rhuston/624/0/"
DaVinci().HistogramFile = temp + "RealHistos_" + nEvts + "_" + type + "_" + date + ".root"    # Histogram file
DaVinci().TupleFile     = temp + "Real_Tree_" + nEvts + "_" + type + "_" + date + ".root"     # Tuple file
DaVinci().EvtMax        = -1#110#00#121939                                                             # Number of events
#DaVinci().DataType      = "2009" #"MC09"                                                         # Default is "MC09" for v24 onwards (was "DC06" before)
#DaVinci().Simulation    = True                                                                   # It's MC
#DaVinci().SkipEvents    = 4000 #20000

# TEST
DaVinci().DDDBtag       = "head-20100119"
DaVinci().CondDBtag     = "head-20100303"
DaVinci().DataType      = "2010"
DaVinci().Simulation    = False
DaVinci().Lumi          = True # TEST : Haven't tried this before.. I think it should calculate the luminosity of the data analysed
#importOptions("$APPCONFIGOPTS/DisableLFC.py") #commented out for interactive run or when UseOracle==False
#MessageSvc().OutputLevel = 5
#
# Add our own stuff
#
DaVinci().UserAlgorithms = [ realseq]#, tuple]#, etuple ]
#DaVinci().HltType = 'Hlt1' # Changed in DaVinci v24r4
#DaVinci().Hlt = True
DaVinci().MainOptions  = ""                                                                # None
########################################################################
#
# To run in shell :
# gaudirun.py options/Z2TauTau2MuNuNuOptions.py /options/Z2TauTau2MuNuNu-AllEvts.py
#
########################################################################
