#######################################################################
#import ROOT
from Gaudi.Configuration import *
#######################################################################
#
# 1) Let's define a sequence
#
from Configurables import GaudiSequencer
simseq = GaudiSequencer("SimSeq")
#######################################################################
#
# 2) Add the Tutorial algorithm
#
from Configurables import Z2TauTauTupleAlgorithm
sim = Z2TauTauTupleAlgorithm("Sim")
#sim.OutputLevel = 2
simseq.Members += [ sim ]
sim.InputLocations = [ "StdTightMuons" ] #["StdLooseMuons" ] #[ "StdNoPIDsMuons" ]#[ "StdLoosePions" ]#
from GaudiKernel.SystemOfUnits import MeV
#sim.MassWindow = 30*MeV
#sim.MaxChi2 = 100
sim.Particle = "Z0"
#######################################################################
#
# 3) Needed to analyse real data?? -- Doesnt seem to work anynmore?
#
from Configurables import CondDB
#CondDB(UseOracle = True)
CondDB().IgnoreHeartBeat = True
#######################################################################
#
# 4) Configure the application
#
from Configurables import DaVinci
nEvts                   = "_"
type                    = "MC10"
mag                     = "_Down" #"_Up"
pids                    = "Tight" #"Loose" #"NoPIDs"#"LoosePions" 
date                    = "15Apr11"
temp                    = "/tmp/rhuston/DaVinci/" 

#DaVinci().HistogramFile = temp + "Z2TauTauHistos_" + type + mag + "_" + pids + "_" + date + ".root"    # Histogram file.. trying to remove
DaVinci().TupleFile     = temp + "Z2TauTau_Tree_" + type + mag + "_" + pids + "_OLD_" + date + ".root"     # Tuple file

# Tags
DaVinci().DDDBtag       = "head-20100119"
DaVinci().CondDBtag     = "head-20100303"
DaVinci().DataType      = "2010"
DaVinci().EvtMax        = -1 #121939
DaVinci().Simulation    = False
DaVinci().Lumi          = False # Need to be able to set this to true!
importOptions("$APPCONFIGOPTS/DisableLFC.py")
#MessageSvc().OutputLevel = 2

DaVinci().UserAlgorithms = [ simseq ]
DaVinci().MainOptions  = ""                                                                # None
########################################################################
#
# To run in shell :
# gaudirun.py options/2010Data_Options.py /options/(list_of_files).py
#
########################################################################
