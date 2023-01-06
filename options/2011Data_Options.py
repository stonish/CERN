#######################################################################
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
from Configurables import Z2TauTau
real = Z2TauTau("Real")
#real.OutputLevel = 2
realseq.Members += [ real ]
Inputs = "StdTightMuons" #"StdLooseMuons" #"StdNoPIDsMuons" #"StdLoosePions" 
#real.InputLocations = [ "StdTightMuons" ] #["StdLooseMuons" ] #[ "StdNoPIDsMuons" ]#[ "StdLoosePions" ]#
#real.InputLocations = [ "Phys/" + Inputs + "/Particles" ]
real.Inputs = [ "Phys/" + Inputs + "/Particles" ]
from GaudiKernel.SystemOfUnits import MeV
#real.MassWindow = 30*MeV
#real.MaxChi2 = 100
real.Particle = "Z0"
#######################################################################
#
# 3) Needed to analyse real data?? -- Doesnt seem to work anynmore?
#
from Configurables import CondDB
#CondDB(UseOracle = True)         # Only use if running on very recent data (taken after the last release of SQLDDDB); then use the online Oracle DB
#CondDB().IgnoreHeartBeat = True  # Ignore the check that the online DB used is valid for the data being processed. Not recomended!
#CondDB().UseLatestTags=["2010"]   # For real data, it should always be safe to use the latest LHCBCOND and DDDB global tags. Not true for MC.

#######################################################################
#
# 4) Configure the application
#
from Configurables import DaVinci
nEvts                   = "_"
type                    = "EWK_2011_R12S17"
mag                     = "_Down" #_Up"
pids                    = "Tight" #"Loose" #"NoPIDs"#"LoosePions" 
date                    = "04Dec11"
temp                    = "/tmp/rhuston/DaVinci/New" 

#DaVinci().HistogramFile = temp + "RealHistos_" + type + mag + "_" + pids + "_" + date + ".root"    # Histogram file
DaVinci().TupleFile     = temp + "Real_" + type + mag + "_" + pids + "_" + date + ".root"     # Tuple file

# Tags
#DaVinci().DDDBtag       = "head-20101026" #"head-20100119" # Should only need to specify for MC. Find out which was used in Feicim
#DaVinci().CondDBtag     = "head-20101112" #"head-20100303" # Should only need to specify for MC. Find out which was used in Feicim
DaVinci().DataType      = "2011"
DaVinci().EvtMax        = -1 #100000
DaVinci().Simulation    = False
DaVinci().Lumi          = True #False 
#importOptions("$APPCONFIGOPTS/DisableLFC.py") # Temporary workaround when using Oracle on the Grid. Should no longer be needed
#MessageSvc().OutputLevel = 5

DaVinci().UserAlgorithms = [ realseq ]
DaVinci().MainOptions  = ""                                                                # None
########################################################################
#
# To run in shell :
# gaudirun.py options/2010Data_Options.py /options/(list_of_files).py
#
########################################################################
