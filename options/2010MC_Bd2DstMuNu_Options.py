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
from Configurables import B2DMuNu
sim = B2DMuNu("Sim")
#sim.OutputLevel = 2
simseq.Members += [ sim ]
Inputs = "StdTightMuons" #"StdLooseMuons" #"StdNoPIDsMuons" #"StdLoosePions"
#sim.InputLocations = [ "StdTightMuons" ] #["StdLooseMuons" ] #[ "StdNoPIDsMuons" ]#[ "StdLoosePions" ]#
sim.InputLocations = [ "Phys/" + Inputs + "/Particles" ]
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
#CondDB().IgnoreHeartBeat = True
#######################################################################
#
# 4) Configure the application
#
from Configurables import DaVinci
nEvts                   = "_"
type                    = "MC10"
mag                     = "_Down" #"_Up"
pids                    = "Tight" #"Loose" #"NoPIDs"#"LoosePions" 
date                    = "01Jun11"
temp                    = "" #"/tmp/rhuston/DaVinci/" 

DaVinci().TupleFile     = temp + "Bd2DstMuNu_Tree_" + type + mag + "_" + pids + "_" + date + ".root"     # Tuple file

# Tags
DaVinci().DDDBtag       = "head-20101206" #"head-20100119"
if      mag == "_Down": DaVinci().CondDBtag = "sim-20101210-vc-md100" #"head-20100303"
elif    mag == "_Up"  : DaVinci().CondDBtag = "sim-20101210-vc-mu100"
else                  : print "WARNING: No CondDBtag specified. Code will fail buddy."
DaVinci().DataType      = "2010"
DaVinci().EvtMax        = -1 #121939
DaVinci().Simulation    = True
DaVinci().Lumi          = True # Need to be able to set this to true!
#importOptions("$APPCONFIGOPTS/DisableLFC.py")
#MessageSvc().OutputLevel = 2

DaVinci().UserAlgorithms = [ simseq ]
DaVinci().MainOptions  = ""                                                                # None
########################################################################
#
# To run in shell :
# gaudirun.py options/2010Data_Options.py /options/(list_of_files).py
#
########################################################################
