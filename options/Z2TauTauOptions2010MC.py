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
z2tautau.OutputLevel = 2 # Debug
z2tautauseq.Members += [ z2tautau ]
z2tautau.InputLocations = [ "Phys/StdLooseMuons" ]
from GaudiKernel.SystemOfUnits import MeV
#z2tautau.MassWindow = 30*MeV
#z2tautau.MaxChi2 = 100
z2tautau.Particle = "Z0"
#######################################################################
#
# 3) Configure the application
#
from Configurables import DaVinci
nEvts                   = "All"#"500"#0"
type                    = "Loose" #"Tight"#"NoPIDs"
date                    = "25Nov10"
temp                    = "/tmp/rhuston/DaVinci/" #""
DaVinci().HistogramFile = temp + "Z2TauTauHistos_MC2010_HLT_" + nEvts + "evt_" + type + "_" + date + ".root"   # Histogram file
DaVinci().TupleFile     = temp + "Z2TauTau_Tree_MC2010_HLT_" + nEvts + "evt_" + type + "_" + date + ".root"    # Tuple file
DaVinci().EvtMax        = 500 #-1 #5000                                                                            # Number of events
DaVinci().DataType      = "2010"                                                                           # Default is "MC09"??
DaVinci().Simulation    = True                                                                             # It's MC
DaVinci().Hlt = True      
DaVinci().HltThresholdSettings = "Physics_3000Vis_200L0_20Hlt1_ExpressHlt2_Oct10"#"Physics_320Vis_300L0_10Hlt1_ExpressHlt2_Feb10"#"Physics_320Vis_300L0_10Hlt1_Feb10"   # TEST #some settings. See HltConf for more.
DaVinci().ReplaceL0BanksWithEmulated = True  # to get L0 compatible with Hlt # TEST!
DaVinci().CondDBtag     = "sim-20101026-vc-mu100" # TEST!
#DaVinci().CondDBtag     = "sim-20100715-vc-mu100" # TEST!

# TEST
#DaVinci().Lumi = False # Need to be able to set this to true!
#MessageSvc().OutputLevel = 5
# DDDB / CondDB tags? Disable LFC?

#
# Add our own stuff
#
DaVinci().UserAlgorithms = [ z2tautauseq ]
DaVinci().MainOptions  = ""                    # None
########################################################################
#
# To run in shell :
# gaudirun.py options/Z2TauTauOptions2010MC.py /options/Z2TauTau2MuNuNu-AllEvts.py
#
########################################################################
