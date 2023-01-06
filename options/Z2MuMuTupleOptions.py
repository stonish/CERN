########################################################################
from Gaudi.Configuration import *
#######################################################################
#
# 1) Let's define a sequence
#
from Configurables import GaudiSequencer
z2mumuseq = GaudiSequencer("Z2MuMuSeq")
#######################################################################
#
# 2) Add the Tutorial algorithm
#
from Configurables import Z2MuMuTupleAlgorithm, PhysDesktop
z2mumu = Z2MuMuTupleAlgorithm("Z2MuMu")
z2mumuseq.Members += [ z2mumu ]
z2mumu.addTool( PhysDesktop() )
z2mumu.PhysDesktop.InputLocations = [ "StdLooseMuons" ]
from GaudiKernel.SystemOfUnits import MeV
#z2tautau.MassWindow = 30*MeV
#z2tautau.MaxChi2 = 100
z2mumu.Particle = "Z0"
#######################################################################
#
# 3) Configure the application
#
from Configurables import DaVinci
nEvts = "5000"
DaVinci().HistogramFile = "Z2MuMuHistos_" + nEvts + "evt.root"     # Histogram file
DaVinci().TupleFile = "Z2MuMu_Tree_" + nEvts + "evt.root"          # Tuple file
DaVinci().EvtMax = 5000                                            # Number of events
DaVinci().DataType = "DC06"                                        # Default is "DC06"
DaVinci().Simulation   = True                                      # It's MC
#
# Add our own stuff
#
DaVinci().UserAlgorithms = [ z2mumuseq ]
DaVinci().MainOptions  = ""                    # None
########################################################################
#
# To run in shell :
# gaudirun.py options/Z2MuMuTupleOptions.py /options/Z2TauTau2MuNuNu-AllEvts.py
#
########################################################################
