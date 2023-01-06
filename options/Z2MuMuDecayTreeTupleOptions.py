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
from Configurables import Z2MuMuTupleAlgorithm#, PhysDesktop
z2mumu = Z2MuMuTupleAlgorithm("Z2MuMu")
z2mumuseq.Members += [ z2mumu ]
#z2mumu.addTool( PhysDesktop() )
#z2mumu.PhysDesktop.InputLocations = [ "StdTightMuons" ]
z2mumu.InputLocations = [ "StdTightMuons" ]
from GaudiKernel.SystemOfUnits import MeV
#z2mumu.MassWindow = 30*MeV
#z2mumu.MaxChi2 = 100
z2mumu.Particle = "Z0"
#######################################################################
#
# 3) Decay Tree Tuple
#
from Configurables import DecayTreeTuple#, PhysDesktop
tuple = DecayTreeTuple()
#tuple.addTool( PhysDesktop() )
#tuple.PhysDesktop.InputLocations = [ "Z2MuMu" ]
tuple.InputLocations = [ "Z2MuMu" ]
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
# Print Reconstructed Bs
#
from Configurables import PrintDecayTree, PrintDecayTreeTool#, PhysDesktop
tree = PrintDecayTree("PrintFoundZs")
tree.InputLocations = [ "Z2MuMu" ]
tree.addTool( PrintDecayTreeTool )
tree.PrintDecayTreeTool.Information = "Name M P Px Py Pz Pt chi2"
#######################################################################
#
# Print  Print All True Bs
#
from Configurables import PrintMCTree, PrintMCDecayTreeTool
mctree = PrintMCTree("PrintTrueZs")
mctree.addTool( PrintMCDecayTreeTool )
mctree.PrintMCDecayTreeTool.Information = "Name M P Px Py Pz Pt chi2"
mctree.ParticleNames = [  "p+" ]
mctree.Depth = 10  # down to the K and mu
#######################################################################
#
# 4) Configure the application
#
from Configurables import DaVinci
nEvts = "50000"
tmp = "/tmp/rhuston/DaVinci/"
DaVinci().HistogramFile = "Z2MuMu_Histos_Tight_" + nEvts + "evt_MC09_TruthCone.root"     # Histogram file
DaVinci().TupleFile = "Z2MuMu_Tree_Tight_" + nEvts + "evt_MC09_TruthCone.root"           # Tuple file
DaVinci().EvtMax = 1000 #50000                                                                 # Number of events
DaVinci().DataType = "MC09"                                                              # Default is "MC09" (was "DC06")
DaVinci().Simulation   = True                                                            # It's MC
#
# Add our own stuff
#
DaVinci().UserAlgorithms = [ z2mumuseq, tuple, etuple]#, mctree]#, tree ]
#DaVinci().HltType = 'Hlt1' #Changed in DaVinci v24r4
DaVinci().Hlt = True
DaVinci().MainOptions  = ""                                                        # None
########################################################################
#
# To run in shell :
# gaudirun.py options/Z2MuMuDecayTreeTupleOptions.py /options/Z2MuMu-AllEvts.py
#
########################################################################
