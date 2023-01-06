########################################################################
from Gaudi.Configuration import *
#######################################################################
#
# 1) Let's define a sequence
#
from Configurables import GaudiSequencer
b2dmunuseq = GaudiSequencer("B2DMuNuSeq")
#######################################################################
#
# 2) Add the Tutorial algorithm
#
from Configurables import B2DMuNuTupleAlgorithm#, PhysDesktop
b2dmunu = B2DMuNuTupleAlgorithm("B2DMuNu")
b2dmunuseq.Members += [ b2dmunu ]
#b2dmunu.addTool( PhysDesktop() )
#b2dmunu.PhysDesktop.InputLocations = [ "StdTightMuons" ]
b2dmunu.InputLocations = [ "StdTightMuons" ]
from GaudiKernel.SystemOfUnits import MeV
#z2mumu.MassWindow = 30*MeV
#z2mumu.MaxChi2 = 100
b2dmunu.Particle = "Z0"
#######################################################################
#
# 3) Decay Tree Tuple
#
from Configurables import DecayTreeTuple#, PhysDesktop
tuple = DecayTreeTuple()
#tuple.addTool( PhysDesktop() )
#tuple.PhysDesktop.InputLocations = [ "B2DMuNu" ]
tuple.InputLocations = [ "B2DMuNu" ]
tuple.ToolList += [ "TupleToolTrigger", "TupleToolMCTruth" ]
tuple.Decay = "[B? -> D? ^mu+  ]cc"
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
tree = PrintDecayTree("PrintFoundBs")
tree.InputLocations = [ "B2DMuNu" ]
tree.addTool( PrintDecayTreeTool )
tree.PrintDecayTreeTool.Information = "Name M P Px Py Pz Pt chi2"
#######################################################################
#
# Print  Print All True Bs
#
from Configurables import PrintMCTree, PrintMCDecayTreeTool
mctree = PrintMCTree("PrintTrueBs")
mctree.addTool( PrintMCDecayTreeTool )
mctree.PrintMCDecayTreeTool.Information = "Name M P Px Py Pz Pt chi2"
#mctree.ParticleNames = [
#    "B*0", "B*~0", "B*_2+", "B*_2-", "B*_20", "B*_2~0", "B*_s0", "B*_s~0",
#    "B*_s20", "B*_s2~0", "B*_s00", "B*_s0~0", "B*_00", "B*_0~0", "B*_0+", "B*_0-" ] 
mctree.ParticleNames = [
    "B0",
    "B~0",
    "B+",
    "B-"#,
#    "B*0",
#    "B*~0",
#    "B*-",
#    "B*+",
#    "B*_00",
#    "B*_0~0",
#    "B*_0+",
#    "B*_0-",
#    "B_1(H)0",
#    "B_1(H)~0",
#    "B_1(H)+",
#    "B_1(H)-",
#    "B_1(L)0",
#    "B_1(L)~0",
#    "B_1(L)+",
#    "B_1(L)-",
#    "B*_20",
#    "B*_2~0",
#    "B*_2+",
#    "B*_2-",
#    "B_s0",
#    "B_s~0",
#    "B*_s0",
#    "B*_s~0",
#    "B*_s00",
#    "B*_s0~0",
#    "B_s1(H)0",
#    "B_s1(H)~0",
#    "B_s1(L)0",
#    "B_s1(L)~0",
#    "B*_s20",
#    "B*_s2~0",
#    "B_c+",
#    "B_c-",
#    "B_c*+",
#    "B_c*-",
#    "B_c0*+",
#    "B_c0*-",
#    "B_c1(L)+",
#    "B_c1(L)-",
#    "B_c1(H)+",
#    "B_c1(H)-",
#    "B_c2*+",
#    "B_c2*-"
    ] 
mctree.Depth = 5  # down to the K and mu
#######################################################################
#
# 4) Configure the application
#
from Configurables import DaVinci
nEvts = "30000"
tmp = "/tmp/rhuston/DaVinci/"
DaVinci().HistogramFile = "B2DMuNuHistos_Tight_" + nEvts + "evt_TruthCone.root"    # Histogram file
DaVinci().TupleFile = "B2DMuNu_Tree_Tight_" + nEvts + "evt_TruthCone.root"         # Tuple file
DaVinci().EvtMax = 520000                                                          # Number of events
DaVinci().DataType = "MC09"                                                        # Default is now "MC09", was DC06 before v24
DaVinci().Simulation = True                                                        # It's MC
#
# Add our own stuff
#
DaVinci().UserAlgorithms = [ b2dmunuseq, etuple]#, mctree, tree ]#, tuple]
#DaVinci().HltType = 'Hlt1' #Changed in DaVinci v24r4
DaVinci().Hlt = True
DaVinci().MainOptions  = ""                                         # None
########################################################################
#
# To run in shell :
# gaudirun.py options/B2DMuNuDecayTreeTupleOptions.py /options/B2DlX-AllEvts.py
#
########################################################################
