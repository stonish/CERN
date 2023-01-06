########################################################################
#
# $Id: DaVinci-WriteSelDST.py,v 1.8 2010/01/07 08:18:23 panmanj Exp $
#
# Options for a DaVinci job creating DSTs
#
# It looks for a J/psi and then saves the events in a DST
#
# For a more complicated example with several DSTs
# corresponding to various mass ranges look in
# DaVinciSys/tests/options
#
# @author Patrick Koppenburg
# @date 2009-02-03
#
##############################################################################
from Gaudi.Configuration import *
from Configurables import DaVinci, SelDSTWriter
##############################################################################

# Define a sequence
#
from Configurables import GaudiSequencer
tutorialseq = GaudiSequencer("TutorialSeq")
#######################################################################


#from Configurables import FilterDesktop
MessageSvc().Format = "% F%60W%S%7W%R%T %0W%M"

#
# Get a Z
#
from Configurables import CombineParticles
z2mumu = CombineParticles("Z2MuMu")
z2mumu.InputLocations = [ "Phys/StdLooseMuons" ] 
z2mumu.DecayDescriptor = "Z0 -> mu+ mu-" 
z2mumu.DaughtersCuts = { "mu+" : "PT>0*GeV" } 
#z2mumu.CombinationCut = "AM>2.0*GeV"
z2mumu.MotherCut = "MM>2.5*GeV" 
#z2mumu.MotherCut = "(VFASPF(VCHI2/VDOF)<100)" 
tutorialseq.Members += [ z2mumu ]


#from PhysSelPython.Wrappers import Selection, DataOnDemand, SelectionSequence

#MyLooseJpsi = DataOnDemand('DODStdLooseJpsi2MuMu',
 #                          Location = 'Phys/StdLooseJpsi2MuMu')

#MyLooseZ = DataOnDemand('DODStdLooseMuons',
#                           Location = 'Phys/StdLooseMuons')



from Configurables import InputCopyStream

InputCopyStream().Output = "DATAFILE='PFN:ZCands.dst' TYP='POOL_ROOTTREE' OPT='REC'"
tutorialseq.Members.append(InputCopyStream())



#SelZ = Selection('SelZ',
#                    Algorithm = z2mumu,
#                    RequiredSelections = [MyLooseZ])


#SeqZ = SelectionSequence('SeqJpsi', TopSelection = SelZ)

#from MicroDSTExample.Selections import SeqBs2Jpsi2MuMuPhi2KK
#SeqBs = SeqBs2Jpsi2MuMuPhi2KK.SeqBs2Jpsi2MuMuPhi2KK

#
# This the bit that declares the sequence to the Dst writer
# seq will bee passed to DaVinci(). You do not need to pass it yourself.
# The second line allows to save particles (not yet 100% functional)
#
#dstWriter = SelDSTWriter("ZDSTWriter",
#                         SelectionSequences = [SeqZ],
#                         SaveCandidates = True,
#                       #  CopyMCTruth = True,
#                         WriteFSR = True,
#                         OutputFileSuffix = "Test",
               #          )
##############################################################################
DaVinci().EvtMax = -1
DaVinci().PrintFreq = 10
DaVinci().UserAlgorithms = [ tutorialseq ]
# DaVinci().SkipEvents = 0
DaVinci().DataType = "2010" 
#DaVinci().WriteFSR = True

# DaVinci().Simulation   = False

# input file
#DaVinci().Input =  [




#  "   DATAFILE='PFN:LumiData_062910_0000000001_fsr.dst' TYP='POOL_ROOTTREE' OPT='READ'",
#"DATAFILE='PFN:castor:/castor/cern.ch/grid/lhcb/production/DC06/phys-v4-lumi2/00002146/DST/0000/00002146_00000001_5.dst' TYP='POOL_ROOTTREE' OPT='READ'",
## "DATAFILE='PFN:castor:/castor/cern.ch/grid/lhcb/production/DC06/phys-v4-lumi2/00002146/DST/0000/00002146_00000002_5.dst' TYP='POOL_ROOTTREE' OPT='READ'",
## "DATAFILE='PFN:castor:/castor/cern.ch/grid/lhcb/production/DC06/phys-v4-lumi2/00002146/DST/0000/00002146_00000003_5.dst' TYP='POOL_ROOTTREE' OPT='READ'",
## "DATAFILE='PFN:castor:/castor/cern.ch/grid/lhcb/production/DC06/phys-v4-lumi2/00002146/DST/0000/00002146_00000004_5.dst' TYP='POOL_ROOTTREE' OPT='READ'",
## "DATAFILE='PFN:castor:/castor/cern.ch/grid/lhcb/production/DC06/phys-v4-lumi2/00002146/DST/0000/00002146_00000005_5.dst' TYP='POOL_ROOTTREE' OPT='READ'",
## "DATAFILE='PFN:castor:/castor/cern.ch/grid/lhcb/production/DC06/phys-v4-lumi2/00002146/DST/0000/00002146_00000006_5.dst' TYP='POOL_ROOTTREE' OPT='READ'",
## "DATAFILE='PFN:castor:/castor/cern.ch/grid/lhcb/production/DC06/phys-v4-lumi2/00002146/DST/0000/00002146_00000007_5.dst' TYP='POOL_ROOTTREE' OPT='READ'",
## "DATAFILE='PFN:castor:/castor/cern.ch/grid/lhcb/production/DC06/phys-v4-lumi2/00002146/DST/0000/00002146_00000009_5.dst' TYP='POOL_ROOTTREE' OPT='READ'",
## "DATAFILE='PFN:castor:/castor/cern.ch/grid/lhcb/production/DC06/phys-v4-lumi2/00002146/DST/0000/00002146_00000011_5.dst' TYP='POOL_ROOTTREE' OPT='READ'",
## "DATAFILE='PFN:castor:/castor/cern.ch/grid/lhcb/production/DC06/phys-v4-lumi2/00002146/DST/0000/00002146_00000013_5.dst' TYP='POOL_ROOTTREE' OPT='READ'"
#]

##############################################################################
