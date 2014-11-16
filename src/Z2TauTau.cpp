// $Id: Z2TauTau.cpp, 28/09/2009 rhuston Exp $

//#####################//
//### Include files ###//
//#####################// 

//### from Gaudi ###//
#include "GaudiKernel/AlgFactory.h" 
//### local ###//
#include "Z2TauTau.h"
#include "MC.h"
#include "Rec.h"
#include "LoopOnDaughters.h"
//### New Filter ###//
#include "Kernel/ParticleFilters.h"
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
//### MC ###//
#include "Event/MCParticle.h"
//### Determine if MC reconstructible ###//
#include "MCInterfaces/IMCReconstructible.h"
//### Extract L0 Decision ###//
#include "Event/L0DUReport.h"
//### TEST for extrapolator ###//
#include "TrackInterfaces/ITrackExtrapolator.h"
//### Extract Run Number and Event (L0) ID ###//
#include "Event/ODIN.h"

using namespace Gaudi::Units;
using namespace boost::lambda;

//-----------------------------------------------------------------------------
// Implementation file for class : Z2TauTau
//
// 07-07-2009 : Shane Huston
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( Z2TauTau );

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
// Initialise variables defined in header file
Z2TauTau::Z2TauTau( const std::string& name,
                    ISvcLocator* pSvcLocator)
  : DaVinciTupleAlgorithm ( name , pSvcLocator ),
    m_runMC(1),                          // Set to 1 to run over MC data, 0 to skip
    m_passEvent(1),                      // Automatically run code on this event unless later overridden 
    m_errorCode(-1234567),
    m_motherID(0),     
    m_motherMass(0.),
    m_coneMax(1),
    m_numDivs(2), //5
    m_muHits(0),                         // Set to 1 to run over muon hits, calculate distances etc
    m_pi(3.1415926535897932384626434)
{
  declareProperty("Particle", m_motherName = "Undefined" );
  declareProperty("MassWindow", m_motherMassWin = 85.*GeV);
  declareProperty("MaxChi2", m_motherChi2 = 10.);
    
  debug() << "######################" << pSvcLocator << ", " << name << endmsg;
  m_local = pSvcLocator;

}

//=============================================================================
// Destructor
//=============================================================================
Z2TauTau::~Z2TauTau() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode Z2TauTau::initialize() {
  StatusCode sc = DaVinciTupleAlgorithm::initialize(); 
  if ( sc.isFailure() ) return sc;

  if (msgLevel(MSG::DEBUG)) debug() << "==> Initialize" << endmsg;

  const LHCb::ParticleProperty* mother = ppSvc()->find( m_motherName );
  if ( !mother )
  {
    err() << "Cannot find particle property for " << m_motherName << endmsg;
    return StatusCode::FAILURE;
  }
 
  m_motherID = LHCb::ParticleID(mother->pdgID());
  m_motherMass = mother->mass();
  info() << "Will reconstruct " << mother->particle() << " (ID=" << m_motherID.pid()
         << ") with mass " << m_motherMass << " buddy" << endmsg;
  info() << "mass window is " << m_motherMassWin << " Mev" << " buddy" << endmsg;
  info() << "Max chi^2 is " << m_motherChi2 << " buddy" << endmsg;

  //### Define tool to retrieve if MCParticle is reconstructible ###//
  m_recoTool = tool<IMCReconstructible>("MCReconstructible");

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode Z2TauTau::execute() {
  
  if (msgLevel(MSG::DEBUG)) debug() << "==> Execute" << endmsg;
  StatusCode sc = StatusCode::SUCCESS ;
  
  //### Extract Run number and Event number ###//
  const LHCb::ODIN* odin = get<LHCb::ODIN>(LHCb::ODINLocation::Default);
  const unsigned int runNum = odin->runNumber(), evNum = odin->eventNumber();
  debug() << "Run number is " << runNum << " and Event ID is " << evNum 
          << " at time " << odin->eventTime() << endmsg; //odin->gpsTime() << endmsg;
  counter("EventID")++;
  debug() << "### Processing event number " << counter("EventID").nEntries() << " in job ###" << endmsg;
  
  setFilterPassed(true);   // Mandatory. Set to true if event is accepted. 

  //### Only run on MC data if requested ###// 
  if (m_runMC==1)
  {
    //### Test MC data to see if event warrants reconstruction ###// 
    m_passEvent = 0;
    debug() << "Performing Analysis on MC data" << endmsg;
    sc = loopOnMC();
    if (!sc) return sc;
  }
    
  if (m_passEvent==1)
  {
    debug() << "### Event number " << counter("EventID").nEntries() << " accepted for analysis buddy ###" << endmsg;
    counter("EventPasses")++;
    
    debug() << "Extracting daughters..." << endmsg;
    LHCb::Particle::ConstVector daughters = this->i_particles();

    if (daughters.size()!=0)
    {
      debug() << "### There are " << daughters.size() << " possible Muons in this event (EventID "
              << counter("EventID").nEntries() << ") ###" << endmsg;
      counter("MuonsPresent")++;
      
      //### Extract muon hits and muon misid if requested ###//
      if (m_muHits) {
        const std::string& strExtra = "Extra";
        Extra* extra = new Extra(strExtra, m_local);
        Tuple fakeTuple = nTuple("fakeMuonTuple");
        sc = extra->fakeMuon(daughters, fakeTuple);
        Tuple hitTuple = nTuple("muonHitTuple");
        sc = extra->muonChamberHits(hitTuple);
        delete extra;
        if (!sc) return sc;
      }
      
      //### Call loopOnRec ###//
      const std::string& strRec = "Rec";
      Rec* rec = new Rec(strRec, m_local);
      Tuple recAllTuple=nTuple("recAllTuple"), hitDistTuple=nTuple("muHitDistTuple"), motherTuple=nTuple("motherTuple"), 
        candTuple=nTuple("candidateTuple");
      const LHCb::RecVertex::Range prims = this->primaryVertices();
      LHCb::Particle::Vector mothers = 
        rec->loopOnRec(daughters, prims, m_motherID, recAllTuple, hitDistTuple, motherTuple, candTuple, m_runMC);

      delete rec;

      if (!sc) return sc;
    }
    else debug() << "### There are no possible muons in the event!! ###" << endmsg;
    debug() << "### Finishing analysis for Event number " << counter("EventID").nEntries() << " ###" << endmsg;
  }
  else
  {
    debug() << "Event rejected" << endmsg;
    counter("EventFails")++;
  }
  debug() << "Event fully analysed." << endmsg;
  return StatusCode::SUCCESS;
}

//============================================================================= 
// Look at MC data                                                             
//=============================================================================
StatusCode Z2TauTau::loopOnMC()
{
  StatusCode sc = StatusCode::SUCCESS;

  //### Declare vectors used later to store relevant particles, filtered by charge ###//
  LHCb::MCParticle::Vector MuPlus, MuMinus, TauPlus, TauMinus;
 
  //### Obtain MCParticles from the default location and store in MCParticle vector ###//
  m_mcparts.clear();
  LHCb::MCParticles* kmcparts;

  if (exist<LHCb::MCParticles>(LHCb::MCParticleLocation::Default))
  {
    kmcparts = get<LHCb::MCParticles>(LHCb::MCParticleLocation::Default);
    LHCb::MCParticles::iterator j;
    for (j = kmcparts->begin(); j != kmcparts->end(); ++j) m_mcparts.push_back(*j);
  } 
  else debug() << "There are no MC particles in the default location" << endmsg;
  debug() << "There are " << m_mcparts.size() << " MC particles" << endmsg;
  
  if (m_mcparts.size()==0) 
  {
    m_passEvent = 1; // Automatically run event if MC data is unavailable
    m_runMC     = 0; // No need to try looking at MC data any more if it is unavailable
  }
  else 
  {
    //### Declare Tuple to store data ###//
    Tuple mcTuple = nTuple("mcTuple");
    
    //### Test ###//
    int loop = 0; // Is this really needed anymore??
    
    //### Loop over MCParticle vector ###//
    for ( unsigned int k = 0; k<m_mcparts.size(); k++)
    {
      LHCb::MCParticle *part = m_mcparts[k];
      int pid = part->particleID().pid();
      verbose() << "MC Particle is a " << pid << endmsg;
      
      //### Identify the mother of the MCParticle ###//
      const LHCb::MCParticle *mcMother = part->mother();
      
      if (mcMother != NULL)
      {
        //### Filter particles into vector of mu+ and mu- ###//
        if ( abs(pid) == 13 )
        {
          //info() << "MCParticle reconstructible..." << m_recoTool->text(m_recoTool->reconstructible(part)) << endmsg;
          sc = filterMC(part, mcMother, MuMinus, MuPlus, TauMinus, TauPlus);
          if (!sc) return sc;
          loop++;
        }
        //### Repeat for taus ###// 
        if ( abs(pid) == 15 )
        {
          sc = filterMC(part, mcMother, MuMinus, MuPlus, TauMinus, TauPlus);
          if (!sc) return sc;
          loop++;
        }
      }
    }
    
    //### There should be one of each particle from the correct mother per event ###//
    if ( loop == 0 || MuPlus.size() != 1 || MuMinus.size() != 1 || TauPlus.size() != 1 || TauMinus.size() != 1 )
    {
      debug() << "There are too many / few daughters in this event." << endmsg;
      debug() << "There are " << MuPlus.size() << " MC anti-muons" << endmsg;
      debug() << "There are " << MuMinus.size() << " MC muons" << endmsg;
      debug() << "There are " << TauPlus.size() << " MC anti-tauons" << endmsg;
      debug() << "There are " << TauMinus.size() << " MC tauons" << endmsg;
      counter("TooManyMCDaughters")++;
    }
    else
    {
      //### Accept Event ###// 
      m_passEvent = 1;
      
      //### Find DiMuon and DiTauon Invariant Masses, plot along with other quantites ###//
      const std::string& strMC = "mc";
      MC* mc = new MC(strMC, m_local);
      sc = mc->invMassMC(MuPlus, MuMinus, m_mcparts, mcTuple, "Mu");
      if (!sc) return sc;
      sc = mc->invMassMC(TauPlus, TauMinus, m_mcparts, mcTuple, "Tau");
      delete mc;
      if (!sc) return sc;
      mcTuple->write();
    }
  }
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// Filter MC data by particle type and charge
//=============================================================================
StatusCode Z2TauTau::filterMC(LHCb::MCParticle* part, const LHCb::MCParticle* mother,
                              LHCb::MCParticle::Vector& MuMinus, LHCb::MCParticle::Vector& MuPlus, 
                              LHCb::MCParticle::Vector& TauMinus, LHCb::MCParticle::Vector& TauPlus)
{
  StatusCode sc = StatusCode::SUCCESS;
  
  int MotherID = mother->particleID().pid();
  
  //### Mu +/- must come from Tau +/- which comes from Z0 or another Tau +/- ###//
  if ( part->particleID().pid()==13 && MotherID==15)
  {
    const LHCb::MCParticle *gMother = mother->mother();
    int gMotherID = gMother->particleID().pid();
    if (abs(gMotherID)==15 || abs(gMotherID)==23) MuMinus.push_back(part);
  }
  else if ( part->particleID().pid()==-13 && MotherID==-15)
  {
    const LHCb::MCParticle *gMother = mother->mother();
    int gMotherID = gMother->particleID().pid();
    if (abs(gMotherID)==15 || abs(gMotherID)==23) MuPlus.push_back(part); 
  }

  //### Tau +/- must come from Z0 ###//
  if ( part->particleID().pid()==15 && MotherID==23 )         { TauMinus.push_back(part); }
  else if ( part->particleID().pid()==-15 && MotherID==23 )   { TauPlus.push_back(part);  }

  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode Z2TauTau::finalize() {

  if (msgLevel(MSG::DEBUG)) debug() << "==> Finalize" << endmsg;

  return DaVinciTupleAlgorithm::finalize(); 
} 

//=============================================================================
