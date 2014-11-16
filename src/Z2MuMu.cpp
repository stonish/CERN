// $Id: Z2MuMu.cpp, 11/08/2009 rhuston Exp $

//#####################//
//### Include files ###//
//#####################// 

//### from Gaudi ###//
#include "GaudiKernel/AlgFactory.h" 
//### local ###//
#include "Z2MuMu.h"
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
//### Extract HLT decision ###//
//#include "Event/HltDecReports.h"
//### TEST for extrapolator ###//
#include "TrackInterfaces/ITrackExtrapolator.h"
//### Extract Run Number and Event (L0) ID ###//
#include "Event/ODIN.h"

using namespace Gaudi::Units;
using namespace boost::lambda; // For new filter

//-----------------------------------------------------------------------------
// Implementation file for class : Z2MuMu
//
// 11-08-2009 : Shane Huston
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( Z2MuMu );

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
// Initialise variables defined in header file
Z2MuMu::Z2MuMu( const std::string& name,
                ISvcLocator* pSvcLocator)
  //: DVAlgorithm ( name , pSvcLocator ),
  : DaVinciTupleAlgorithm ( name , pSvcLocator ),
  //: DaVinciAlgorithm ( name , pSvcLocator ),
    m_runMC(1),                           // Set to 1 to run over MC Data, 0 to skip
    //m_passEvent(1),                      // Automatically run code on this event unless later overridden
    m_errorCode(-1234567),
    m_motherID(0),     
    m_motherMass(0.),
    m_coneMax(1),
    m_numDivs(5),
    m_muHits(0),                         // Set to 1 to run over muon hits, calculate distances etc
    m_pi(3.1415926535897932384626434)
{
  declareProperty("Particle", m_motherName = "Undefined" );
  declareProperty("MassWindow", m_motherMassWin = 85.*GeV); // was 85
  declareProperty("MaxChi2", m_motherChi2 = 10.); // was 1000

  debug() << "######################" << pSvcLocator << ", " << name << endmsg;
  m_local = pSvcLocator;
}

//=============================================================================
// Destructor
//=============================================================================
Z2MuMu::~Z2MuMu() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode Z2MuMu::initialize() {
  //StatusCode sc = DVAlgorithm::initialize(); 
  StatusCode sc = DaVinciTupleAlgorithm::initialize(); 
  //StatusCode sc = DaVinciAlgorithm::initialize(); 
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
         << ") with mass " << m_motherMass << endmsg;
  info() << "mass window is " << m_motherMassWin << " Mev" << endmsg;
  info() << "Max chi^2 is " << m_motherChi2 << endmsg;

  //### Define tool to retrieve if MCParticle is reconstructible ###//
  m_recoTool = tool<IMCReconstructible>("MCReconstructible");

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode Z2MuMu::execute() {
    
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
    //m_passEvent = 0;
    debug() << "Performing Analysis on MC data" << endmsg;
    sc = loopOnMC();
    if (!sc) return sc;
  }
  
  //if (m_passEvent==1)
  debug() << "### Event number " << counter("EventID").nEntries() << " accepted for analysis buddy ###" << endmsg;
  counter("EventPasses")++;
  
  //### Get Particles ###//
  debug() << "Extracting daughters..." << endmsg;
  //LHCb::Particle::ConstVector daughters = desktop()->particles();
  LHCb::Particle::ConstVector daughters = this->i_particles(); // As of v27r0, PhysDesktop no longer exists!                
  //LHCb::Particle::Range foo = this->particles(); // An alternative to the above?

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

    //### Test - Call loopOnRec ###//
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
  //}
  //else
  //{
  //debug() << "Event rejected" << endmsg;
  //counter("EventFails")++;
  //}

  debug() << "Event fully analysed." << endmsg;
  return StatusCode::SUCCESS;
}

//============================================================================= 
// Look at MC data                                                             
//=============================================================================
StatusCode Z2MuMu::loopOnMC()
{
  StatusCode sc = StatusCode::SUCCESS;

  //### Declare vectors used later to store relevant particles, filtered by charge ###//
  LHCb::MCParticle::Vector MuPlus, MuMinus;
 
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
    //m_passEvent = 1; // Automatically run event if MC data is unavailable
    m_runMC     = 0;  // No need to try looking at MC data any more if it is unavailable
  }
  else
  {
    //### Declare Tuple to store data ###//
    Tuple mcTuple = nTuple("mcTuple");
    
    //### Loop over MCParticle vector ###//
    for ( unsigned int k = 0; k<m_mcparts.size(); k++)
    {
      LHCb::MCParticle *part = m_mcparts[k];
      int pid = part->particleID().pid();
      verbose() << "MC Particle is a " << pid << endmsg;
      
      //### Identify the mother of the MCParticle ###//
      const LHCb::MCParticle *mcMother = part->mother();
      
      //### Filter particles into vector of mu+ and mu- ###//
      //### MC Intermediate Vector Bosons erased in DC06 data so Z0 id cannot be retrieved, use NULL value instead ###//  
      //### If using MC09 Data, need to revert to using Z0 ID ###//
      //### For MC10, just use Z0 id and require mother is not null?? ###//
      if ( (mcMother == NULL || mcMother->particleID().pid() == 23) && abs(pid) == 13 )
      {
        //info() << "MCParticle reconstructible..." << m_recoTool->text(m_recoTool->reconstructible(part)) << endmsg;
        sc = filterMC(part, mcMother, MuMinus, MuPlus);
        if (!sc) return sc;
      }
    }
    
    //### There should be one of each particle from the correct mother per event ###//
    if ( MuPlus.size() != 1 || MuMinus.size() != 1 )
    {
      warning() << "There are too many / few daughters in this event." << endmsg;
      info() << "There are " << MuPlus.size() << " MC anti-muons" << endmsg;
      info() << "There are " << MuMinus.size() << " MC muons" << endmsg;
      counter("TooManyMCDaughters")++;
    }
    //else{
    //m_passEvent = 1;

    //### Find DiMuon Invariant Masses, plot along with other quantites ###//
    const std::string& strMC = "mc";
    MC* mc = new MC(strMC, m_local);
    sc = mc->invMassMC(MuPlus, MuMinus, m_mcparts, mcTuple, "Mu");
    if (!sc) return sc;
    mcTuple->write();
    delete mc;
  }
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// Filter MC data by particle type and charge
//=============================================================================
StatusCode Z2MuMu::filterMC(LHCb::MCParticle* part, const LHCb::MCParticle* mother,
                            LHCb::MCParticle::Vector& MuMinus, LHCb::MCParticle::Vector& MuPlus)
{
  StatusCode sc = StatusCode::SUCCESS;
  
  //### Mu +/- must come from Z0 ###//
  if ( part->particleID().pid()==13 )        MuMinus.push_back(part);
  else if ( part->particleID().pid()==-13 )  MuPlus.push_back(part); 
  
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode Z2MuMu::finalize() {

  if (msgLevel(MSG::DEBUG)) debug() << "==> Finalize" << endmsg;

  //return DVAlgorithm::finalize(); 
  return DaVinciTupleAlgorithm::finalize(); 
  //return DaVinciAlgorithm::finalize(); 
} 

//=============================================================================
