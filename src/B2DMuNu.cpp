// $Id: B2DMuNu.cpp, 14/08/2009 rhuston Exp $

//#####################//
//### Include files ###//
//#####################// 

//### from Gaudi ###//
#include "GaudiKernel/AlgFactory.h" 
//### local ###//
#include "B2DMuNu.h"
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
using namespace boost::lambda; // For new filter

//-----------------------------------------------------------------------------
// Implementation file for class : B2DMuNu
//
// 14-08-2009 : Shane Huston
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( B2DMuNu );

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
// Initialise variables defined in header file
B2DMuNu::B2DMuNu( const std::string& name,
                  ISvcLocator* pSvcLocator)
  : DaVinciTupleAlgorithm ( name , pSvcLocator ),
    m_runMC(1),                         // Set to 1 to run over MC data, 0 to skip
    m_passEvent(1),                     // Automatically run code on this event unless later overridden
    m_errorCode(-1234567),
    m_motherID(0),     
    m_motherMass(0.), 
    m_coneMax(1),
    m_numDivs(5),
    m_muHits(1),                         // Set to 1 to run over muon hits, calculate distances etc
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
B2DMuNu::~B2DMuNu() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode B2DMuNu::initialize() {
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
         << ") with mass " << m_motherMass << endmsg;
  info() << "mass window is " << m_motherMassWin << " Mev" << endmsg;
  info() << "Max chi^2 is " << m_motherChi2 << endmsg;

  //### Define tool to retrieve if MCParticle is reconstructible ###//
  m_recoTool = tool<IMCReconstructible>("MCReconstructible");
  //### Define function to retrieve MC particle associated with reconstructed ###//
  //m_assoc = new Particle2MCLinker(this, Particle2MCMethod::Links, "");

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode B2DMuNu::execute() {
  
  if (msgLevel(MSG::DEBUG)) debug() << "==> Execute" << endmsg;
  StatusCode sc = StatusCode::SUCCESS;
  
  //### Extract Run number and Event number ###//
  const LHCb::ODIN* odin = get<LHCb::ODIN>(LHCb::ODINLocation::Default);
  const unsigned int runNum = odin->runNumber(), evNum = odin->eventNumber();
  debug() << "Run number is " << runNum << " and Event ID is " << evNum
          << " at time " << odin->eventTime() << endmsg; //odin->gpsTime() << endmsg;
  counter("EventID")++;
  debug() << "### Processing event number " << counter("EventID").nEntries() << " in job ###" << endmsg;
  
  setFilterPassed(true); // Mandatory. Set to true if event is accepted.

  //### Only run on MC data if requested ###//
  if (m_runMC==1)
  {
    //### Test MC data to see if event warrants reconstruction (Require B0/B- OR B0~/B+ both decaying to muons) ###//
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
StatusCode B2DMuNu::loopOnMC()
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
      
      if ( abs(pid) == 13 )
      {
        //### Identify the mother of the MCParticle ###//
        const LHCb::MCParticle *mcMother = part->mother();
                    
      //### TEST ###//
      //if ( abs(pid) == 513 || abs(pid) == 525 || abs(pid) == 515 || abs(pid) == 533 || abs(pid) == 535 || abs(pid) == 10531
      //     || abs(pid) == 10511 || abs(pid) == 10521 || abs(pid) == 20513 || abs(pid) == 20523 || abs(pid) == 20533 ) 
      /*
        if (abs(pid == 511) || abs(pid == 521) )
        {
        //info() << "Found mother B (" << pid << ")" <<  endmsg;
        if (mcMother != NULL)
        {
        info() << "Found B (" << pid << ")" <<  endmsg;
        info() << "B Mother:" << mcMother->particleID().pid() << endmsg;
        const LHCb::MCParticle *gMother = mcMother->mother();
        if ( gMother == NULL ) info() << "B has no GrandMother" << endmsg;
        else info() << "B GrandMother: " << gMother->particleID().pid() << endmsg;
        }
        //else info() << "B has no mother." << endmsg;
        loop++;
        }
      */
      /*
        if (mcMother==NULL)
        {
        info() << "Found particle with NULL mother: " << pid << endmsg;
        loop++;
        }
      */
      
      /*    if (mcMother!=NULL)
            {
            if ( abs(pid) == 13 && (abs(mcMother->particleID().pid()) == 511 || abs(mcMother->particleID().pid()) == 521 ) )
            {
            loop++;
            }
            }*/
      
      
        //### Filter particles into vector of mu+ and mu- ###//
        if ( abs( mcMother->particleID().pid() ) == 511 || abs( mcMother->particleID().pid() ) == 521 )
        {
          //info() << "MCParticle reconstructible..." << m_recoTool->text(m_recoTool->reconstructible(part)) << endmsg;
          sc = filterMC(part, mcMother, MuMinus, MuPlus);
          if (!sc) return sc;
          loop++;
        }
      }
    }
  
    //### There should be one of each particle from the correct mother per event ###//
    if ( loop <= 1 || MuMinus.size() != 1 || MuPlus.size() != 1 ) 
    {
      debug() << "There are too many / few daughters in this event." << endmsg;
      debug() << "There are " << MuMinus.size() << " MC muons" << endmsg;
      debug() << "There are " << MuPlus.size() << " MC anti-muons" << endmsg;
      counter("TooManyMCDaughters")++;
    }
    else
    {
      //### Accep Event ###//
      m_passEvent = 1;
      
      //### Find DiMuon Invariant Masses, plot along with other quantites ###//
      const std::string& strMC = "mc";
      MC* mc = new MC(strMC, m_local);
      sc = mc->invMassMC(MuPlus, MuMinus, m_mcparts, mcTuple, "Mu");
      if (!sc) return sc;
      mcTuple->write();
      delete mc;
    }
  }
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// Filter MC data by particle type and charge
//=============================================================================
StatusCode B2DMuNu::filterMC(LHCb::MCParticle* part, const LHCb::MCParticle* mother,
                             LHCb::MCParticle::Vector& MuMinus, LHCb::MCParticle::Vector& MuPlus)
{
  StatusCode sc = StatusCode::SUCCESS;
  
  //### Mu +/- ###//
  if ( part->particleID().pid() == 13 )          MuMinus.push_back(part);
  else if ( part->particleID().pid() == -13 )    MuPlus.push_back(part); 
  
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode B2DMuNu::finalize() {

  if (msgLevel(MSG::DEBUG)) debug() << "==> Finalize" << endmsg;

  return DaVinciTupleAlgorithm::finalize(); 
} 

//=============================================================================
