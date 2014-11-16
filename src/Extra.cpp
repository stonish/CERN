// $Id: LoopOnDaughters.cpp, 12/04/2011 rhuston Exp $

//#####################//
//### Include files ###//
//#####################// 

//### from Gaudi ###//
#include "GaudiKernel/AlgFactory.h" 
//### local ###//
#include "Extra.h"
//### Work with Muon Detector ###//
#include "MuonDet/DeMuonDetector.h"

using namespace Gaudi::Units;

//-----------------------------------------------------------------------------
// Implementation file for class : Extra
//
// 12-04-2011 : Shane Huston
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( Extra );

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
Extra::Extra( const std::string& name, ISvcLocator* pSvcLocator)
  : DaVinciTupleAlgorithm ( name , pSvcLocator ),
    m_errorCode(-1234567)
{
  m_local = pSvcLocator;
  //m_distTool = tool<IDistanceCalculator>("LoKi::DistanceCalculator");
}

//=============================================================================
// Destructor
//=============================================================================
Extra::~Extra() { } 

//=============================================================================
// Initialization
//=============================================================================
StatusCode Extra::initialize() {
  StatusCode sc = DaVinciTupleAlgorithm::initialize(); 
  if ( sc.isFailure() ) return sc;
  if (msgLevel(MSG::DEBUG)) debug() << "==> Initialize" << endmsg;
  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode Extra::execute() {
  if (msgLevel(MSG::DEBUG)) debug() << "==> Execute" << endmsg;
  StatusCode sc = StatusCode::SUCCESS ;
  setFilterPassed(true);   // Mandatory. Set to true if event is accepted. 
  err() << "Execute phase of Extra.cpp has been called." << endmsg; // Execute phase should not be called 
  return StatusCode::SUCCESS;
}

//=============================================================================
// Calculate total IP for two particles
//=============================================================================
StatusCode Extra::ImpactParameterSum(const LHCb::Particle* particle1, const LHCb::Particle* particle2,
                                     const LHCb::VertexBase* vertex,
                                     double& ip1, double& ip2, double &ipe1, double& ipe2, double& ipTot, double& ipeTot)
{
  StatusCode sc = StatusCode::SUCCESS;

  sc = SignedImpactParameter(particle1, vertex, ip1, ipe1);
  if (!sc) 
  {
    ip1 = m_errorCode*mm;
    ipe1 = m_errorCode*mm;
    ipTot = m_errorCode*mm;
    ipeTot = m_errorCode*mm;
    debug() << "SHUST: Error calculating SignedImpactParameter for first particle. Setting all values to " << ip1 << endmsg;
    return sc;
  }

  sc = SignedImpactParameter(particle2, vertex, ip2, ipe2);
  if (!sc) 
  {
    debug() << "SHUST: Error calculating SignedImpactParameter for second particle. Setting all values to " << m_errorCode*mm 
           << endmsg;
    ip2 = m_errorCode*mm;
    ipe2 = m_errorCode*mm;
    ipTot = m_errorCode*mm;
    ipeTot = m_errorCode*mm;
    return sc;
  }

  ipTot = ip1 + ip2;
  ipeTot = sqrt(pow(ipe1, 2) + pow(ipe2, 2));

  //SHUST: What if ipe1 or ipe2 are less than zero? Previously, set this to error code
  if (ipe1 == 0)
  {
    debug() << "SHUST: Error calculating IPE for first particle" << endmsg;
    ipe1 = m_errorCode*mm;
    ipeTot = m_errorCode*mm;
  }
  if (ipe2 == 0) 
  {
    debug() << "SHUST: Error calculating IPE for second particle" << endmsg;
    ipe2 = m_errorCode*mm;
    ipeTot = m_errorCode*mm;
  }

  return sc;
}

//=============================================================================
// Given a particle and a vertex, calculate ImpactParameter and determine sign
//=============================================================================
StatusCode Extra::SignedImpactParameter(const LHCb::Particle* particle, const LHCb::VertexBase* vertex, double& ip,
                                        double& ipError)
{
  StatusCode sc;
  Gaudi::XYZVector ipVector = Gaudi::XYZVector(0, 0, 0);

  sc = m_distTool->distance(particle, vertex, ipVector);
  if (!sc) return sc;

  // Determine sign here

  sc = m_distTool->distance(particle, vertex, ip, ipError);

  return sc;
}

//=============================================================================
// Given a particle and a point, calculate ImpactParameter and determine sign
//=============================================================================
StatusCode Extra::SignedImpactParameter(const LHCb::Particle* particle, const Gaudi::XYZPoint& point, double& ip,
                                        double& ipError)
{
  StatusCode sc;
  Gaudi::XYZVector ipVector = Gaudi::XYZVector(0, 0, 0);

  sc = m_distTool->distance(particle, point, ipVector);
  if (!sc) return sc;

  // Determine sign here

  sc = m_distTool->distance(particle, point, ip, ipError);

  return sc;
}

//=============================================================================
// TEST - Find probability of pion punch thru / decay in flight
//=============================================================================
StatusCode Extra::fakeMuon(const LHCb::Particle::ConstVector& daughters, Tuple fakeTuple)
{
  StatusCode sc = StatusCode::SUCCESS ;
  
  //### For 2009 data, assume all muon candidates are "fake" muons                ###//  
  //### Need to work out teh probability of having a fake muon from the 2009 data ###//

  //### Pre-loop to establish total numbers of Muon and MuonLoose candidates ###//
  int numMu=0, numMuLoose=0;
  for(LHCb::Particle::ConstVector::const_iterator ida = daughters.begin(); ida!=daughters.end(); ++ida)
  {
    const LHCb::Particle* da = *ida;
    if(da->proto()->muonPID()->IsMuonLoose()) numMuLoose++;
    if(da->proto()->muonPID()->IsMuon())      numMu++;
  }

  debug() << "In this event we have " << daughters.size() << " particles, " << numMu
          << " muon candidates and " << numMuLoose << " loose muon candidates." << endmsg;
  
  //### Loop on all daughters and plot some details ###//
  int da_i=0, mu_i=0, muLoose_i=0;
  int test1=0, test2=0;
  for(LHCb::Particle::ConstVector::const_iterator ida = daughters.begin(); ida!=daughters.end(); ++ida)
  {
    fakeTuple->column("EventID",          counter("EventID").nEntries());
    fakeTuple->column("NumLongTracks",    daughters.size());
    fakeTuple->column("NumMuonCand",      numMu);
    fakeTuple->column("NumMuonLooseCand", numMuLoose);
    
    const LHCb::Particle* da = *ida;
    da_i++;
    fakeTuple->column("LoopID", da_i);
    
    const LHCb::MuonPID* mu = da->proto()->muonPID();
    if(mu->IsMuonLoose())
    {
      muLoose_i++;
      debug() << "Looking at loose muon candidate #" << muLoose_i << endmsg;
      fakeTuple->column("MuLooseID", muLoose_i);
    }
    else fakeTuple->column("MuLooseID", 0);
    
    if(mu->IsMuon())
    {
      mu_i++;
      debug() << "Looking at muon candidate #" << mu_i << endmsg;
      fakeTuple->column("MuID", mu_i);
    }
    else fakeTuple->column("MuID", 0);
    
    if(da_i==daughters.size()) fakeTuple->column("IsFinalLoop",    1);
    else                       fakeTuple->column("IsFinalLoop",    0);
    
    if(numMuLoose==0)
    {
      if (da_i==daughters.size())
      {
        debug() << "This is the Final Loose Muon Candidate" << endmsg;
        fakeTuple->column("IsFinalMuLoose", 1);
      }
    }
    else
    {
      if(muLoose_i==numMuLoose && test1==0)
      {
        debug() << "This is the final Loose Muon Candidate" << endmsg;
        test1++;
        fakeTuple->column("IsFinalMuLoose", 1);
      }
      else                       fakeTuple->column("IsFinalMuLoose", 0);
    }
    
    if(numMu==0)
    {
      if (da_i==daughters.size())
      {
        debug() << "This is the final Muon Candidate" << endmsg;
        if (da_i==daughters.size()) fakeTuple->column("IsFinalMu", 1);
      }
    }
    else
    {
      if(mu_i==numMu && test2==0)
      {
        debug() << "This is the final Muon Candidate" << endmsg;
        test2++;
        fakeTuple->column("IsFinalMu",      1);
      }
      else                       fakeTuple->column("IsFinalMu",      0);
    }
    
    fakeTuple->column("IsMuon",      da->proto()->muonPID()->IsMuon());
    fakeTuple->column("IsMuonLoose", da->proto()->muonPID()->IsMuonLoose());
    fakeTuple->column("P",           da->p());
    fakeTuple->column("Pt",          da->pt());
    fakeTuple->column("Eta",         da->momentum().Eta());
    fakeTuple->column("Phi",         da->momentum().phi());
    
    fakeTuple->write();
  }
  return StatusCode::SUCCESS;
}

//=============================================================================
// TEST - Store all available info for muon chamber hits                                                                  
//=============================================================================
StatusCode Extra::muonChamberHits(Tuple hitTuple)
{
  StatusCode sc = StatusCode::SUCCESS ;
  
  //### TEST - extract all muonPIDs, where possible, find muon hits and store coordinates ##//
  //Tuple hitTuple = nTuple("muonHitTuple");
  std::vector<Gaudi::XYZPoint> allMuonHits;
  //Vector to store coordinates of all hits retrieved
  LHCb::MuonPIDs* AllMuPIDs = get<LHCb::MuonPIDs>(LHCb::MuonPIDLocation::Default);
  
  for (LHCb::MuonPIDs::iterator m = AllMuPIDs->begin(); m!=AllMuPIDs->end(); ++m)
  {
    LHCb::MuonPID* mupid = *m;
    if (mupid->IsMuonLoose())
    {
      debug() << "MuonPID IsMuonLoose() = " << mupid->IsMuonLoose() << endmsg;
      const LHCb::Track* muT = mupid->muonTrack();
      const LHCb::Track* longT = mupid->idTrack();
      debug() << "Have retrieved MuonTrack with momentum " << muT->p() << endmsg;
      
      LHCb::Track::LHCbIDContainer ids = muT->lhcbIDs();
      for (LHCb::Track::LHCbIDContainer::const_iterator itIDs = ids.begin(); itIDs!=ids.end(); ++itIDs)
      {
        if (itIDs->isMuon()) // Just to double check
        {
          debug() << "LHCbIDs for this muon track is of muon type." << endmsg;
          LHCb::MuonTileID tileID = itIDs->muonID();
          debug()  << "Muon passes through muon station " << tileID.station()+1 << endmsg;
          debug() << "*** tile position ***" << tileID << endmsg;
          
          //### Access DeMuonDetector class to convert muonTileID to XYZ coordinates ###//
          double x, y, z, dx, dy, dz;
          DeMuonDetector* muonDet;
          muonDet = getDet<DeMuonDetector>(DeMuonLocation::Default);
          muonDet->Tile2XYZ( tileID, x, dx, y, dy, z, dz );
          
          debug() << "The tile's position is (" << x << ", " << y << ", " << z << ")" << endmsg;
          debug() << "and deltas (" << dx << ", " << dy << ", " << dz << ")" << endmsg;
          
          hitTuple->column("Station",   tileID.station()+1);
          hitTuple->column("Region",    tileID.region());
          hitTuple->column("Quarter",   tileID.quarter());
          hitTuple->column("XIndex",    tileID.nX());
          hitTuple->column("YIndex",    tileID.nY());
          hitTuple->column("isValid",   tileID.isValid());
          hitTuple->column("isDefined", tileID.isDefined());
          hitTuple->column("X",         x);
          hitTuple->column("Y",         y);
          hitTuple->column("Z",         z);
          hitTuple->column("dX",        dx);
          hitTuple->column("dY",        dy);
          hitTuple->column("dZ",        dz);
          hitTuple->column("isGaudi",   0);
          
          hitTuple->write();
          
          //### Store vector of XYZ coordinates for each muonTileID ###//
          Gaudi::XYZPoint test = Gaudi::XYZPoint(x,y,z);
          allMuonHits.push_back(test);
        }
        else debug() << "LHCbID for this track is not a muon!" << endmsg;
      }
    }
    else debug() << "Particle does not satisfy looseMuon criteria." << endmsg;
  }
  
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode Extra::finalize() {

  if (msgLevel(MSG::DEBUG)) debug() << "==> Finalize" << endmsg;

  return DaVinciTupleAlgorithm::finalize(); 
} 
//=============================================================================
