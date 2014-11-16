// $Id: LoopOnDaughters.cpp, 30/09/2009 rhuston Exp $

//#####################//
//### Include files ###//
//#####################// 

//#include "Kernel/DVAlgorithm.h"
//### from Gaudi ###//
#include "GaudiKernel/AlgFactory.h" 
//### local ###//
#include "LoopOnDaughters.h"
//### MC ###//
#include "Event/MCParticle.h"
//### Determine if MC reconstructible ###//
#include "MCInterfaces/IMCReconstructible.h"
//### Extract L0 Decision ###//
#include "Event/L0DUReport.h"
//### Extract HLT decision ###//
//#include "Event/HltDecReports.h"
#include "HltDecReports.h"
//### Extrapolate tracks ###//
#include "TrackInterfaces/ITrackExtrapolator.h"
//### Work with Muon Detector ###//
#include "MuonDet/DeMuonDetector.h"
//### TEST ###//
#include "Event/MuonCoord.h"
//### Extract Run Number and Event (L0) ID ###//
#include "Event/ODIN.h"
//### TEST trying to access Global Event Cut variables ###//
#include "Event/VeloCluster.h"
#include "Event/STMeasurement.h"
#include "Event/OTTime.h"

#include "IsolationMeasurement.h"
 
using namespace Gaudi::Units;

//-----------------------------------------------------------------------------
// Implementation file for class : LoopOndaughters
//
// 30-09-2009 : Shane Huston
//-----------------------------------------------------------------------------

int LoopOnDaughters::constructorCalls;
int LoopOnDaughters::destructorCalls;

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( LoopOnDaughters );

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
// Initialise variables defined in header file
LoopOnDaughters::LoopOnDaughters( const std::string& name, ISvcLocator* pSvcLocator)
  //: DVAlgorithm ( name , pSvcLocator ),
  : DaVinciTupleAlgorithm ( name , pSvcLocator ),
  //: DaVinciAlgorithm ( name , pSvcLocator ),
    m_errorCode(-1234567),
    m_coneMax(1),
    m_numDivs(2), //(5)
    m_muHits(0),           // Set to 1 to run on muon hits, determine distances etc
    m_pi(3.1415926535897932384626434)
  
{
  debug() << "SHUST: LoopOnDaughters constructor called" << endmsg;
  m_local = pSvcLocator;
  //const std::string& strExtra = "Extra";
  //m_extra = new Extra(strExtra, m_local);
  LoopOnDaughters::constructorCalls++;
  if (LoopOnDaughters::constructorCalls - LoopOnDaughters::destructorCalls != 1)
    info() << "SHUST: LoopOnDaughters Constructor misbalance: " << LoopOnDaughters::constructorCalls 
           << " vs " << LoopOnDaughters::destructorCalls << endmsg;
  
}

//=============================================================================
// Destructor
//=============================================================================
LoopOnDaughters::~LoopOnDaughters() 
{
  LoopOnDaughters::destructorCalls++;
  if (LoopOnDaughters::constructorCalls - LoopOnDaughters::destructorCalls != 0)
    info() << "SHUST: LoopOnDaughters destructor misbalance: " << LoopOnDaughters::constructorCalls 
           << " vs " << LoopOnDaughters::destructorCalls << endmsg;
  
  //delete m_extra;
} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode LoopOnDaughters::initialize() {
  //StatusCode sc = DVAlgorithm::initialize(); 
  StatusCode sc = DaVinciTupleAlgorithm::initialize(); 
  //StatusCode sc = DaVinciAlgorithm::initialize(); 
  if ( sc.isFailure() ) return sc;

  if (msgLevel(MSG::DEBUG)) debug() << "==> Initialize" << endmsg;

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode LoopOnDaughters::execute() {

  if (msgLevel(MSG::DEBUG)) debug() << "==> Execute" << endmsg;
  StatusCode sc = StatusCode::SUCCESS ;

  setFilterPassed(true);   // Mandatory. Set to true if event is accepted. 

  // code goes here
  err() << "Execute phase of LoopOnDaughters.cpp has been called." << endmsg;
  //sc = loopOnDaughters(daughters);
  //if (!sc) return sc;
 
  return StatusCode::SUCCESS;
}

//=============================================================================
// loop on daughters
//=============================================================================
StatusCode LoopOnDaughters::loopOnDaughters(const LHCb::Particle::ConstVector& daughters, 
                                            const LHCb::RecVertex::Range prims, Tuple recAllTuple, 
                                            Tuple hitDistTuple, int runMC)const 
{
  StatusCode sc = StatusCode::SUCCESS ;

  int daNum = 0;
  
  for ( LHCb::Particle::ConstVector::const_iterator im =  daughters.begin() ;
        im != daughters.end() ; ++im )
  {
    daNum++;
    debug() << "### Looking at daughter #" << daNum << " ###" << endmsg;
    //### Plot data for reconstructed and associated MC particle ###//
    sc = plotDaughter(*im, prims, recAllTuple, hitDistTuple, 
                      "All", runMC);
    if (!sc) return sc;
    recAllTuple->write();
  }

  return StatusCode::SUCCESS;
}


//=============================================================================
// Plot reconstructed data
//=============================================================================
StatusCode LoopOnDaughters::plotReconstructedData(const LHCb::Particle* da, Tuple daTuple, const std::string& head) const
{
  daTuple->column("rec" + head + "ID",              da->particleID().pid());
  daTuple->column("rec" + head + "Mass",            da->momentum().M());
  daTuple->column("rec" + head + "P",               da->p());
  daTuple->column("rec" + head + "Pt",              da->pt());
  daTuple->column("rec" + head + "E",               da->momentum().E());
  daTuple->column("rec" + head + "Phi",             da->momentum().phi());
  daTuple->column("rec" + head + "Eta",             da->momentum().Eta());
  daTuple->column("rec" + head + "sigmaP",          sqrt(da->proto()->track()->firstState().errP2()));
  daTuple->column("rec" + head + "IsMuon",          da->proto()->muonPID()->IsMuon());
  daTuple->column("rec" + head + "IsMuonLoose",     da->proto()->muonPID()->IsMuonLoose());
  daTuple->column("rec" + head + "TrackChi2",       da->proto()->track()->chi2());
  daTuple->column("rec" + head + "TrackNumDoF",     da->proto()->track()->nDoF());
  daTuple->column("rec" + head + "TrackChi2PerDoF", da->proto()->track()->chi2PerDoF());
  daTuple->column("rec" + head + "TrackProbChi2",   da->proto()->track()->probChi2());
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// Plot LHCbID details
//=============================================================================
StatusCode LoopOnDaughters::plotLHCbIDs(LHCb::Track::LHCbIDContainer lhcbIDs, Tuple daTuple, const std::string& head) const
{
  //### Extract details on on LHCbIDs ###//
  //LHCb::Track::LHCbIDContainer lhcbIDs = da->proto()->track()->lhcbIDs();
  int numVeloIDs=0, numTTIDs=0, numITIDs=0, numOTIDs=0, numRichIDs=0, numCaloIDs=0;
  
  for ( LHCb::Track::LHCbIDContainer::const_iterator ihit = lhcbIDs.begin(); ihit != lhcbIDs.end(); ++ihit )
  {
    if      ( ihit->isVelo() ) numVeloIDs++;
    else if ( ihit->isTT() )   numTTIDs++;
    else if ( ihit->isIT() )   numITIDs++;
    else if ( ihit->isOT() )   numOTIDs++;
    else if ( ihit->isRich() ) numRichIDs++;
    else if ( ihit->isCalo() ) numCaloIDs++;
  }
  
  daTuple->column("rec" + head + "TrackNumLHCbIDs", (unsigned long long)lhcbIDs.size());
  daTuple->column("rec" + head + "TrackNumVeloIDs", numVeloIDs);
  daTuple->column("rec" + head + "TrackNumTTIDs",   numTTIDs);
  daTuple->column("rec" + head + "TrackNumITIDs",   numITIDs);
  daTuple->column("rec" + head + "TrackNumOTIDs",   numOTIDs);
  daTuple->column("rec" + head + "TrackNumRichIDs", numRichIDs);
  daTuple->column("rec" + head + "TrackNumCaloIDs", numCaloIDs);

  return StatusCode::SUCCESS;
}

//=============================================================================
// Plot Global Event Cut Variables
//=============================================================================
StatusCode LoopOnDaughters::plotGlobalEventCutVariables(Tuple daTuple, const std::string& head) const
{
  //### Global Event Cut variables ###//
  LHCb::VeloClusters* allVeloClusters = get<LHCb::VeloClusters>(LHCb::VeloClusterLocation::Default);
  LHCb::STClusters*   allTTClusters   = get<LHCb::STClusters>(LHCb::STClusterLocation::TTClusters);
  LHCb::STClusters*   allITClusters   = get<LHCb::STClusters>(LHCb::STClusterLocation::ITClusters);
  LHCb::OTTimes*      allOTTimes      = get<LHCb::OTTimes>(LHCb::OTTimeLocation::Default);
  LHCb::CaloDigits*   allSpdDigits    = get<LHCb::CaloDigits>(LHCb::CaloDigitLocation::Spd);
  LHCb::CaloDigits*   allPrsDigits    = get<LHCb::CaloDigits>(LHCb::CaloDigitLocation::Prs);
  LHCb::CaloDigits*   allEcalDigits   = get<LHCb::CaloDigits>(LHCb::CaloDigitLocation::Ecal);
  LHCb::CaloDigits*   allHcalDigits   = get<LHCb::CaloDigits>(LHCb::CaloDigitLocation::Hcal);
  
  int numVeloTracks = 0;
  if ( exist<LHCb::Tracks>(LHCb::TrackLocation::Velo) ) {
    LHCb::Tracks*       allVeloTracks   = get<LHCb::Tracks>(LHCb::TrackLocation::Velo);
    numVeloTracks = allVeloTracks->size();
  }
  
  daTuple->column("rec" + head + "NumVeloClusters", (unsigned long long)allVeloClusters->size());
  daTuple->column("rec" + head + "NumOTClusters",   (unsigned long long)allOTTimes->size());
  daTuple->column("rec" + head + "NumITClusters",   (unsigned long long)allITClusters->size());
  daTuple->column("rec" + head + "NumTTClusters",   (unsigned long long)allTTClusters->size());
  daTuple->column("rec" + head + "NumSpdHits",      (unsigned long long)allSpdDigits->size());
  daTuple->column("rec" + head + "NumPrsHits",      (unsigned long long)allPrsDigits->size());
  daTuple->column("rec" + head + "NumEcalHits",     (unsigned long long)allEcalDigits->size());
  daTuple->column("rec" + head + "NumHcalHits",     (unsigned long long)allHcalDigits->size());
  daTuple->column("rec" + head + "NumVeloTracks",   numVeloTracks);

  return StatusCode::SUCCESS;
}
//=============================================================================
// Plot Impact Parameters
//  1) If there are no PVs, store an errorCode          
//  2) If there is 1 PV, store the IP for this PV       
//  3) If there is more than 1 PV, store the minimum IP 
//  4) Store the PV coords? Or an ID for the PV?        
//=============================================================================
StatusCode LoopOnDaughters::plotIPs(const LHCb::Particle* da, const LHCb::RecVertex::Range prims, 
                                    Tuple daTuple, const std::string& head) const
{
  double IP=m_errorCode*mm, IPE=m_errorCode*mm; // Initialise to errorCode, only change if possible
  int pvID = -1;                                // Store ID of which PV was used to calculate IP
  if (prims.size()>0)                           // A requirement for Real Data... some events may not have any PVs
  {
    double tmpIP, tmpIPE; // Placeholder for each calculation of IP
    int tmpPvID=0;
    
    for (LHCb::RecVertex::Range::const_iterator ipv = prims.begin(); ipv != prims.end(); ++ipv )
    {
      tmpPvID++;
      if (msgLevel(MSG::DEBUG)) debug() << (*ipv)->position() << endmsg;
      /*
      StatusCode sc = m_extra->SignedImpactParameter(da, (*ipv), tmpIP, tmpIPE); // Calculate IP
      if (sc && tmpIP<abs(IP)) // If the calculation worked and the new IP is smaller than the previous one, update!
      {
        IP = tmpIP;
        if      (tmpIPE>0)         IPE = tmpIPE;
        else if (IPE!=m_errorCode) IPE = m_errorCode*mm; // Reset the IPE to an errorCode if no value was returned
        pvID = tmpPvID;
      }
      */
    }
  }
  
  daTuple->column("rec" + head + "minIP_whichPV", pvID);
  daTuple->column("rec" + head + "minIP",         IP);
  daTuple->column("rec" + head + "minIP_Error",   IPE); //Q: Previously stored IP/IPE, may need to do this again?

  return StatusCode::SUCCESS;
}

//=============================================================================
// Plot IP calculation (LEGACY: Unused since Apr 2011) 
//=============================================================================
StatusCode LoopOnDaughters::testIPs(const LHCb::Particle* da, const LHCb::RecVertex::Range prims,
                                    Tuple daTuple, const std::string& head) const
{
  //vector_IP = vector_ReferencePoint - vector_Vertex;
  //vector_IP -= vector_Momentum*(vector_ReferencePoint-vector_Vertex).vector_momentum/magnitude(vector_Momentum);
  Gaudi::XYZPoint vecP      = Gaudi::XYZPoint(da->momentum().x(), da->momentum().y(), da->momentum().z());
  Gaudi::XYZPoint refPoint  = da->referencePoint();

  daTuple->column("rec" + head + "Px",    da->momentum().x());
  daTuple->column("rec" + head + "Py",    da->momentum().y());
  daTuple->column("rec" + head + "Pz",    da->momentum().z());
  daTuple->column("rec" + head + "Ptest", sqrt(pow(da->momentum().x(),2)+pow(da->momentum().y(),2)+pow(da->momentum().z(),2)));
  daTuple->column("rec" + head + "refX",  da->referencePoint().x());  // Coord where momentum has been measured
  daTuple->column("rec" + head + "refY",  da->referencePoint().y());
  daTuple->column("rec" + head + "refZ",  da->referencePoint().z());
  
  for (LHCb::RecVertex::Range::const_iterator ipv = prims.begin(); ipv != prims.end(); ++ipv )
  {
    const LHCb::RecVertex *pv = *(ipv);
    Gaudi::XYZPoint pvPoint   = pv->position();
    daTuple->column("rec" + head + "PVx", pv->position().x());
    daTuple->column("rec" + head + "PVy", pv->position().y());
    daTuple->column("rec" + head + "PVz", pv->position().z());

    //double c = (refPoint.x()-pvPoint.x())*vecP.x();
    //c       += (refPoint.y()-pvPoint.y())*vecP.y();
    //c       += (refPoint.z()-pvPoint.z())*vecP.z();
    //c       /= pow(vecP.x(),2) + pow(vecP.y(),2) + pow(vecP.z(),2);

    //double ipX = refPoint.x() - pvPoint.x() - vecP.x()*c;
    //double ipY = refPoint.y() - pvPoint.y() - vecP.y()*c;
    //double ipZ = refPoint.z() - pvPoint.z() - vecP.z()*c;

    //Gaudi::XYZPoint vecIP = Gaudi::XYZPoint(ipX, ipY, ipZ);
    //double ip = sqrt(pow(ipX,2) + pow(ipY,2) + pow(ipZ,2));
    //info() << "My calculation of the impact parameter: " << ip << "mm." << endmsg;
    //Gaudi::XYZPoint vecIP = refPoint - pvPoint - vecP/(pow(vecP.x(),2) + pow(vecP.y(),2) + pow(vecP.z(),2));
  }
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// Calculate Isolation based on tracks
//=============================================================================
StatusCode LoopOnDaughters::trackIsolation(const LHCb::Particle* da, Tuple daTuple, const std::string& head, 
                                           const bool runTest) const
{
  //bool runTest = true;
  
  //### Identify the protoParticle and track associated with the current Particle ###//
  const LHCb::ProtoParticle* daProto = da->proto();
  const LHCb::Track* daTrack = daProto->track();

  //### Store Calorimeter data for the protoParticle associated with the given Particle                       ###//
  //### The info() function searches a vector of ExtraInfo for the given keyID and returns the data           ###//
  //### If it reaches the end of the vector without finding the keyID, it will return the error code supplied ###//
  daTuple->column("rec" + head + "EcalE", daProto->info( daProto->CaloEcalE, m_errorCode) );
  daTuple->column("rec" + head + "HcalE", daProto->info( daProto->CaloHcalE, m_errorCode) );

  //### Extract and loop over all other tracks in event ###//
  LHCb::Track::Vector otherTracks;
  otherTracks.clear();
  LHCb::Tracks* kotherTracks;

  if (exist<LHCb::Tracks>(LHCb::TrackLocation::Default))
  {
    kotherTracks = get<LHCb::Tracks>(LHCb::TrackLocation::Default);
    LHCb::Tracks::iterator it;
    for (it = kotherTracks->begin(); it != kotherTracks->end(); ++it) otherTracks.push_back(*it);
  }
  else info() << "There are no Tracks in the default location" << endmsg;

  double dist_PhiEta = 0;
  
  //### Define quantites to be recorded in tuple and set all values to zero ###//
  if (runTest) 
  {
    //IsolationMeasurement* isol = new IsolationMeasurement(m_numDivs, m_coneMax, "Isol", m_local);
    IsolationMeasurement* isol = new IsolationMeasurement("Isol", m_local);
    isol->Setup(m_numDivs, m_coneMax);
    
    //### Loop on tracks in event ###//
    for ( LHCb::Track::Vector::iterator itt = otherTracks.begin(); itt != otherTracks.end(); ++itt)
    {
      LHCb::Track* trk = *itt;
      if (trk == daTrack) continue;

      //### Calculate distance in Phi-Eta space ###//
      dist_PhiEta = da->momentum().phi() - trk->phi();
      if ( fabs(dist_PhiEta) > m_pi) dist_PhiEta = 2*m_pi - fabs(dist_PhiEta);
      dist_PhiEta = sqrt(pow(dist_PhiEta,2) + pow(da->momentum().Eta() - trk->momentum().Eta(), 2));
      debug() << "Distance in Phi-Eta space from this track to the Muon is " << dist_PhiEta << endmsg;

      isol->UpdateMeasurement(dist_PhiEta, trk);
    }
    isol->PlotMeasurement(daTuple, head, "Tracks");
    delete isol;
  }
  else
  {
    //### Identify different cone sizes based on m_coneDist ###//
    double test1 = double(m_coneMax)/m_numDivs;
    double coneSize[m_numDivs];

    for (int i = 0; i<m_numDivs; i++) coneSize[i] = m_coneMax - (i * test1);

    int num_tracks_tot[m_numDivs], num_tracks_T[m_numDivs], num_tracks_Velo[m_numDivs], num_tracks_TT[m_numDivs];
    double tot_meas[m_numDivs], tot_ghost[m_numDivs], tot_likelihood[m_numDivs], dist_PhiEta = 0;
    Gaudi::XYZVector ConeMomentum[m_numDivs];
    
    for (int i = 0; i<m_numDivs; i++)
    {
      num_tracks_tot[i] = 0;
      num_tracks_T[i]    = 0;
      num_tracks_Velo[i] = 0;
      num_tracks_TT[i]   = 0;
      tot_meas[i]        = 0;
      tot_ghost[i]       = 0;
      tot_likelihood[i]  = 0;
      
      ConeMomentum[i] = Gaudi::XYZVector(0, 0, 0);
    }
  
    //### Loop on tracks in event ###//
    for ( LHCb::Track::Vector::iterator itt = otherTracks.begin(); itt != otherTracks.end(); ++itt)
    {
      LHCb::Track* trk = *itt;
      if (trk == daTrack) continue;
      
      //### Calculate distance in Phi-Eta space ###//
      dist_PhiEta = da->momentum().phi() - trk->phi();
      if ( fabs(dist_PhiEta) > m_pi) dist_PhiEta = 2*m_pi - fabs(dist_PhiEta);
      dist_PhiEta = sqrt(pow(dist_PhiEta,2) + pow(da->momentum().Eta() - trk->momentum().Eta(), 2));
      debug() << "Distance in Phi-Eta space from this track to the Muon is " << dist_PhiEta << endmsg;
      
      //### Update variables if track is within the cone (but is not related to the particle itself!) ###//
      for (int i=0; i<m_numDivs; i++)
      {
        if (dist_PhiEta < coneSize[i] && trk != daTrack)
        {
          //### NB: hasT, hasVelo etc only check that track passes thru the station ###//
          //### Does not check that their are any hits on the track                 ###//
          //### To do this, would need to start by extracting trk->lhcbIDs()        ###//
          //### Then loop over these IDs and check what type of hit each one is!    ###//
          debug() << "Distance in Phi-Eta space from this track to the Muon is " << dist_PhiEta << endmsg;
          debug() << "Track has hit in T station?\t "    << trk->hasT() << endmsg;
          debug() << "Track has hit in Velo?\t "         << trk->hasVelo() << endmsg;
          debug() << "Track has hit in TT?\t "           << trk->hasTT() << endmsg;
          debug() << "Track has likelihood of: "         << trk->likelihood() << endmsg;
          debug() << "Track has ghostProbability of: "   << trk->ghostProbability() << endmsg;
          debug() << "Distance in Phi-Eta space from this track to the Muon is " << dist_PhiEta << endmsg;
          
          num_tracks_tot[i]++;
          num_tracks_T[i]    += trk->hasT();
          num_tracks_Velo[i] += trk->hasVelo();
          num_tracks_TT[i]   += trk->hasTT();
          tot_meas[i]        += trk->nMeasurements();
          if (trk->ghostProbability()!=999) tot_ghost[i] += trk->ghostProbability();
          tot_likelihood[i]  += trk->likelihood();
          
          ConeMomentum[i] += trk->momentum();
        }
      }
    }
    
    //### Plot data for tracks inside cones ###//
    char num [] = "123456789";
    for (int i=0; i<m_numDivs; i++)
    {
      daTuple->column("rec" + head + "Tracks_" + num[i] + "_ConeSize",   coneSize[i]);
      daTuple->column("rec" + head + "Tracks_" + num[i] + "_Num",        num_tracks_tot[i]);
      daTuple->column("rec" + head + "Tracks_" + num[i] + "_hasT",       num_tracks_T[i]);
      daTuple->column("rec" + head + "Tracks_" + num[i] + "_hasVelo",    num_tracks_Velo[i]);
      daTuple->column("rec" + head + "Tracks_" + num[i] + "_hasTT",      num_tracks_TT[i]);
      daTuple->column("rec" + head + "Tracks_" + num[i] + "_NumMeasure", tot_meas[i]);
      daTuple->column("rec" + head + "Tracks_" + num[i] + "_Likelihood", tot_likelihood[i]);
      daTuple->column("rec" + head + "Tracks_" + num[i] + "_GhostProb",  tot_ghost[i]);
      daTuple->column("rec" + head + "Tracks_" + num[i] + "_P",          ConeMomentum[i].R());
      daTuple->column("rec" + head + "Tracks_" + num[i] + "_Pt",         ConeMomentum[i].R()*sin(ConeMomentum[i].theta()));
      daTuple->column("rec" + head + "Tracks_" + num[i] + "_Phi",        ConeMomentum[i].phi());
      daTuple->column("rec" + head + "Tracks_" + num[i] + "_Eta",        ConeMomentum[i].Eta());
    }
  }
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// TEST - Use ProtoParticle instead of Tracks for Isolation Criteria
//=============================================================================
StatusCode LoopOnDaughters::protoParticleIsolation(const LHCb::Particle* da, Tuple daTuple, const std::string& head, 
                                                   bool runTest) const
{
  //bool runTest = false;
  
  //### Identify the protoParticle associated with the current Particle ###//
  const LHCb::ProtoParticle* daProto = da->proto();
  
  //### Extract and loop over all other tracks in event ###//
  LHCb::ProtoParticle::Vector otherProtos = getAllProtos();

  double dist_PhiEta = 0;
  
  if (runTest)
  {
    //IsolationMeasurement* isol = new IsolationMeasurement(m_numDivs, m_coneMax, "Isol", m_local);
    IsolationMeasurement* isol = new IsolationMeasurement("Isol", m_local);
    isol->Setup(m_numDivs, m_coneMax);
        
    //### Loop on tracks in event ###//
    for ( LHCb::ProtoParticle::Vector::iterator itpp = otherProtos.begin(); itpp != otherProtos.end(); ++itpp)
    {
      LHCb::ProtoParticle* prt = *itpp;
      if (prt == daProto) continue;

      //### Calculate distance in Phi-Eta space                                               ###//
      //### No Momentum data etc stored for ProtoParticle so need to extract associated track ###//
      //### So need to make sure there is an associated track                                 ###//
      if (prt->track())
      {
        dist_PhiEta = da->momentum().phi()-prt->track()->phi();
        if ( fabs(dist_PhiEta) > m_pi) dist_PhiEta = 2*m_pi - fabs(dist_PhiEta);
        dist_PhiEta = sqrt(pow(dist_PhiEta,2) + pow(da->momentum().Eta() - prt->track()->momentum().Eta(), 2));
        debug() << "Distance in Phi-Eta space from this protoParticle to the Muon is " << dist_PhiEta << endmsg;

        isol->UpdateMeasurement(dist_PhiEta, prt);
      }
      else debug() << "No Track associated to this protoParticle" << endmsg;
    }
    isol->PlotMeasurement(daTuple, head, "Protos");
    delete isol;
  }
  else
  {
    double test1 = double(m_coneMax)/m_numDivs;
    double coneSize[m_numDivs];
    for (int i = 0; i<m_numDivs; i++) coneSize[i] = m_coneMax - (i * test1);
    
    //### Define quantites to be recorded in tuple and set all values to zero ###//
    int num_protos_tot[m_numDivs], num_protos_T[m_numDivs], num_protos_Velo[m_numDivs], num_protos_TT[m_numDivs];
    double tot_protos_meas[m_numDivs], tot_protos_ghost[m_numDivs], tot_protos_likelihood[m_numDivs],
      tot_EcalE[m_numDivs], tot_HcalE[m_numDivs];
    Gaudi::XYZVector protos_ConeMomentum[m_numDivs];
    double dist_PhiEta = 0;
    
    for (int i = 0; i<m_numDivs; i++)
    {
      num_protos_tot[i]        = 0;
      num_protos_T[i]          = 0;
      num_protos_Velo[i]       = 0;
      num_protos_TT[i]         = 0;
      tot_protos_meas[i]       = 0;
      tot_protos_ghost[i]      = 0;
      tot_protos_likelihood[i] = 0;
      tot_EcalE[i]             = 0;
      tot_HcalE[i]             = 0;
      
      protos_ConeMomentum[i] = Gaudi::XYZVector(0,0,0);
    }
    
    //### Loop on tracks in event ###//
    for ( LHCb::ProtoParticle::Vector::iterator itpp = otherProtos.begin(); itpp != otherProtos.end(); ++itpp)
    {
      LHCb::ProtoParticle* prt = *itpp;
      if (prt == daProto) continue;
      
      //### Calculate distance in Phi-Eta space                                               ###//
      //### No Momentum data etc stored for ProtoParticle so need to extract associated track ###//
      //### So need to make sure there is an associated track                                 ###//
      if (prt->track())
      {
        dist_PhiEta = da->momentum().phi()-prt->track()->phi();
        if ( fabs(dist_PhiEta) > m_pi) dist_PhiEta = 2*m_pi - fabs(dist_PhiEta);
        dist_PhiEta = sqrt(pow(dist_PhiEta,2) + pow(da->momentum().Eta() - prt->track()->momentum().Eta(), 2));
        debug() << "Distance in Phi-Eta space from this protoParticle to the Muon is " << dist_PhiEta << endmsg;
        
        //### Update variables if track is within the cone (but is not related to the particle itself!) ###//
        for (int i=0; i<m_numDivs; i++)
        {
          if (dist_PhiEta < coneSize[i] && prt != daProto)
          {
            debug() << "Distance in Phi-Eta space from protoParticle to the Muon is less than" << coneSize[i] << endmsg;
            num_protos_tot[i]++;
            num_protos_T[i]          += prt->track()->hasT();
            num_protos_Velo[i]       += prt->track()->hasVelo();
            num_protos_TT[i]         += prt->track()->hasTT();
            tot_protos_meas[i]       += prt->track()->nMeasurements();
            tot_protos_ghost[i]      += prt->track()->ghostProbability();
            tot_protos_likelihood[i] += prt->track()->likelihood();
            debug() << "Ecal Energy for this protoParticle: " << prt->info( prt->CaloEcalE, m_errorCode ) << endmsg;
            debug() << "Hcal Energy for this protoParticle: " << prt->info( prt->CaloHcalE, m_errorCode ) << endmsg;
            tot_EcalE[i]   += prt->info( prt->CaloEcalE, 0);
            tot_HcalE[i]   += prt->info( prt->CaloHcalE, 0);
            
            protos_ConeMomentum[i] += prt->track()->momentum();
          }
        }
      }
      else debug() << "No Track associated to this protoParticle" << endmsg;
    }

    //### Plot data for tracks inside cones ###//
    char num [] = "123456789";
    for (int i=0; i<m_numDivs; i++)
    {
      daTuple->column("rec" + head + "Protos_" + num[i] + "_ConeSize",   coneSize[i]);
      daTuple->column("rec" + head + "Protos_" + num[i] + "_Num",        num_protos_tot[i]);
      daTuple->column("rec" + head + "Protos_" + num[i] + "_hasT",       num_protos_T[i]);
      daTuple->column("rec" + head + "Protos_" + num[i] + "_hasVelo",    num_protos_Velo[i]);
      daTuple->column("rec" + head + "Protos_" + num[i] + "_hasTT",      num_protos_TT[i]);
      daTuple->column("rec" + head + "Protos_" + num[i] + "_NumMeasure", tot_protos_meas[i]);
      daTuple->column("rec" + head + "Protos_" + num[i] + "_Likelihood", tot_protos_likelihood[i]);
      daTuple->column("rec" + head + "Protos_" + num[i] + "_GhostProb",  tot_protos_ghost[i]);
      daTuple->column("rec" + head + "Protos_" + num[i] + "_P",          protos_ConeMomentum[i].R()*MeV);
      daTuple->column("rec" + head + "Protos_" + num[i] + "_Pt",
                      protos_ConeMomentum[i].R()*sin(protos_ConeMomentum[i].theta())*MeV);
      daTuple->column("rec" + head + "Protos_" + num[i] + "_Phi",        protos_ConeMomentum[i].phi());
      daTuple->column("rec" + head + "Protos_" + num[i] + "_Eta",        protos_ConeMomentum[i].Eta());
      daTuple->column("rec" + head + "Protos_" + num[i] + "_EcalE",      tot_EcalE[i]);
      daTuple->column("rec" + head + "Protos_" + num[i] + "_HcalE",      tot_HcalE[i]);
    }
  }
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// plot daughters in given tuple                                                                  
//=============================================================================
StatusCode LoopOnDaughters::plotDaughter(const LHCb::Particle* da, const LHCb::RecVertex::Range prims, 
                                         Tuple daTuple, Tuple hitDistTuple, 
                                         const std::string& head, int runMC) const
{
  StatusCode sc = StatusCode::SUCCESS;
  debug() << "plotDaughters has been called for " << head << " with " << prims.size() << " vertices." << endmsg;
  char num [] = "123456789"; // Used several times to convert an integer to a string 
    
  //### Store Run Number and Event (L0) ID ###//
  const LHCb::ODIN* odin = get<LHCb::ODIN>(LHCb::ODINLocation::Default);
  const unsigned int runNum = odin->runNumber(), evNum = odin->eventNumber();
  daTuple->column("RunNumber", runNum);
  daTuple->column("EventID",   evNum);
  
  //##############################################//
  //### Check the trigger lines for this event (May need to update to include other lines?) ###//
  //##############################################//
  // Initialise decisions to errorCode, only update if possible
  int L0GlobalDecision=m_errorCode, L0MuonDecision=m_errorCode, L0MuonHighDecision=m_errorCode, L0DiMuonDecision=m_errorCode;
  int Hlt1GlobalDecision=m_errorCode, Hlt2GlobalDecision=m_errorCode, Hlt2ZmmDecision=m_errorCode; 
  
  if ( exist<LHCb::L0DUReport>(LHCb::L0DUReportLocation::Default) )
  {
    LHCb::L0DUReport *L0Rep = get<LHCb::L0DUReport>(LHCb::L0DUReportLocation::Default);
    /*
      unsigned int channelID=0;
      std::string channelName = "";
      while (1>0) {
      channelName = L0Rep->channelName(channelID);
      if (channelName.size()>20) break;
      if (L0Rep->chanelDecisionByName(channelName)) 
        info()<<"Trigger "<<channelName<<", decision "<<L0Rep->channelDecisionByName(channelName)<<endmsg;
      channelID++; }
    */
    L0GlobalDecision   = L0Rep->decision();
    L0MuonDecision     = L0Rep->channelDecisionByName("Muon");
    L0MuonHighDecision = L0Rep->channelDecisionByName("MuonHigh");
    L0DiMuonDecision   = L0Rep->channelDecisionByName("DiMuon");
  }
  else Warning("Can't get LHCb::L0DUReportLocation::Default (" + LHCb::L0DUReportLocation::Default + ")" );
  
  if ( exist<LHCb::HltDecReports>(LHCb::HltDecReportsLocation::Default) )
  {
    LHCb::HltDecReports *HltRep = get<LHCb::HltDecReports>(LHCb::HltDecReportsLocation::Default);
    std::vector< std::string > decNames = HltRep->decisionNames();
    //for (int i=0; i<decNames.size(); i++) {
    //  if ( HltRep->decReport(decNames[i])->decision()==1 ) info() << "Decision Name: " << decNames[i] << "; Decision: " 
    //                                                              << HltRep->decReport(decNames[i])->decision() << endmsg;
    //}
    if (HltRep->decReport("Hlt1Global")) Hlt1GlobalDecision = HltRep->decReport("Hlt1Global")->decision();
    if (HltRep->decReport("Hlt2Global")) Hlt2GlobalDecision = HltRep->decReport("Hlt2Global")->decision();
    //### TEST - Not sure if this is the right/only Z2MuMu line to consider! ###//
    if (HltRep->decReport("Hlt2DiMuonUnbiasedZmmDecision")) 
      Hlt2ZmmDecision = HltRep->decReport("Hlt2DiMuonUnbiasedZmmDecision")->decision();
  }
  else Warning("No HltDecReports at " + LHCb::HltDecReportsLocation::Default, StatusCode::FAILURE, 1);
      
  daTuple->column("rec" + head + "L0GlobalDecision",   L0GlobalDecision);
  daTuple->column("rec" + head + "L0MuonDecision",     L0MuonDecision);
  daTuple->column("rec" + head + "L0MuonHighDecision", L0MuonHighDecision);
  daTuple->column("rec" + head + "L0DiMuonDecision",   L0DiMuonDecision);
  daTuple->column("rec" + head + "Hlt1GlobalDecision", Hlt1GlobalDecision);
  daTuple->column("rec" + head + "Hlt2GlobalDecision", Hlt2GlobalDecision);
  daTuple->column("rec" + head + "Hlt2ZmmDecision",    Hlt2ZmmDecision);
    
  //##########################//
  //### Reconstructed Data ###//
  //##########################//

  debug() << "Daughter has P=" << da->p() << "MeV" << endmsg;
  if (msgLevel(MSG::DEBUG)) debug() << da->momentum() << endmsg;
    
  //### Plot reconstructed particle data ###//
  sc = plotReconstructedData(da, daTuple, head);
  if (!sc) return sc;
    
  //### Extract details on on LHCbIDs ###//
  sc = plotLHCbIDs(da->proto()->track()->lhcbIDs(), daTuple, head);
  if (!sc) return sc;

  //### Global Event Cut variables ###//
  sc = plotGlobalEventCutVariables(daTuple, head);
  if (!sc) return sc;
      
  //### Impact Parameters ###//
  sc = plotIPs(da, prims, daTuple, head);
  if (!sc) return sc;

  //sc = testIPs(da, prims, daTuple, head); //(LEGACY: Unused since Apr 2011)
  //if (!sc) return sc;
  
  //##########################//  
  //### Isolation Criteria ###//
  //##########################//

  bool tmpRunTest = true;
  
  sc = trackIsolation(da, daTuple, head, tmpRunTest);
  if (!sc) return sc;

  sc = protoParticleIsolation(da, daTuple, head, tmpRunTest);
  if (!sc) return sc;

  //##########################//
  //### Associated MC Data ###//
  //##########################//

  if (runMC==1)
  {
    //### Identify associated MC particle and plot data ###//
    Particle2MCLinker* assoc = new Particle2MCLinker(this, Particle2MCMethod::Links, "");
    const LHCb::MCParticle *mcda = assoc->firstMCP(da);
    
    if ( mcda!=NULL)
    {
      debug() << "AssocMCParticle is a " << mcda->particleID().pid() << endmsg;
      const LHCb::MCParticle *mcmo = mcda->mother();
      
      daTuple->column("assoc" + head + "ID",                                mcda->particleID().pid());
      daTuple->column("assoc" + head + "Mass",                              mcda->momentum().M());
      daTuple->column("assoc" + head + "P",                                 mcda->p());
      daTuple->column("assoc" + head + "Pt",                                mcda->pt());
      daTuple->column("assoc" + head + "E",                                 mcda->momentum().E());
      daTuple->column("assoc" + head + "Phi",                               mcda->momentum().phi());
      daTuple->column("assoc" + head + "Eta",                               mcda->momentum().Eta());
      if (mcmo != NULL) 
      {
        daTuple->column("assoc" + head + "MotherID",                        mcmo->particleID().pid());
        const LHCb::MCParticle *mcgm = mcmo->mother();
        if (mcgm != NULL) daTuple->column("assoc" + head + "GrandMotherID", mcgm->particleID().pid());
        else daTuple->column("assoc" + head + "GrandMotherID",              m_errorCode);
      }
      else 
      {
        daTuple->column("assoc" + head + "MotherID",                        m_errorCode);
        daTuple->column("assoc" + head + "GrandMotherID",                   m_errorCode);
      }
      IMCReconstructible* recoTool = tool<IMCReconstructible>("MCReconstructible");
      daTuple->column("assoc" + head + "Reco",                              recoTool->reconstructible(mcda));
      
      //### Calculate Impact Parameter - Make LHCb:Particle from LHCb:MCParticle so distanceCalculator can be used ###//
      const LHCb::MCVertex *vert = mcda->primaryVertex(), *vertRef = mcda->originVertex();
      Gaudi::XYZPoint point = vert->position(), pointRef = vertRef->position();
      
      // SHUST: This no longer compiles as of v35r0
      //LHCb::Particle ipPart = LHCb::Particle::Particle(mcda->particleID());
      LHCb::Particle ipPart = LHCb::Particle(mcda->particleID());
      
      ipPart.setMeasuredMass(mcda->momentum().M());
      ipPart.setMomentum(mcda->momentum());
      ipPart.setReferencePoint(pointRef);
      
      double assocIP, assocIPE;
      if (msgLevel(MSG::DEBUG)) debug() << vert->position() << endmsg;
      /*
      sc = m_extra->SignedImpactParameter(&ipPart, point, assocIP, assocIPE);
      if (!sc)
      {
        daTuple->column("assoc" + head + "IP",      m_errorCode*mm);
        daTuple->column("assoc" + head + "IPchi2",  m_errorCode*mm);
        return sc;
      }
      else
      {
        daTuple->column("assoc" + head + "IP", assocIP);
        daTuple->column("assoc" + head + "IPchi2", assocIPE);
      }
      */

      //### Get L0 decision ###//
      if ( exist<LHCb::L0DUReport>(LHCb::L0DUReportLocation::Default) )
      {
        LHCb::L0DUReport *L0Rep = get<LHCb::L0DUReport>(LHCb::L0DUReportLocation::Default);
        daTuple->column("assoc" + head + "L0Decision", L0Rep->decision());
      }
      else
      {
        Warning("Can't get LHCb::L0DUReportLocation::Default (" + LHCb::L0DUReportLocation::Default + ")" );
        daTuple->column("assoc" + head + "L0Decision", m_errorCode);
      }
      
      //### Get HLT decision ###//
      if ( exist<LHCb::HltDecReports>(LHCb::HltDecReportsLocation::Default) )
      {
        LHCb::HltDecReports *HltRep = get<LHCb::HltDecReports>(LHCb::HltDecReportsLocation::Default);
        if ( HltRep->decReport("Hlt1Global") ) daTuple->column("assoc" + head + "Hlt1Decision", 
                                                               HltRep->decReport("Hlt1Global")->decision());
        else                                   daTuple->column("assoc" + head + "Hlt1Decision", 
                                                               m_errorCode);
        if ( HltRep->decReport("Hlt2Global") ) daTuple->column("assoc" + head + "Hlt2Decision", 
                                                               HltRep->decReport("Hlt2Global")->decision());
        else                                   daTuple->column("assoc" + head + "Hlt2Decision", 
                                                               m_errorCode);
      }
      else
      {
        Warning("No HltDecReports at " + LHCb::HltDecReportsLocation::Default, StatusCode::FAILURE, 1);
        daTuple->column("assoc" + head + "Hlt1Decision", m_errorCode);
        daTuple->column("assoc" + head + "Hlt2Decision", m_errorCode);
      }
      
      counter("plotDaughterTest")++;
    }
    // runMC will not change during a job.. no need to set errorcodes if not running mc?
    else  //### output error code ###//
    {
      debug() << "There was no AssocMCParticle retrieved." << endmsg;
      daTuple->column("assoc" + head + "ID",            m_errorCode);
      daTuple->column("assoc" + head + "Mass",          m_errorCode*MeV);
      daTuple->column("assoc" + head + "P",             m_errorCode*MeV);
      daTuple->column("assoc" + head + "Pt",            m_errorCode*MeV);
      daTuple->column("assoc" + head + "E",             m_errorCode*MeV);
      daTuple->column("assoc" + head + "Phi",           double(m_errorCode)/100000);
      daTuple->column("assoc" + head + "Eta",           double(m_errorCode)/100000);
      daTuple->column("assoc" + head + "MotherID",      m_errorCode);
      daTuple->column("assoc" + head + "GrandMotherID", m_errorCode);
      daTuple->column("assoc" + head + "Reco",          m_errorCode);
      daTuple->column("assoc" + head + "IP",            m_errorCode*mm);
      daTuple->column("assoc" + head + "IPchi2",        m_errorCode*mm);
      daTuple->column("assoc" + head + "L0Decision",    m_errorCode);
      daTuple->column("assoc" + head + "Hlt1Decision",  m_errorCode);
      daTuple->column("assoc" + head + "Hlt2Decision",  m_errorCode);
    }

    delete assoc;
  }

  //##################################################################################################//
  //### TEST - For Jpsi - Loop over all protos and see if there are any hits in the muon detectors ###//
  //### Not exatcly what i'm supposed to do.. should be extrapolating the "muon" track and finding ###//
  //### the distance from this track to any hit in the moun detectors (for each station)           ###//
  //##################################################################################################//

  if (m_muHits) 
  {
    //### Extract all muonPIDs, where possible, find muon hits and store coordinates ##//
    LHCb::MuonPIDs* AllMuPIDs = get<LHCb::MuonPIDs>(LHCb::MuonPIDLocation::Default);
    std::vector<Gaudi::XYZPoint> allMuonHits; //Vector to store coordinates of all hits retrieved
    std::vector<Gaudi::XYZPoint> allMuonHitErrs; // TEST - Vector to store uncertainties on hit coords
    
    for (LHCb::MuonPIDs::iterator m = AllMuPIDs->begin(); m!=AllMuPIDs->end(); ++m)
    {
      LHCb::MuonPID* mupid = *m;
      if (mupid->IsMuonLoose())
      {
        debug() << "MuonPID IsMuonLoose() = " << mupid->IsMuonLoose() << endmsg;
        
        const LHCb::Track* muT = mupid->muonTrack();
        //const LHCb::Track* longT = mupid->idTrack();
        
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
            debug() << "and deltas (" <<  dx << ", " << dy << ", " << dz << ")" << endmsg;
            
            //### Store vector of XYZ coordinates for each muonTileID ###//
            Gaudi::XYZPoint point    = Gaudi::XYZPoint(x,y,z);
            Gaudi::XYZPoint errPoint = Gaudi::XYZPoint(dx,dy,dz); // TEST - store uncertainties on hits as a vector
            allMuonHits.push_back(point);
            
            // i^th member of allMuonHitErrs should be uncertainties on i^th member of allMuonHits
            // Maybe create a temporary tuple to test this?
            allMuonHitErrs.push_back(errPoint); 

            if (dx!=errPoint.x() || dy!=errPoint.y() || dz!=errPoint.z())
            {
              info() << "Uncertaintes on this point should be (" << dx << ", " << dy << ", " << dz << ")" << endmsg;
              info() << "Stored uncertainties are (" << errPoint.x() <<", "<< errPoint.y() <<", "<< errPoint.z() <<")"<< endmsg;
            }
          }
          else debug() << "LHCbID for this track is not a muon!" << endmsg;
        }
      }
    }
    
    if (allMuonHits.size()==0) debug() << "No hits found in the muon stations for this event." << endmsg;
    else                       debug() << "Retrieved " << allMuonHits.size() << " hits in the muon stations." << endmsg;
    
    //### TEST - make sure vector of hits and vector of uncertainties are the same size ###//
    if (allMuonHits.size()!=allMuonHitErrs.size()) 
    {
      info() << "WARNING: The uncertainites on muon hits have not been stored as expected!" << endmsg;
    }
        
    //### Loop over all muon hits and extract unique z coordinates ###//
    std::vector<double> allMuonHitsZ;
    for (unsigned int i=0; i<allMuonHits.size(); i++)
    {
      //### Extract z-coordinate of hit ###//
      debug()<<"Looking at hit coord ("<<allMuonHits[i].x()<<", "<<allMuonHits[i].y()<<", "<<allMuonHits[i].z()<<")"<<endmsg;
      int store=1;
      
      //### If this is the first z-coordinate to be extracted, it will automatically be stored ###//
      //### Otherwise, check that this z-coordinate has not already been stored ###//
      if (allMuonHitsZ.size()!=0)
      {
        debug() << "Checking against " << allMuonHitsZ.size() << " previously stored coordinates." << endmsg;
        //### Reject this z-coord if it has already been stored ###//
        for (unsigned int j=0; j<allMuonHitsZ.size(); j++) if(allMuonHits[i].z()==allMuonHitsZ[j]) store=0; 
      }
      
      if (store==1) allMuonHitsZ.push_back(allMuonHits[i].z());
      else          debug() << "Coordinate has previously been stored." << endmsg;
    }
    
    //### Define z-coordinate of Muon Stations ###//
    double M2_z = 15250, M3_z = 16450, M4_z = 17650, M5_z = 18850; // Note: M1_z = 12100
    
    //### Store the closest muon hit to the track, if any for each station ###//
    double minDist[5], dist;
    for (int i=0; i<5; i++) minDist[i]=-1*m_errorCode; // Initiate each minDist to a large value
    
    //### Extrapolate track to the z-coordinate of every hit and find the distance to that hit ###//
    int m;
    std::string det;
    for (unsigned int i=0; i<allMuonHitsZ.size(); i++)
    {
      dist=-1*m_errorCode; // initialise to some large value
      
      //### Identify which station the hit is in by its z-coordinate 
      if (allMuonHitsZ[i]<M2_z-1000)                            m=1;     
      if (allMuonHitsZ[i]>M2_z-1000&&allMuonHitsZ[i]<M3_z-1000) m=2;     
      if (allMuonHitsZ[i]>M3_z-1000&&allMuonHitsZ[i]<M4_z-1000) m=3;    
      if (allMuonHitsZ[i]>M4_z-1000&&allMuonHitsZ[i]<M5_z-1000) m=4;    
      if (allMuonHitsZ[i]>M5_z-1000)                            m=5; 
      
      //### Convert int to string ###//
      std::stringstream out;
      out << m;
      det = out.str();
      
      //### Find minimum distance from this hit to a track ###//
      debug() << "Looking for distance from track to hit with z = " << allMuonHitsZ[i] << " in M" << det << endmsg;
      dist = getMinDist(da->proto()->track(), allMuonHitsZ[i], allMuonHits, hitDistTuple, head);
      debug() << "M" << det << "; dist = " << dist << endmsg;
      
      //### Compare distance for this hit to the current smallest distance for this station. Update if needed ###//
      if (dist<minDist[m-1]) minDist[m-1]=dist;
    }
    
    for (int i=0; i<5; i++)
    {
      debug() << "Minimum distance found from hit to track for M" << num[i] << " is " << minDist[i] << endmsg;
      daTuple->column("rec" + head + "ClosestHitDist_M" + num[i],  minDist[i]);
    }
    
    //###########################################################################################//

    int muFlag[6];                        // Flags for track with at least one hit in a muon chamber
    for (int i=0; i<6; i++) muFlag[i]=0;  // Initialise to zero
    
    //### Should already have all protos in the event stored in otherProtos. Otherwise, need to extract all protos now ###//
    LHCb::ProtoParticle::Vector otherProtos = getAllProtos();  
    debug() << "Size of otherProtos vector is " << otherProtos.size() << endmsg;
    
    debug() << "###########################################################################################" << endmsg;
    for ( LHCb::ProtoParticle::Vector::iterator itpr = otherProtos.begin(); itpr != otherProtos.end(); ++itpr)
    {
      LHCb::ProtoParticle* prt = *itpr;
      if (prt->muonPID())
      {
        const LHCb::MuonPID* muPID = prt->muonPID();
        
        if(muPID->IsMuonLoose()!=0) 
        {
          debug() << "### Number of tracks which share hits is " << muPID->nShared() << " ###" << endmsg;
          
          const LHCb::Track* muTrack = muPID->muonTrack();
          debug() << "Have retrieved the muon track with momentum " << muTrack->p() << endmsg;
          
          if(muTrack->lhcbIDs().size()>0||muTrack->lhcbIDs().size()==0) 
          {
            debug() << "Can recover LHCbIDs container" << endmsg;
            debug() << " with size " << endmsg;
            debug() << muTrack->lhcbIDs().size() << endmsg;
          }
          else info() << "Cannot recover LHCbIDs container" << endmsg;
          
          LHCb::Track::LHCbIDContainer lhcbids = muTrack->lhcbIDs();
          debug() << "### Size of LHCbIDs Vector: " << muTrack->lhcbIDs().size() << " ###" << endmsg;
          
          for (LHCb::Track::LHCbIDContainer::const_iterator itID = lhcbids.begin(); itID!=lhcbids.end(); ++itID)
          {
            if(itID->isMuon())
            {
              debug() << "MUON type" << endmsg;
              LHCb::MuonTileID muID = itID->muonID();
              
              if(muID.station()==3&&(muID.quarter()==1||muID.region()==1)){
                debug() << "String " << muID.toString() << " ###" << endmsg;
              }
              
              debug()  << "### Muon passes through muon station " << muID.station()+1 << " ###" << endmsg;
              debug() << "### " << muID << " ###" << endmsg;
              
              muFlag[0]++;
              muFlag[muID.station()+1]++;
            }
          }
        }
        else 
        {
          /*
            LHCb::MuonCoords* coords;
            coords = get<LHCb::MuonCoords>(LHCb::MuonCoordLocation::MuonCoords);
            if ( coords==0 ) {
            info() << " Cannot retrieve MuonCoords " << endreq;
            //return StatusCode::FAILURE;
            }
            else 
            {
            //          for (LHCb::MuonCoords::iterator c = coords->begin(); c!=coords->end(), ++c) 
            //{
            //### Loop over coords and plot ###//
            debug() << "Success!!" << endmsg;
            //coordsTuple->column("rec" + head + "MuonCoordX", c->x());
            //}
            }
          */
          
          //### Extract the original track ###//
          const LHCb::Track* idTrack = prt->muonPID()->idTrack();
          debug() << "### Have found the original track ###" << endmsg;
          
          ///### Extract vector of "states" on track ###//
          LHCb::Track::StateContainer vecStates = idTrack->states();
          debug() << "Number of states on this track is " << vecStates.size() << endmsg;
          
          for (LHCb::Track::StateContainer::const_iterator itST=vecStates.begin(); itST!=vecStates.end(); ++itST)
          {
            LHCb::State* st = *itST;
            debug() << "This states has coordinates: (" << st->x() << ", " << st->y() << ", " << st->z() << ")" << endmsg;
          }
          
          //### Find state furthest into the detector ###//
          LHCb::State state = idTrack->closestState(20000);
          debug() << "Selecting state at (" <<  state.x() << ", " << state.y() << ", " << state.z() << ")" << endmsg;
          
          for (LHCb::Track::LHCbIDContainer::const_iterator itID = idTrack->lhcbIDs().begin();
               itID!=idTrack->lhcbIDs().end(); ++itID)
          {
            debug() << "Is muon??" << itID->isMuon() << endmsg;
            LHCb::MuonTileID testing = itID->muonID();
            debug() << "String " << testing.toString() << " ###" << endmsg;
          }
          
          //### Try to propogate the track instead of the state ###//
          LHCb::State state2 = idTrack->closestState(20000);
          debug() << "Closest state was (" << state2.x() << ", " << state2.y() << ", " << state2.z() << ")" << endmsg;
          
        }
      }
    }
    
    if(muFlag[0]!=0)
    {
      debug() << "Total number of muon hits in this event is " << muFlag[0] << " with" << endmsg;
      debug() << muFlag[1] << " in Station 1, " << muFlag[2] << " in Station 2, " << muFlag[3] << " in Station 3, " 
              << muFlag[4] << " in Station 4 and " << muFlag[5] << " in Station 5." << endmsg;
    }
    
    daTuple->column("rec" + head + "TotalMuonHitsInEvent", muFlag[0]);
    daTuple->column("rec" + head + "M1HitsInEvent",        muFlag[1]);
    daTuple->column("rec" + head + "M2HitsInEvent",        muFlag[2]);
    daTuple->column("rec" + head + "M3HitsInEvent",        muFlag[3]);
    daTuple->column("rec" + head + "M4HitsInEvent",        muFlag[4]);
    daTuple->column("rec" + head + "M5HitsInEvent",        muFlag[5]);
  }
  
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Get minimum distance from track to hit in a muon detector
//=============================================================================
double LoopOnDaughters::getMinDist(const LHCb::Track* daTrk, double Mz, std::vector<Gaudi::XYZPoint> allMuonHits, 
                                   Tuple hitDistTuple, const std::string& head) const
{
  StatusCode sc = StatusCode::SUCCESS;

  //### Need extrapolator tool to extrapolate this state ###//
  ITrackExtrapolator* trackExtrap = tool<ITrackExtrapolator>("TrackMasterExtrapolator");
  if(!trackExtrap)
  {
    err() << "Unable to get TrackMasterExtrapolator" << endmsg;
    return  StatusCode::FAILURE;
  }
  else debug() << "Track extrapolator tool retrieved." << endmsg;
  
  //### Since this function may be called several times in one event for the same track, need to record some id ###//
  int id=0;
  if (head=="All") id=1;                              // 1st call to loopOnDaughters for all protoparticle retrieved
  else if (head=="_MuPlus_"||head=="_MuMinus_") id=2; // 2nd call to loopOnDaughters after mother has been reconstructed
  else id=3;                                          // Just in case I add any other headings or there is some error 
  
  //### Find state closest to this muon chamber hit ###//
  LHCb::State daState = daTrk->closestState(Mz);
  debug() << "Selecting state at (" << daState.x() << ", " << daState.y() << ", " << daState.z() << ")" << endmsg;
  
  //### Extrapolate state from last point to Muon Station ###//
  sc = trackExtrap->propagate(daState, Mz);
  debug() << "State propogated to (" << daState.x() << ", " << daState.y() << ", " << daState.z() << ")" << endmsg;
  double distToHit = 0, minDist = -1*m_errorCode, minZ = -1*m_errorCode;

  //### Find number of hits with the current z-coordinate ###//
  int size=0;
  for (unsigned int i=0; i<allMuonHits.size(); i++)  if (allMuonHits[i].z()==Mz)  size++;
  debug() << "There are " << size << " hits which share this z coordinate." << endmsg;
  
  //### Loop over all retrieved muon hits ###//
  int num=0;
  for (unsigned int i=0; i<allMuonHits.size(); i++)
  {
    debug() << "Looping on hit " << i+1 << " from a total of " << allMuonHits.size() << " hits" << endmsg;
    Gaudi::XYZPoint hit = allMuonHits[i];
    debug() << "This hit has coords (" << hit.x() << ", " << hit.y() << ", " << hit.z() << ")" << endmsg;

    //### Only look at hits with this particular z-coordinate ###//
    if (hit.z()==Mz)
    {
      num++;
      debug() << "Looking at hit " << num << " of " << size << " hits for this z coordinate." << endmsg;
      debug() << "This hit has coords (" << hit.x() << ", " << hit.y() << ", " << hit.z() << ")" << endmsg;
      distToHit = sqrt( pow(daState.x() - hit.x(), 2) + pow(daState.y() - hit.y(), 2) + pow(daState.z() - hit.z(), 2) );
      debug() << "Distance from state to hit is " << distToHit << endmsg;
      
      //### Update minimum Distance if appropriate ###//
      if (distToHit<minDist)
      {
        minDist = distToHit;
        minZ = hit.z();
      }
      
      //### For now, store all values in a tuple ###//
      if(abs(hit.z()-daState.z())>0.00000000001) 
      {
        debug()<< "### Track extrapolator attempted to make state outside of detector! " 
               << "Returns previously iterated state at z = " << daState.z() << " ###" <<endmsg;
      }
      hitDistTuple->column("LoopID",      id);              //Record how function was called (see above) 
      hitDistTuple->column("TrackCharge", daTrk->charge()); //Can check if there is any difference between mu+ and mu-
      hitDistTuple->column("TrackP",      daTrk->p());
      hitDistTuple->column("TrackPt",     daTrk->pt());
      hitDistTuple->column("state_x",     daState.x());
      hitDistTuple->column("state_y",     daState.y());
      hitDistTuple->column("state_z",     daState.z());
      hitDistTuple->column("state_err_x", pow(daState.errX2(), 0.5));
      hitDistTuple->column("state_err_y", pow(daState.errY2(), 0.5));
      hitDistTuple->column("dist",        distToHit);
      hitDistTuple->column("numHits",     num);
      hitDistTuple->column("hit_x",       hit.x());
      hitDistTuple->column("hit_y",       hit.y());
      hitDistTuple->column("hit_z",       hit.z());

      //### Store minimum distance if appropriate, otherwise store an error code ###//      
      if (num==size)
      {
        hitDistTuple->column("minDist", minDist);
        hitDistTuple->column("minZ",    minZ);
      }
      else
      {
        hitDistTuple->column("minDist", m_errorCode*mm);
        hitDistTuple->column("minZ",    m_errorCode*mm);
      }
      hitDistTuple->write();
    }
  }
  
  debug() << "Minimum distance for this hit is " << minDist << endmsg;
  return minDist;
}

//=============================================================================
// Retrieve all Tracks in an event 
//=============================================================================
LHCb::ProtoParticle::Vector LoopOnDaughters::getAllProtos() const 
{
  StatusCode sc = StatusCode::SUCCESS;
 
  LHCb::ProtoParticle::Vector otherProtos;
  otherProtos.clear();
  LHCb::ProtoParticles* kotherProtos;

  //### No default location for ProtoParticles ###//
  //### Need to choose between: Charged, Upstream , Neutrals, HltCharged and HltNeutrals ###//
  if (exist<LHCb::ProtoParticles>(LHCb::ProtoParticleLocation::Charged))
  {
    kotherProtos = get<LHCb::ProtoParticles>(LHCb::ProtoParticleLocation::Charged);
    LHCb::ProtoParticles::iterator it;
    for (it = kotherProtos->begin(); it != kotherProtos->end(); ++it) otherProtos.push_back(*it);
  }
  else
  {
    info() << "There are no ProtoParticles in the charged location" << endmsg;
    counter("NoProtosCharged")++;
  }

  if (exist<LHCb::ProtoParticles>(LHCb::ProtoParticleLocation::Neutrals))
  {
    kotherProtos = get<LHCb::ProtoParticles>(LHCb::ProtoParticleLocation::Neutrals);
    LHCb::ProtoParticles::iterator it;
    for (it = kotherProtos->begin(); it != kotherProtos->end(); ++it) otherProtos.push_back(*it);
  }
  else
  {
    info() << "There are no ProtoParticles in the neutrals location" << endmsg;
    counter("NoProtosNeutral")++;
  }

  if (exist<LHCb::ProtoParticles>(LHCb::ProtoParticleLocation::Upstream))
  {
    kotherProtos = get<LHCb::ProtoParticles>(LHCb::ProtoParticleLocation::Upstream);
    LHCb::ProtoParticles::iterator it;
    for (it = kotherProtos->begin(); it != kotherProtos->end(); ++it) otherProtos.push_back(*it);
  }
  else
  {
    debug() << "There are no ProtoParticles in the upstream location" << endmsg;
    counter("NoProtosUpstream")++;
  }

  if (exist<LHCb::ProtoParticles>(LHCb::ProtoParticleLocation::HltCharged))
  {
    kotherProtos = get<LHCb::ProtoParticles>(LHCb::ProtoParticleLocation::HltCharged);
    LHCb::ProtoParticles::iterator it;
    for (it = kotherProtos->begin(); it != kotherProtos->end(); ++it) otherProtos.push_back(*it);
  }
  else
  {
    debug() << "There are no ProtoParticles in the HLTcharged location" << endmsg;
    counter("NoProtosHLTcharged")++;
  }

  if (exist<LHCb::ProtoParticles>(LHCb::ProtoParticleLocation::HltNeutrals))
  {
    kotherProtos = get<LHCb::ProtoParticles>(LHCb::ProtoParticleLocation::HltNeutrals);
    LHCb::ProtoParticles::iterator it;
    for (it = kotherProtos->begin(); it != kotherProtos->end(); ++it) otherProtos.push_back(*it);
  }
  else
  {
    debug() << "There are no ProtoParticles in the HLTneutrals location" << endmsg;
    counter("NoProtosHLTneutrals")++;
  }
  
  return(otherProtos);
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode LoopOnDaughters::finalize() {

  if (msgLevel(MSG::DEBUG)) debug() << "==> Finalize" << endmsg;

  //return DVAlgorithm::finalize(); 
  return DaVinciTupleAlgorithm::finalize(); 
  //return DaVinciAlgorithm::finalize(); 
} 

//=============================================================================
