// $Id: Rec.cpp, 12/04/2011 rhuston Exp $

//#####################//
//### Include files ###//
//#####################// 

//### from Gaudi ###//
#include "GaudiKernel/AlgFactory.h" 
//### local ###//
#include "Rec.h"
#include "LoopOnDaughters.h"
//### New Filter ###//
#include "Kernel/ParticleFilters.h"
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
//### MC ###//
#include "Event/MCParticle.h"
//### Extract L0 Decision ###//
#include "Event/L0DUReport.h"
#include "HltDecReports.h"
//### Extract Run Number and Event (L0) ID ###//
#include "Event/ODIN.h"
//### TEST - Does this make vertexFitter work? ###//
#include "Kernel/IVertexFit.h"
#include "MCInterfaces/IMCReconstructible.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

using namespace Gaudi::Units;
using namespace boost::lambda;

//-----------------------------------------------------------------------------
// Implementation file for class : Rec
//
// 12-04-2011 : Shane Huston
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( Rec );

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
// Initialise variables defined in header file
Rec::Rec( const std::string& name, ISvcLocator* pSvcLocator)
  : DaVinciTupleAlgorithm ( name , pSvcLocator ),
    m_errorCode(-1234567)
{
  m_local = pSvcLocator;
}

//=============================================================================
// Destructor
//=============================================================================
Rec::~Rec() { } 

//=============================================================================
// Initialization
//=============================================================================
StatusCode Rec::initialize() {
  StatusCode sc = DaVinciTupleAlgorithm::initialize(); 
  if ( sc.isFailure() ) return sc;
  if (msgLevel(MSG::DEBUG)) debug() << "==> Initialize" << endmsg;
  info() << "Initialize called for Rec" << endmsg;

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode Rec::execute() {
  if (msgLevel(MSG::DEBUG)) debug() << "==> Execute" << endmsg;
  StatusCode sc = StatusCode::SUCCESS ;
  setFilterPassed(true);   // Mandatory. Set to true if event is accepted. 
  err() << "Execute phase of Rec.cpp has been called." << endmsg; // Execute phase should not be called 
  return StatusCode::SUCCESS;
}

//=============================================================================
// Extract daughters and PVs, activate LoopOnDaughter and makeMother
//=============================================================================
LHCb::Particle::Vector Rec::loopOnRec(const LHCb::Particle::ConstVector& daughters, const LHCb::RecVertex::Range prims, 
                                      LHCb::ParticleID motherID, Tuple recAllTuple, Tuple hitDistTuple, Tuple motherTuple, 
                                      Tuple candTuple,  int runMC)
{

  StatusCode sc = StatusCode::SUCCESS ;
  LHCb::Particle::Vector mothers;
  
  //### Declare variables needed to call LoopOnDaughters::loopOnDaughters ###// 
  if (daughters.size()!=0)
  {
    debug() << "### There are " << daughters.size() << " possible Muons in this event (EventID "
            << counter("EventID").nEntries() << ") ###" << endmsg;
    counter("MuonsPresent")++;
    
    const std::string& strLoop = "Loop";
    LoopOnDaughters* Loop = new LoopOnDaughters(strLoop, m_local);
    sc = Loop->loopOnDaughters(daughters, prims, recAllTuple, hitDistTuple, runMC);
    delete Loop;
    if(!sc) return mothers;
        
    mothers = makeMother(daughters, prims, motherID, motherTuple, candTuple, hitDistTuple, runMC);
    if (!sc) return mothers;
  }
  else debug() << "### There are no possible muons in the event!! ###" << endmsg;
  
  return mothers;
}

//=============================================================================
// Make Mother
//=============================================================================
LHCb::Particle::Vector Rec::makeMother(const LHCb::Particle::ConstVector& daughters, const LHCb::RecVertex::Range prims, 
                                       LHCb::ParticleID motherID, Tuple motherTuple, Tuple candTuple, Tuple hitDistTuple, 
                                       int runMC)
{
  StatusCode sc = StatusCode::SUCCESS;
  LHCb::Particle::Vector mothers;
  
  getAndStoreEventNumber(motherTuple);
  getAndStoreRunNumberAndL0EventID(motherTuple);
  
  //### Seperate daughters into positive and negative ###//
  LHCb::Particle::ConstVector DaPluses, DaMinuses;
  size_t nDaughters = DaVinci::filter(daughters, bind(&LHCb::Particle::charge,_1)<0, DaMinuses);
  debug() << "Number of muMinus is " << nDaughters << endmsg;
  if (nDaughters>0) // TODO: Is this necessary? Why would this be negative? If this can happen, shouldn't it be treated as an error instead?
  {
    nDaughters += DaVinci::filter(daughters, bind(&LHCb::Particle::charge,_1)>0, DaPluses);
  }
  debug() << "Total number of muons is " << nDaughters << endmsg;
  debug() << "Total number of dimuons is " << DaPluses.size()*DaMinuses.size() << "." << endmsg;

  storeNumberOfMuonsAndPrimaryVertices(nDaughters, prims, motherTuple);
  
  int numCandidates = 0;
  LHCb::Track::Vector longTracks = extractAllLongTracksForEvent();
  
  //### Loop over all dimuons ###//
  for (LHCb::Particle::ConstVector::const_iterator imp = DaPluses.begin() ;
       imp != DaPluses.end(); ++imp )
  {
    const LHCb::Particle* daPlus = *imp;
    
    for (LHCb::Particle::ConstVector::const_iterator imm = DaMinuses.begin();
         imm != DaMinuses.end(); ++imm)
    {
      const LHCb::Particle* daMinus = *imm;

      storeNumberOfMuonsAndLongTracksPerEvent(nDaughters, longTracks, motherTuple);
      getAndStoreDiMuonInvariantMass(daPlus, daMinus, motherTuple);
      
      //### Identify the tracks associated with the current dimuon pair ###// 
      const LHCb::Track* muTrack1 = daPlus->proto()->track();
      const LHCb::Track* muTrack2 = daMinus->proto()->track();
      
      //### Loop over all other tracks in event ###//
      Gaudi::XYZVector otherLongTracks_P = Gaudi::XYZVector(0, 0, 0);
      for ( LHCb::Track::Vector::iterator itt = longTracks.begin(); itt != longTracks.end(); ++itt)
      {
        LHCb::Track* trk = *itt;
        if (trk != muTrack1 && trk != muTrack2) otherLongTracks_P += trk->momentum();
      }
      motherTuple->column("otherLongTracks_P",  otherLongTracks_P.R());
      motherTuple->column("otherLongTracks_Pt", otherLongTracks_P.R()*sin(otherLongTracks_P.theta()));
      
      //### Count number of candidates per event ###//
      numCandidates++;
      motherTuple->column("CandidateNumber", numCandidates);
      
      //### Find and plot diparticle mass of associated MC particles (if MC data is available) ###//
      if (runMC==1)
      {
        debug() << "Finding Invariant mass of MC DiMuon associated with reconstructed signal." << endmsg;
        Particle2MCLinker* assoc = new Particle2MCLinker(this, Particle2MCMethod::Links, ""); // Retrieves associated MC particle
        const LHCb::MCParticle *mcp = assoc->firstMCP(daPlus), *mcm = assoc->firstMCP(daMinus);
        Gaudi::LorentzVector twoDaMC;
        if (mcp !=NULL && mcm != NULL)
        {
          twoDaMC = (mcp)->momentum() + (mcm)->momentum();
          motherTuple->column("mc_DiMuon_InvMass", twoDaMC.M());
        }
        else motherTuple->column("mc_DiMuon_InvMass", m_errorCode*MeV);
        
        delete assoc;
      }
      
      //calculateImpactParametersWithReconstructedPrimaryVertices(prims, daPlus, daMinus, muTrack1, myTrack2, motherTuple);
      //getAndStoreDiMuonDistanceOfClosestApproach(daPlus, daMinus, motherTuple);
      fitVertexAndStoreImpactParameterData(motherID, daPlus, daMinus, motherTuple, hitDistTuple);
      
      //### Mandatory. Set to true if event is accepted. ###//
      setFilterPassed(true);
      
      //### Declare the mother to the PhysDesktop ###//
      mothers.push_back( Mother.clone() );

      if (msgLevel(MSG::DEBUG)) debug() << "Saved mother " << Mother.particleID().pid()
                                        << " to desktop" << endmsg;
      
      bool plottedDaughters = plotDaughters("plotDaughter", daPlus, daMinus, prims, motherTuple, hitDistTuple, runMC);
      if (!plottedDaughters) return mothers;
      
      motherTuple->write();
      counter ("Mothers")++; //Booked just like plots, counts # of candidates 
    }
  }

  storeAndWriteNumberOfCandidatesPerEvent(numCandidates, candTuple);
  
  return mothers;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode Rec::finalize() {

  if (msgLevel(MSG::DEBUG)) debug() << "==> Finalize" << endmsg;

  return DaVinciTupleAlgorithm::finalize(); 
} 
//=============================================================================

void getAndStoreEventNumber(Tuple tuple)
{
  counter("EventNumber")++;
  tuple->column("EventNumber", (unsigned long long)counter("EventNumber").nEntries());
}

void getAndStoreRunNumberAndL0EventID(Tuple tuple)
{
  const LHCb::ODIN* odin = get<LHCb::ODIN>(LHCb::ODINLocation::Default);
  const unsigned int runNum = odin->runNumber(), evNum = odin->eventNumber();
  tuple->column("RunNumber", runNum);
  tuple->column("EventID",   evNum);
}

void storeNumberOfMuonsAndPrimaryVertices(size_t nDaughters, const LHCb::RecVertex::Range prims, Tuple tuple)
{
  tuple->column("numMuons", (unsigned long long)nDaughters);
  tuple->column("numPVs",   (unsigned long long)prims.size());
}

void getAndStoreDiMuonInvariantMass(const LHCb::Particle* muPlus, const LHCb::Particle* muMinus, Tuple tuple)
{
  Gaudi::LorentzVector diMuonMomentum = muPlus->momentum() + muMinus->momentum();
  debug() << "Rec two daughter mass is " << diMuonMomentum.M()/GeV << " GeV" << endmsg;
  tuple->column("rec_DiMuon_InvMass", diMuonMomentum.M());
}

// This will be used to find the total energy and Pt of the event
LHCb::Track::Vector extractAllLongTracksForEvent()
{
  LHCb::Track::Vector longTracks;
  longTracks.clear();

  LHCb::Tracks* allTracks;

  if (exist<LHCb::Tracks>(LHCb::TrackLocation::Default))
  {
    allTracks = get<LHCb::Tracks>(LHCb::TrackLocation::Default);
    LHCb::Tracks::iterator it;
    for (it = allTracks->begin(); it != allTracks->end(); ++it){
      LHCb::Track* track = *it;
      if (track->type()==3) longTracks.push_back(track);
    }
  }
  else info() << "There are no Tracks in the default location" << endmsg;

  return longTracks;
}

void storeNumberOfMuonsAndLongTracksPerEvent(size_t nDaughters, LHCb::Track::Vector longTracks, Tuple motherTuple)
{
  tuple->column("numMuonsPerEvent",      (unsigned long long)nDaughters);
  tuple->column("numLongTracksPerEvent", (unsigned long long)longTracks.size());
}

void storeImpactParameterData(double fitIPplus, double fitIPminus, double fitIPEplus, double fitIPEminus,
                              double fitIPtot, double fitIPEtot, Tuple tuple)
{
  tuple->column("rec_FittedVertex_IP_Total",    fitIPtot);
  tuple->column("rec_FittedVertex_IPE_Total",   fitIPEtot);
  tuple->column("rec_FittedVertex_IP_MuPlus",   fitIPplus);
  tuple->column("rec_FittedVertex_IPE_MuPlus",  fitIPEplus);
  tuple->column("rec_FittedVertex_IP_MuMinus",  fitIPminus);
  tuple->column("rec_FittedVertex_IPE_MuMinus", fitIPEminus);
}

// TODO: This should really return an object with all of the IP data
// Then we would have a separate method to write that data
void calculateImpactParametersWithReconstructedPrimaryVertices(
        const LHCb::RecVertex::Range prims, const LHCb::Particle* muPlus, const LHCb::Particle* muMinus,
        const LHCb::Track* muPlusTrack, const LHCb::Track* muMinusTrack, Tuple tuple)
{
  //#################################################################//
  //### Total IP for DiMuons with reconstructed PV:               ###//
  //###  1) If there are no PVs, store an errorCode               ###//
  //###  2) If there is 1 PV, store the totalIP for this PV       ###//
  //###  3) If there is more than 1 PV, store the minimum totalIP ###//
  //###  4) Store the PV coords? Or an ID for the PV used?        ###//
  //### Refit PV :                                                ###//
  //###  1) Refit the PV used above, excluding the two daughters  ###//
  //###  2) Calculate new IP for the refitted vertex              ###//
  //#################################################################//

  //### Declare some variables ###//
  IDistanceCalculator* distTool = tool<IDistanceCalculator>("LoKi::DistanceCalculator");// Tool to calculate distances
  double recIPplus = m_errorCode*mm,recIPEplus = m_errorCode*mm, recIPminus = m_errorCode*mm,recIPEminus = m_errorCode*mm;
  double recIPtot = m_errorCode*mm, recIPEtot = m_errorCode*mm;
  double refitIPplus  = m_errorCode*mm, refitIPminus  = m_errorCode*mm, refitIPtot  = m_errorCode*mm;
  double refitIPEplus = m_errorCode*mm, refitIPEminus = m_errorCode*mm, refitIPEtot = m_errorCode*mm;
  double recPV_x = m_errorCode*mm, recPV_y = m_errorCode*mm, recPV_z = m_errorCode*mm,
    refitPV_x = m_errorCode*mm, refitPV_y = m_errorCode*mm, refitPV_z = m_errorCode*mm;

  int pvID=-1;                                                // Store ID of which PV was used to calculate IP
  LHCb::RecVertex refittedPV;                                 // Placeholder for refitted PV
  int refitComplete = 1;                                      // Record whether refit of PV has succeeded (=1) or failed (=0)

  //### Access and loop over reconstructed PVs ###//
  //### Note: In Real Data, some events may not have any PVs
  if (prims.size()>0)
  {
    double tmpIPplus, tmpIPEplus, tmpIPminus, tmpIPEminus, tmpIPtot, tmpIPEtot;
    int tmpPvID=0;
    for (LHCb::RecVertex::Range::const_iterator ipv = prims.begin(); ipv != prims.end(); ++ipv )
    {
      tmpPvID++; // Keep track of which PV is being used

      //### Find IP for muPlus for this PV ###//
      StatusCode scIPS = m_extra->ImpactParameterSum(muPlus, muMinus, *ipv,
                                                     tmpIPplus, tmpIPminus, tmpIPEplus, tmpIPEminus, tmpIPtot, tmpIPEtot);
      //StatusCode scPlus  = m_extra->SignedImpactParameter(muPlus, (*ipv), tmpIPplus,  tmpIPEplus);
      //StatusCode scMinus = m_extra->SignedImpactParameter(muMinus, (*ipv), tmpIPminus, tmpIPEminus);

      //### If IP's found for both muons, check total IP for this PV and update final result if needed ###//
      //if (scPlus && scMinus && tmpIPplus+tmpIPminus < abs(recIPplus+recIPminus) )
      if (scIPS && abs(tmpIPtot) < abs(recIPtot))
      {
        recIPplus  = tmpIPplus;
        recIPminus = tmpIPminus;
        recIPtot   = tmpIPtot;
        recIPEtot  = tmpIPEtot;

        if (tmpIPEplus>0)
          recIPEplus  = tmpIPEplus;
        else if (recIPEplus!=m_errorCode)
          recIPEplus  = m_errorCode*mm; // Reset the IPE to an errorCode if no value was returned
        if (tmpIPEminus>0)
          recIPEminus = tmpIPEminus;
        else if (recIPEminus!=m_errorCode)
          recIPEminus = m_errorCode*mm; // Reset the IPE to an errorCode if no value was returned

        pvID = tmpPvID;
      }
    }

    //### Extract coords of desired PV ###//
    recPV_x = prims[pvID-1]->position().x();
    recPV_y = prims[pvID-1]->position().y();
    recPV_z = prims[pvID-1]->position().z();

    //### Now that the correct PV has been identified, refit it without the two daughters ###//
    std::vector<const LHCb::Track*> twoDaughters;
    twoDaughters.push_back(muPlusTrack);
    twoDaughters.push_back(muMinusTrack);
    m_pvtool = tool<IPVOfflineTool>("PVOfflineTool");
    StatusCode scRefit = m_pvtool->reDoSinglePV(prims[pvID-1]->position(), twoDaughters, refittedPV);
    if (!scRefit) refitComplete = 0;

    //### Store coords of refitted vertex ###//
    refitPV_x = refittedPV.position().x();
    refitPV_y = refittedPV.position().y();
    refitPV_z = refittedPV.position().z();

    //### calculate new IP for each daughter based on the refitted PV ###//
    //m_extra->ImpactParameterSum(muPlus, muMinus, &refittedPV, refitIPplus, refitIPminus,
    //                         refitIPEplus, refitIPEminus, refitIPtot, refitIPEtot);
    //StatusCode scRefitPlus  = m_extra->SignedImpactParameter(muPlus, &refittedPV, refitIPplus,  refitIPEplus);
    //StatusCode scRefitMinus = m_extra->SignedImpactParameter(muMinus, &refittedPV, refitIPminus, refitIPEminus);
    //if (!scRefitPlus) {
    //  refitIPplus = m_errorCode*mm;
    //  refitIPEplus = m_errorCode*mm;
    //}
    //else if (refitIPEplus==0) refitIPEplus = m_errorCode*mm;
    //if (!scRefitMinus) {
    //  refitIPminus = m_errorCode*mm;
    //  refitIPEminus = m_errorCode*mm;
    //}
    //else if (refitIPEminus==0) refitIPEminus = m_errorCode*mm;
  }

  tuple->column("rec_DiMuon_recIP_whichPV",    pvID);
  tuple->column("rec_DiMuon_recPV_X",          recPV_x);
  tuple->column("rec_DiMuon_recPV_Y",          recPV_y);
  tuple->column("rec_DiMuon_recPV_Z",          recPV_z);
  tuple->column("rec_DiMuon_recIP_Total",      recIPtot);
  tuple->column("rec_DiMuon_recIPE_Total",     recIPEtot);
  tuple->column("rec_DiMuon_recIP_MuPlus",     recIPplus);
  tuple->column("rec_DiMuon_recIP_MuMinus",    recIPminus);
  tuple->column("rec_DiMuon_recIPE_MuPlus",    recIPEplus);
  tuple->column("rec_DiMuon_recIPE_MuMinus",   recIPEminus);
  tuple->column("rec_DiMuon_refitComplete",    refitComplete); //Refit can fail while still returning a new PV
  tuple->column("rec_DiMuon_refitPV_X",        refitPV_x);
  tuple->column("rec_DiMuon_refitPV_Y",        refitPV_y);
  tuple->column("rec_DiMuon_refitPV_Z",        refitPV_z);
  tuple->column("rec_DiMuon_refitIP_Total",    refitIPtot);
  tuple->column("rec_DiMuon_refitIP_MuPlus",   refitIPplus);
  tuple->column("rec_DiMuon_refitIP_MuMinus",  refitIPminus);
  tuple->column("rec_DiMuon_refitIPE_Total",   refitIPEtot);
  tuple->column("rec_DiMuon_refitIPE_MuPlus",  refitIPEplus);
  tuple->column("rec_DiMuon_refitIPE_MuMinus", refitIPEminus);
}

void getAndStoreDiMuonDistanceOfClosestApproach(const LHCb::Particle* muPlus, const LHCb::Particle* muMinus, Tuple tuple)
{
  double doca=m_errorCode*mm, docaChi2=m_errorCode*mm;
  sc = distTool->distance(muPlus, muMinus, doca, docaChi2);
  tuple->column("rec_DiMuon_DOCA",      doca);
  tuple->column("rec_DiMuon_DOCA_Chi2", docaChi2);
}

void fitVertexAndStoreImpactParameterData(LHCb::ParticleID motherID, const LHCb::Particle* muPlus,
                                          const LHCb::Particle* muMinus, Tuple motherTuple, Tuple hitDistTuple)
{
  //### Now make the vertex by calling the Vertex Fitter (returns vertex and mother particle) ###//
  LHCb::Vertex DaDaVertex;
  LHCb::Particle Mother(motherID); // Fixed address in stack of mother to be created - Will be replaced by pointer to heap
  double fitIPplus=m_errorCode*mm, fitIPminus=m_errorCode*mm, fitIPEplus=m_errorCode*mm, fitIPEminus=m_errorCode*mm;
  double fitIPtot=m_errorCode*mm, fitIPEtot = m_errorCode*mm;
  IVertexFit* testTool = tool<IVertexFit>("LoKi::VertexFitter"); // TEST - Tool to fit vertices
  StatusCode scFit = testTool->fit(*(muPlus),*(muMinus),DaDaVertex,Mother); // Seems to work but need to verify
  //StatusCode scFit = vertexFitter()->fit(*(muPlus),*(muMinus),DaDaVertex,Mother); // Old method in Z2TauTau implementation

  if (!scFit)
  {
    //### Output error code and plot any remaining data ###//
    Warning("Fit error").ignore();
    err() << "Error while fitting" << endmsg;
    motherTuple->column("rec_DiMuon_Chi2", -m_errorCode*mm);
    storeImpactParameterData(fitIPplus, fitIPminus, fitIPEplus, fitIPEminus, fitIPtot, fitIPEtot, motherTuple);
    bool plottedDaughters = plotDaughters("plotDaughter", muPlus, muMinus, prims, motherTuple, hitDistTuple, runMC);
    if (!plottedDaughters) return mothers;

    motherTuple->write();
    counter ("FitError")++;
    continue;
  }

  if (msgLevel(MSG::DEBUG)) debug() << "Vertex fit at " << DaDaVertex.position()/cm
                                    << " with chi^2 " << DaDaVertex.chi2() << endmsg;
  motherTuple->column("rec_DiMuon_Chi2", DaDaVertex.chi2());

  //### Find IP for dimuon with the fitted vertex  ###//
  //m_extra->ImpactParameterSum(muPlus, muMinus, &DaDaVertex, fitIPplus, fitIPminus,
  //                          fitIPEplus, refitIPminus, fitIPtot, fitIPEtot);
  //sc = m_extra->SignedImpactParameter(muPlus, &DaDaVertex, fitIPplus,  fitIPEplus);  // Calculate for muPlus
  //sc = m_extra->SignedImpactParameter(muMinus, &DaDaVertex, fitIPminus, fitIPEminus); // Calculate for muMinus

  storeImpactParameterData(fitIPplus, fitIPminus, fitIPEplus, fitIPEminus, fitIPtot, fitIPEtot, motherTuple);
}

bool plotDaughters(const std::string& counterName, const LHCb::Particle* muPlus, const LHCb::Particle* muMinus,
                   const LHCb::RecVertex::Range prims, Tuple motherTuple, Tuple hitDistTuple, int runMC)
{
  // Note: It may be inefficient to reinitialise (and delete) this for each DiMuon pair
  // It may be better to create a private m_Da variable and initialise it once
  // It looks like I've tried to do that before for IsolationMeasurement and for Extra and may require usage of #ifndef
  // It may be better not to try to do this until I'm able to access the DaVinci compiler again
  const std::string& strDa = "Loop";
  LoopOnDaughters *Da = new LoopOnDaughters(strDa, m_local);

  sc = Da->plotDaughter(muPlus, prims, motherTuple, hitDistTuple,
                              "_MuPlus_", runMC);
  if(!sc) return false;
  counter(counterName)++;

  sc = Da->plotDaughter(muMinus, prims, motherTuple, hitDistTuple,
                        "_MuMinus_", runMC);
  if (!sc) return return false;
  counter(counterName)++;

  delete Da; // Deallocate memory from the heap for LoopOnDaughters class
  return true;
}

void storeAndWriteNumberOfCandidatesPerEvent(int numCandidates, Tuple tuple)
{
  debug() << "Event Number: " << counter("EventNumber").nEntries() << endmsg;
  debug() << "Number of candidates in this event: " << numCandidates << endmsg;
  tuple->column("EventNumber", (unsigned long long)counter("EventNumber").nEntries());
  tuple->column("numCandidates", numCandidates);
  tuple->write();
}