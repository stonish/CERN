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
  
  //### Declare object to call LoopOnDaughters::PlotDaughters ###//
  const std::string& strDa = "Loop";
  LoopOnDaughters *Da = new LoopOnDaughters(strDa, m_local);
  
  //### Store Event Number in sequence###// 
  counter("EventNumber")++;
  motherTuple->column("EventNumber", (unsigned long long)counter("EventNumber").nEntries());
  
  //### Store Run number and event number (L0 Event ID) ###//
  const LHCb::ODIN* odin = get<LHCb::ODIN>(LHCb::ODINLocation::Default);
  const unsigned int runNum = odin->runNumber(), evNum = odin->eventNumber();
  motherTuple->column("RunNumber", runNum);
  motherTuple->column("EventID",   evNum);
  
  //### Seperate daughters into positive and negative ###//
  LHCb::Particle::ConstVector DaPlus, DaMinus;
  size_t nDaughters = DaVinci::filter(daughters, bind(&LHCb::Particle::charge,_1)<0, DaMinus);
  debug() << "Number of muMinus is " << nDaughters << endmsg;
  if (nDaughters>0) nDaughters += DaVinci::filter(daughters, bind(&LHCb::Particle::charge,_1)>0, DaPlus);
  debug() << "Total number of muons is " << nDaughters << endmsg;
  debug() << "Total number of dimuons is " << DaPlus.size()*DaMinus.size() << "." << endmsg;
  motherTuple->column("numMuons", (unsigned long long)nDaughters);
  motherTuple->column("numPVs",   (unsigned long long)prims.size());
  
  int numCandidates = 0;
  
  //### Extract all tracks in the event (used to find total energy, Pt of event) ###// 
  LHCb::Track::Vector longTracks;
  longTracks.clear();
  LHCb::Tracks* kotherTracks;
  
  if (exist<LHCb::Tracks>(LHCb::TrackLocation::Default))
  {
    kotherTracks = get<LHCb::Tracks>(LHCb::TrackLocation::Default);
    LHCb::Tracks::iterator it;
    for (it = kotherTracks->begin(); it != kotherTracks->end(); ++it){
      LHCb::Track* test = *it;
      if (test->type()==3) longTracks.push_back(*it);
    }
  }
  else info() << "There are no Tracks in the default location" << endmsg;
  
  //### Loop over DaPlus and DaMinus ###// 
  for (LHCb::Particle::ConstVector::const_iterator imp = DaPlus.begin() ;
       imp != DaPlus.end(); ++imp )
  {
    const LHCb::Particle* daPlus = *imp;
    
    for (LHCb::Particle::ConstVector::const_iterator imm = DaMinus.begin();
         imm != DaMinus.end(); ++imm)
    {
      const LHCb::Particle* daMinus = *imm;

      //### Find and plot diparticle mass of reconstructed particles ###// 
      Gaudi::LorentzVector twoDa = daPlus->momentum() + daMinus->momentum();
      debug() << "Rec two daughter mass is " << twoDa.M()/GeV << " GeV" << endmsg;
      motherTuple->column("rec_DiMuon_InvMass",    twoDa.M());
      motherTuple->column("numMuonsPerEvent",      (unsigned long long)nDaughters);
      motherTuple->column("numLongTracksPerEvent", (unsigned long long)longTracks.size());
      
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
      /*
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
          StatusCode scIPS = m_extra->ImpactParameterSum(daPlus, daMinus, *ipv, 
          //                                               tmpIPplus, tmpIPminus, tmpIPEplus, tmpIPEminus, tmpIPtot, tmpIPEtot);
          
          //StatusCode scPlus  = m_extra->SignedImpactParameter(daPlus, (*ipv), tmpIPplus,  tmpIPEplus);
          //StatusCode scMinus = m_extra->SignedImpactParameter(daMinus, (*ipv), tmpIPminus, tmpIPEminus);
          
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
        twoDaughters.push_back(muTrack1);
        twoDaughters.push_back(muTrack2);
        m_pvtool = tool<IPVOfflineTool>("PVOfflineTool");
        StatusCode scRefit = m_pvtool->reDoSinglePV(prims[pvID-1]->position(), twoDaughters, refittedPV);
        if (!scRefit) refitComplete = 0;
        
        //### Store coords of refitted vertex ###//
        refitPV_x = refittedPV.position().x();
        refitPV_y = refittedPV.position().y();
        refitPV_z = refittedPV.position().z();
                
        //### calculate new IP for each daughter based on the refitted PV ###//
        //m_extra->ImpactParameterSum(daPlus, daMinus, &refittedPV, refitIPplus, refitIPminus, 
        //                         refitIPEplus, refitIPEminus, refitIPtot, refitIPEtot);
        //StatusCode scRefitPlus  = m_extra->SignedImpactParameter(daPlus, &refittedPV, refitIPplus,  refitIPEplus);
        //StatusCode scRefitMinus = m_extra->SignedImpactParameter(daMinus, &refittedPV, refitIPminus, refitIPEminus);
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
      motherTuple->column("rec_DiMuon_recIP_whichPV",    pvID);
      motherTuple->column("rec_DiMuon_recPV_X",          recPV_x);
      motherTuple->column("rec_DiMuon_recPV_Y",          recPV_y);
      motherTuple->column("rec_DiMuon_recPV_Z",          recPV_z);
      motherTuple->column("rec_DiMuon_recIP_Total",      recIPtot);
      motherTuple->column("rec_DiMuon_recIPE_Total",     recIPEtot);
      motherTuple->column("rec_DiMuon_recIP_MuPlus",     recIPplus);
      motherTuple->column("rec_DiMuon_recIP_MuMinus",    recIPminus);
      motherTuple->column("rec_DiMuon_recIPE_MuPlus",    recIPEplus);
      motherTuple->column("rec_DiMuon_recIPE_MuMinus",   recIPEminus);
      motherTuple->column("rec_DiMuon_refitComplete",    refitComplete); //Refit can fail while still returning a new PV 
      motherTuple->column("rec_DiMuon_refitPV_X",        refitPV_x);
      motherTuple->column("rec_DiMuon_refitPV_Y",        refitPV_y);
      motherTuple->column("rec_DiMuon_refitPV_Z",        refitPV_z);
      motherTuple->column("rec_DiMuon_refitIP_Total",    refitIPtot);
      motherTuple->column("rec_DiMuon_refitIP_MuPlus",   refitIPplus);
      motherTuple->column("rec_DiMuon_refitIP_MuMinus",  refitIPminus);
      motherTuple->column("rec_DiMuon_refitIPE_Total",   refitIPEtot);
      motherTuple->column("rec_DiMuon_refitIPE_MuPlus",  refitIPEplus);
      motherTuple->column("rec_DiMuon_refitIPE_MuMinus", refitIPEminus);
      
      //### Find Distance of Closest Approach for DiMuons ###//
      double doca=m_errorCode*mm, docaChi2=m_errorCode*mm;
      sc = distTool->distance(*imp, *imm, doca, docaChi2);
      motherTuple->column("rec_DiMuon_DOCA",      doca);
      motherTuple->column("rec_DiMuon_DOCA_Chi2", docaChi2);
      */
      //### Now make the vertex by calling the Vertex Fitter (returns vertex and mother particle) ###//
      LHCb::Vertex DaDaVertex;
      LHCb::Particle Mother(motherID); // Fixed address in stack of mother to be created - Will be replaced by pointer to heap
      IVertexFit* testTool = tool<IVertexFit>("LoKi::VertexFitter"); // TEST - Tool to fit vertices
      StatusCode scFit = testTool->fit(*(*imp),*(*imm),DaDaVertex,Mother); // Seems to work but need to verify
      //StatusCode scFit = vertexFitter()->fit(*(*imp),*(*imm),DaDaVertex,Mother); // Old method in Z2TauTau implementation
           
      if (!scFit)
      {
        //### Output error code and plot any remaining data ###// 
        Warning("Fit error").ignore();
        err() << "Error while fitting" << endmsg;
        motherTuple->column("rec_DiMuon_Chi2",              -m_errorCode*mm);
        motherTuple->column("rec_FittedVertex_IP_Total",     m_errorCode*mm);
        motherTuple->column("rec_FittedVertex_IPE_Total",    m_errorCode*mm);
        motherTuple->column("rec_FittedVertex_IP_MuPlus",    m_errorCode*mm);
        motherTuple->column("rec_FittedVertex_IPE_MuPlus",   m_errorCode*mm);
        motherTuple->column("rec_FittedVertex_IP_MuMinus",   m_errorCode*mm);
        motherTuple->column("rec_FittedVertex_IPE_MuMinus",  m_errorCode*mm);
        
        sc = Da->plotDaughter(*imp, prims, motherTuple, hitDistTuple,
                              "_MuPlus_", runMC);
        if(!sc) return mothers;
        else 
        {
          counter("plotDaughter")++;
          sc = Da->plotDaughter(*imm, prims, motherTuple, hitDistTuple,
                                "_MuMinus_", runMC);
          if (sc) counter("plotDaughter")++;
          else return mothers;
        }
        motherTuple->write();
        counter ("FitError")++;
        continue;
      }
      
      if (msgLevel(MSG::DEBUG)) debug() << "Vertex fit at " << DaDaVertex.position()/cm
                                        << " with chi^2 " << DaDaVertex.chi2() << endmsg;
      motherTuple->column("rec_DiMuon_Chi2", DaDaVertex.chi2());
      
      //### Find IP for dimuon with the fitted vertex  ###//
      double fitIPplus=m_errorCode*mm, fitIPminus=m_errorCode*mm, fitIPEplus=m_errorCode*mm, fitIPEminus=m_errorCode*mm;
      double fitIPtot=m_errorCode*mm, fitIPEtot = m_errorCode*mm;

      //m_extra->ImpactParameterSum(daPlus, daMinus, &DaDaVertex, fitIPplus, fitIPminus,
      //                          fitIPEplus, refitIPminus, fitIPtot, fitIPEtot);
      //sc = m_extra->SignedImpactParameter(*imp, &DaDaVertex, fitIPplus,  fitIPEplus);  // Calculate for muPlus
      //sc = m_extra->SignedImpactParameter(*imm, &DaDaVertex, fitIPminus, fitIPEminus); // Calculate for muMinus 

      motherTuple->column("rec_FittedVertex_IP_Total",    fitIPtot);
      motherTuple->column("rec_FittedVertex_IPE_Total",   fitIPEtot);
      motherTuple->column("rec_FittedVertex_IP_MuPlus",   fitIPplus);
      motherTuple->column("rec_FittedVertex_IPE_MuPlus",  fitIPEplus);
      motherTuple->column("rec_FittedVertex_IP_MuMinus",  fitIPminus);
      motherTuple->column("rec_FittedVertex_IPE_MuMinus", fitIPEminus);
      
      //### Mandatory. Set to true if event is accepted. ###//
      setFilterPassed(true);
      
      //### Declare the mother to the PhysDesktop ###//
      mothers.push_back( Mother.clone() );

      if (msgLevel(MSG::DEBUG)) debug() << "Saved mother " << Mother.particleID().pid()
                                        << " to desktop" << endmsg;
      
      sc =  Da->plotDaughter(*imp, prims, motherTuple, hitDistTuple,
                             "_MuPlus_", runMC);
      if (!sc) return mothers;
      else 
      {
        counter("plotDaughterTest")++;
        sc = Da->plotDaughter(*imm, prims, motherTuple, hitDistTuple,
                              "_MuMinus_", runMC);
        if (sc) counter("plotDaughterTest")++;
        else return mothers;
      }
      
      motherTuple->write();
      counter ("Mothers")++; //Booked just like plots, counts # of candidates 
    }
  }

  //### Store number of candidates per event in a seperate tuple ###// 
  debug() << "Event Number: " << counter("EventNumber").nEntries() << endmsg;
  debug() << "Number of candidates in this event: " << numCandidates << endmsg;
  candTuple->column("EventNumber", (unsigned long long)counter("EventNumber").nEntries());
  candTuple->column("numCandidates", numCandidates);
  candTuple->write();

  if(!sc) return mothers;
  
  delete Da;                        // Deallocate memory from the heap for LoopOnDaughters class 
  
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
