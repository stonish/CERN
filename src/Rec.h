// $Id: Rec.h, 12/04/2011 rhuston Exp $
#ifndef REC_H 
#define REC_H 1 

//#####################//
//### INCLUDE FILES ###//
//#####################//

// Local
#ifndef EXTRA_H // TODO: Assuming we want to include Extra.h, should this be #ifdef instead?
#include "Extra.h" // TODO: It doesn't lok like we're using Extra here. Can this be removed? (Extra is referenced in some commented out code)
#endif

//### from DaVinci, this is a specialized GaudiAlgorithm ###//
#include "Kernel/DaVinciTupleAlgorithm.h"
//### MC data retrieval ###//
#include "Event/MCParticle.h"
//### Is MC Reconstructible? ###//
#include "MCInterfaces/IMCReconstructible.h"
//### Extract L0 Decision ###//
#include "Event/L0DUReport.h"
//### Association between MC and Reconstructed###//
#include "Kernel/Particle2MCLinker.h"
// Convert MCParticle to Particle??
#include "Kernel/MCParticleMakerBase.h"
#include "Kernel/IParticleMaker.h"
#include "GaudiAlg/GaudiTool.h"
//### Tool to refit a PV excluding specified tracks ###//
#include "TrackInterfaces/IPVOfflineTool.h"

/*
 * @class Rec Rec.h
 *  
 *  Master algorithms for handling reconstructed data
 *  Will be used to call LoopOnDaughters and MakeMother
 *  
 *  @author Shane Huston
 *  @date   12-04-2011
 */

class Rec : public DaVinciTupleAlgorithm 
{
public: 
  Rec( const std::string& name, ISvcLocator* pSvcLocator );
  virtual ~Rec( );                   ///< Destructor
  virtual StatusCode initialize();   ///< Algorithm initialization
  virtual StatusCode execute   ();   ///< Algorithm execution
  virtual StatusCode finalize  ();   ///< Algorithm finalization

  // Extract daughters and PVs, run LoopOnDaughters and MakeMother
  LHCb::Particle::Vector loopOnRec(const LHCb::Particle::ConstVector&, const LHCb::RecVertex::Range, 
                                   LHCb::ParticleID, Tuple, Tuple, Tuple, Tuple, int); 
  // Build mother (Trying to recreate Z0 with missing energy due to neutrinos)
  LHCb::Particle::Vector makeMother(const LHCb::Particle::ConstVector& daughters, const LHCb::RecVertex::Range, LHCb::ParticleID, 
                                         Tuple, Tuple, Tuple, int) ;
  
protected:

private:
  int m_errorCode;                       //< Error code for nTuples
  ISvcLocator *m_local;                  //< Test to use pSvcLocator to call external local algorithms 
  IPVOfflineTool* m_pvtool;              //< Tool to refit a PV excluding specified tracks

  void getAndStoreEventNumber(Tuple tuple);
  void getAndStoreRunNumberAndL0EventID(Tuple tuple);
  void storeImpactParameterData(double fitIPplus, double fitIPminus, double fitIPEplus, double fitIPEminus,
                                double fitIPtot, double fitIPEtot, Tuple tuple);
  void calculateImpactParametersWithReconstructedPrimaryVertices(const LHCb::RecVertex::Range prims,
                                const LHCb::Particle* muPlus, const LHCb::Particle* muMinus,
                                const LHCb::Track* muPlusTrack, const LHCb::Track* muMinusTrack, Tuple tuple);
  void getAndStoreDiMuonDistanceOfClosestApproach(const LHCb::Particle* muPlus, const LHCb::Particle* muMinus, Tuple tuple);
  void fitVertexAndStoreImpactParameterData(LHCb::ParticleID motherID, const LHCb::Particle* muPlus,
                                            const LHCb::Particle* muMinus, Tuple motherTuple, Tuple hitDistTuple);
  bool plotDaughters(const std::string& counterName, const LHCb::Particle* muPlus, const LHCb::Particle* muMinus,
                     const LHCb::RecVertex::Range prims, Tuple motherTuple, Tuple hitDistTuple, int runMC);
  void storeAndWriteNumberOfCandidatesPerEvent(int numCandidates, Tuple tuple);
};
#endif // REC_H
