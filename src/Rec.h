// $Id: Rec.h, 12/04/2011 rhuston Exp $
#ifndef REC_H 
#define REC_H 1 

//#####################//
//### INCLUDE FILES ###//
//#####################//

// Local
#ifndef EXTRA_H
#include "Extra.h"
#endif

//### from DaVinci, this is a specialized GaudiAlgorithm ###//
//#include "Kernel/DVAlgorithm.h"
#include "Kernel/DaVinciTupleAlgorithm.h"
//#include "Kernel/DaVinciAlgorithm.h"
//### MC data retrieval ###//
#include "Event/MCParticle.h"
//### Is MC Reconstructible? ###//
#include "MCInterfaces/IMCReconstructible.h"
//### Extract L0 Decision ###//
#include "Event/L0DUReport.h"
//### Extract HLT Decision ###//
//#include "Event/HltDecReports.h"
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

//class Rec : public DVAlgorithm 
class Rec : public DaVinciTupleAlgorithm 
//class Rec : public DaVinciAlgorithm 
{
public: 
  /// Standard constructor
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
  const std::string& m_strExtra = "Extra";
  Extra* m_extra;                        //< Utility for extra analysis
  static int constructorCalls;
  static int destructorCalls;
};
#endif // REC_H
