// $Id: Z2MuMu.h, 07/07/2009 rhuston Exp $
#ifndef Z2MUMU_H 
#define Z2MUMU_H 1 

//#####################//
//### INCLUDE FILES ###//
//#####################//

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

/** @class Z2MuMu Z2MuMu.h
 *  
 *  Simple algorithm that loops on muons and stores some quantities of interest in ntuples. 
 *  (Includes MC truth, Particle2MCLinker, Reconstructibility, and L0 / HLT Decisions)
 *  (No cuts applied, apply cuts to tuple afterwards)
 *  (Output error code to tuple when needed)
 *  (Seperate LoopOnDaughters and plotDaughter into seperate algorithm)
 *
 *  @author Shane Huston
 *  @date   11-08-2009
 */

class Z2MuMu : public DaVinciTupleAlgorithm {
public: 
  Z2MuMu( const std::string& name, ISvcLocator* pSvcLocator );
  virtual ~Z2MuMu( );                  ///< Destructor
  virtual StatusCode initialize();     ///< Algorithm initialization
  virtual StatusCode execute   ();     ///< Algorithm execution
  virtual StatusCode finalize  ();     ///< Algorithm finalization

protected:

private:
  // Look at MC data
  StatusCode loopOnMC();
  // Filter MC particles by type and charge (supplying all possible vectors to be filled)
  StatusCode filterMC(LHCb::MCParticle*, const LHCb::MCParticle*, LHCb::MCParticle::Vector&, LHCb::MCParticle::Vector&);
  
  int m_runMC;                           //< Flag to run on MC data 
  int m_errorCode;                       //< Error code for nTuples 
  double m_motherMassWin;                //< Mass window
  double m_motherChi2;                   //< Max mother chi^2
  LHCb::ParticleID m_motherID;           //< Mother ID
  double m_motherMass;                   //< Mother mass
  std::string m_motherName;              //< Mother name
  LHCb::MCParticle::Vector m_mcparts;    //< Vector to store all MC Particles in event
  const LHCb::L0DUReport *m_L0Rep;       //< Get L0 Report
  IMCReconstructible* m_recoTool;        //< Is MC particle reconstructible?  
  int m_coneMax;                         //< Defines cone around muon for isolation criteria
  int m_numDivs;                         //< Defines the number of different cones to analyse
  int m_muHits;                          //< Flag to run on hits in the muon chamber to determine distances etc
  double m_pi;                           //< Defines pi for calcultion of azimuthal angle difference
  ISvcLocator* m_local;                  //< Test to use pSvcLocator to call external local algorithms
  //### Test
  LHCb::MCParticle::Vector m_McSignalParts; //< Vector to store any MCPartilces belonging to the desired process 

};
#endif // Z2MUMU_H
