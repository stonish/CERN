// $Id: B2DMuNu.h, 14/08/2009 rhuston Exp $
#ifndef B2DMUNU_H 
#define B2DMUNU_H 1 

//#####################//
//### INCLUDE FILES ###//
//#####################//

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

/** @class B2DMuNu B2DMuNu.h
 *  
 *  Simple algorithm that loops on daughters and stores some quantities of interest in ntuples. 
 *  (Includes MC truth, Particle2MCLinker, Reconstructibility, and L0 / HLT Decisions)
 *  (No cuts applied, apply cuts to tuple afterwards)
 *  (Output error code to tuple when needed)
 *  (Seperate LoopOnDaughters and plotDaughter into seperate algorithm)
 *
 *  @author Shane Huston
 *  @date   14-08-2009
 */
//class B2DMuNu : public DVAlgorithm {
class B2DMuNu : public DaVinciTupleAlgorithm {
//class B2DMuNu : public DaVinciAlgorithm {
public: 
  /// Standard constructor
  B2DMuNu( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~B2DMuNu( );              ///< Destructor
  virtual StatusCode initialize();  ///< Algorithm initialization
  virtual StatusCode execute   ();  ///< Algorithm execution
  virtual StatusCode finalize  ();  ///< Algorithm finalization

protected:

private:
  // Look at MC data
  StatusCode loopOnMC();
  // Filter MC particles by type and charge (supplying all possible vectors to be filled)
  StatusCode filterMC(LHCb::MCParticle*, const LHCb::MCParticle*, LHCb::MCParticle::Vector&, LHCb::MCParticle::Vector&);
  // Find invariant mass of the MC Daughters
  //StatusCode invMassMCTest(LHCb::MCParticle::Vector&, LHCb::MCParticle::Vector&, Tuple, const std::string& head);
  // Plot some quantities for pure MC data
  //StatusCode plotMC(const LHCb::MCParticle*, const LHCb::MCParticle*, const std::string& head);
  //StatusCode plotMCTest(const LHCb::MCParticle*, const LHCb::MCParticle*, Tuple, const std::string& head);
  //Make a loop over reconstructed daughters   
  //StatusCode loopOnDaughters(const LHCb::Particle::ConstVector&)const;
  // Plot quantities for daughters (both reconstructed and associated MC)
  //StatusCode plotDaughter(const LHCb::Particle*, const std::string& head)const;
  //StatusCode plotDaughterTest(const LHCb::Particle*, Tuple, const std::string& head)const;
  // Build mother (Trying to recreate Z0 with missing energy due to neutrinos)
  //StatusCode makeMother(const LHCb::Particle::ConstVector&);

  int m_runMC;                          //< Flag to run on MC data
  int m_passEvent;                      //< Test MC to see if event should be analysed
  int m_errorCode;                      //< Error code for nTuples
  double m_motherMassWin;               //< Mass window
  double m_motherChi2;                  //< Max mother chi^2
  LHCb::ParticleID m_motherID;          //< Mother ID
  double m_motherMass;                  //< Mother mass
  std::string m_motherName;             //< Mother name
  LHCb::MCParticle::Vector m_mcparts;   //< Vector to store all MCParticles in event
  const LHCb::L0DUReport *m_L0Rep;      //< Get L0 Report
  //const LHCb::HltDecReports *m_HltRep;  //< Get HLT Report
  //const LHCb::HltDecReports *m_HltReps; //< TEST - New method to extract HLT decision
  IMCReconstructible* m_recoTool;       //< Is MC particle reconstructible?  
  int m_coneMax;                        //< Defines cone around muon for isolation criteria
  int m_numDivs;                        //< Defines the number of different cones to analyse
  int m_muHits;                         //< Flag to run on hits in the muon chamber to determine distances etc
  double m_pi;                          //< Defines pi for calcultion of azimuthal angle difference
  //Particle2MCLinker* m_assoc;           //< Get associated MC particle
  ISvcLocator* m_local;                 //< Test to use pSvcLocator to call external local algorithms
  //### Test
  LHCb::MCParticle::Vector m_McSignalParts; //< Vector to store any MCPartilces belonging to the desired process

};
#endif // B2DMUNU_H
