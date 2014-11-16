// $Id: LoopOnDaughters.h, 30/09/2009 rhuston Exp $
#ifndef LOOPONDAUGHTERS_H 
#define LOOPONDAUGHTERS_H 1 

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
//#include "TrackInterfaces/ITrackExtrapolator.h"

/*
 * @class LoopOnDaughters LoopOnDaughters.h
 *  
 *  Simple algorithm that loops on supplied daughter particles and plots 
 *  some quantities of interest in ntuples. 
 *  (Output error code to tuple when needed)
 *  (Include isolation criteria on muons)
 *
 *  @author Shane Huston
 *  @date   30-09-2009
 */

//class LoopOnDaughters : public DVAlgorithm 
class LoopOnDaughters : public DaVinciTupleAlgorithm 
//class LoopOnDaughters : public DaVinciAlgorithm 
{
public: 
  /// Standard constructor
  LoopOnDaughters( const std::string& name, ISvcLocator* pSvcLocator );
  
  virtual ~LoopOnDaughters( );              ///< Destructor
  virtual StatusCode initialize();          ///< Algorithm initialization
  virtual StatusCode execute   ();          ///< Algorithm execution
  virtual StatusCode finalize  ();          ///< Algorithm finalization

  //Make a loop over reconstructed daughters   
  StatusCode loopOnDaughters(const LHCb::Particle::ConstVector&, const LHCb::RecVertex::Range prims, Tuple, Tuple, int)const;
  //StatusCode loopOnDaughters(const LHCb::Particle::ConstVector&, const LHCb::RecVertex::Container* prims, Tuple, int)const;
  // Plot quantities for daughters (both reconstructed and associated MC)
  StatusCode plotDaughter(const LHCb::Particle*, const LHCb::RecVertex::Range, Tuple, Tuple, const std::string&, int)const;
  //StatusCode plotDaughter(const LHCb::Particle*, const LHCb::RecVertex::Container*, Tuple, const std::string&, int)const;
  StatusCode plotReconstructedData(const LHCb::Particle* da, Tuple daTuple, const std::string& head) const;
  StatusCode plotLHCbIDs(LHCb::Track::LHCbIDContainer lhcbIDs, Tuple daTuple, const std::string& head) const;
  StatusCode plotGlobalEventCutVariables(Tuple daTuple, const std::string& head) const;
  StatusCode plotIPs(const LHCb::Particle* da, const LHCb::RecVertex::Range prims, Tuple daTuple, const std::string& head) const;
  StatusCode testIPs(const LHCb::Particle* da, const LHCb::RecVertex::Range prims, Tuple daTuple, const std::string& head) const;
  StatusCode trackIsolation(const LHCb::Particle* da, Tuple daTuple, const std::string& head, const bool) const;
  StatusCode protoParticleIsolation(const LHCb::Particle* da, Tuple daTuple, const std::string& head, const bool) const;
  // Find the distance from the track to a Muon chamber hit, store the shortest distance
  double getMinDist(const LHCb::Track*, double, std::vector<Gaudi::XYZPoint>, Tuple, const std::string&)const;
  LHCb::ProtoParticle::Vector getAllProtos()const;
  
protected:

private:
  
  int m_errorCode;                       //< Error code for nTuples 
  //const LHCb::L0DUReport *m_L0Rep;       //< Get L0 Report
  //const LHCb::HltDecReports *m_HltRep;   //< Get HLT Report
  int m_coneMax;                         //< Defines maximum cone size around muon for isolation criteria
  int m_numDivs;                         //< Defines the number of differen cones to analyse
  int m_muHits;                          //< Flag to run on hits in the muon chambers and calculate distances etc
  double m_pi;                           //< Defines pi for calcultion of azimuthal angle difference 
  ISvcLocator *m_local;                  //< Test to use pSvcLocator to call external local algorithms 
  //Extra* m_extra;                        //< Utility for extra analysis
  static int constructorCalls;    
  static int destructorCalls;  
};
#endif // LOOPONDAUGHTERS_H
