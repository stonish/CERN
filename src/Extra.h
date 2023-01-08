// $Id: Extra.h, 12/04/2011 rhuston Exp $
#ifndef EXTRA_H
#define EXTRA_H 1

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

/*
 * @class Extra Extra.h
 *
 *  Extra algorithms used when real data was first available
 *  Generally the algorithms are not run anymore but are kept in case they may be needed
 *  Algorithmsa focus on muon misidentification and residuals
 *
 *  @author Shane Huston
 *  @date   12-04-2011
 */

class Extra : public DaVinciTupleAlgorithm
{
public:
  Extra( const std::string& name, ISvcLocator* pSvcLocator );
  virtual ~Extra( );                   ///< Destructor
  virtual StatusCode initialize();     ///< Algorithm initialization
  virtual StatusCode execute   ();     ///< Algorithm execution
  virtual StatusCode finalize  ();     ///< Algorithm finalization

  StatusCode ImpactParameterSum(const LHCb::Particle*, const LHCb::Particle*, const LHCb::VertexBase*,
                   double&, double&, double&, double&, double&, double&);
  StatusCode SignedImpactParameter(const LHCb::Particle*, const LHCb::VertexBase*, double&, double&);
  StatusCode SignedImpactParameter(const LHCb::Particle*, const Gaudi::XYZPoint&, double&, double&);
  StatusCode muonChamberHits(Tuple);                              // Evaluate hits in the muon chambers
  StatusCode fakeMuon(const LHCb::Particle::ConstVector&, Tuple); // Muon Misidentification (Now replaced with GaudiPython code)

protected:

private:
  int m_errorCode;                       //< Error code for nTuples
  //int m_coneMax;                         //< Defines maximum cone size around muon for isolation criteria
  //int m_numDivs;                         //< Defines the number of differen cones to analyse
  //int m_muHits;                          //< Flag to run on hits in the muon chambers and calculate distances etc
  //double m_pi;                           //< Defines pi for calcultion of azimuthal angle difference
  IDistanceCalculator *m_distTool;       //< For IP calculation
  ISvcLocator *m_local;                  //< Test to use pSvcLocator to call external local algorithms 
};
#endif // EXTRA_H
