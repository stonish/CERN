// $Id: $
#ifndef ISOLATIONMEASUREMENT_H 
#define ISOLATIONMEASUREMENT_H 1

// Include files
// from DaVinci, this is a specialized GaudiAlgorithm
//#include "Kernel/DVAlgorithm.h"
#include "Kernel/DaVinciTupleAlgorithm.h"
//#include "Kernel/DaVinciAlgorithm.h"

/** @class IsolationMeasurement IsolationMeasurement.h
 *  
 *
 *  @author  Shane Huston
 *  @date   2012-05-13
 */
//class IsolationMeasurement : public DVAlgorithm {
class IsolationMeasurement : public DaVinciTupleAlgorithm {
//class IsolationMeasurement : public DaVinciAlgorithm {
public: 
  /// Standard constructor
  IsolationMeasurement( const std::string& name, ISvcLocator* pSvcLocator );
  //IsolationMeasurement(int numDivs, int coneMax, const std::string& name, ISvcLocator* pSvcLocator );
  
  virtual ~IsolationMeasurement( ); ///< Destructor

  virtual StatusCode initialize();  ///< Algorithm initialization
  virtual StatusCode execute   ();  ///< Algorithm execution
  virtual StatusCode finalize  ();  ///< Algorithm finalization

  void Setup(int numDivs, int coneMax);
  void UpdateMeasurement(double dist_PhiEta, LHCb::Track* trk);
  void UpdateMeasurement(double dist_PhiEta, LHCb::ProtoParticle* proto);
  void PlotMeasurement(Tuples::Tuple tuple, const std::string& head, const std::string& type);

protected:

private:
  ISvcLocator *m_local;

  int m_numDivs;
  std::vector<int> num_tot;
  std::vector<int> num_T;
  std::vector<int> num_Velo;
  std::vector<int> num_TT;
  std::vector<double> tot_meas;
  std::vector<double> tot_ghost;
  std::vector<double> tot_likelihood;
  std::vector<double> tot_EcalE;
  std::vector<double> tot_HcalE;
  std::vector<Gaudi::XYZVector> ConeMomentum;
  std::vector<double> coneSize;
};
#endif // ISOLATIONMEASUREMENT_H
