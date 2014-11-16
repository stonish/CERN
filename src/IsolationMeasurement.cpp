// $Id: IsolationMeasurement, 13/05/2012 rhuston Exp $

//#####################//
//### Include files ###//
//#####################// 

// from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "IsolationMeasurement.h"

using namespace Gaudi::Units;

//-----------------------------------------------------------------------------
// Implementation file for class : IsolationMeasurement
//
// 2012-05-13 : Robert Shane Huston
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( IsolationMeasurement )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
IsolationMeasurement::IsolationMeasurement(const std::string& name, ISvcLocator* pSvcLocator)
//: DVAlgorithm ( name , pSvcLocator )
  : DaVinciTupleAlgorithm ( name , pSvcLocator )
//: DaVinciAlgorithm ( name , pSvcLocator )
{
  m_local = pSvcLocator;
}
//=============================================================================
// Destructor
//=============================================================================
IsolationMeasurement::~IsolationMeasurement() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode IsolationMeasurement::initialize() 
{
  //StatusCode sc = DVAlgorithm::initialize(); 
  StatusCode sc = DaVinciTupleAlgorithm::initialize(); 
  //StatusCode sc = DaVinciAlgorithm::initialize(); 
  if ( sc.isFailure() ) return sc;

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode IsolationMeasurement::execute() 
{
  err() << "Execute phase of IsolationMeasurement class has been called." << endmsg;
  return StatusCode::SUCCESS;
}

//=============================================================================
// Setup cone size and number of divisions
//=============================================================================
void IsolationMeasurement::Setup(int numDivs, int coneMax)
{
  m_numDivs = numDivs;
  double divSize = double(coneMax)/m_numDivs;

  for (int i=0; i<m_numDivs; i++)
  {
    coneSize.push_back(coneMax - (i * divSize));

    num_tot.push_back(0);
    num_T.push_back(0);
    num_Velo.push_back(0);
    num_TT.push_back(0);

    tot_meas.push_back(0);
    tot_ghost.push_back(0);
    tot_likelihood.push_back(0);
    tot_EcalE.push_back(0);
    tot_HcalE.push_back(0);

    ConeMomentum.push_back(Gaudi::XYZVector(0, 0, 0));
  }
}

//=============================================================================
// Update measurement for track
//=============================================================================
void IsolationMeasurement::UpdateMeasurement(double dist_PhiEta, LHCb::Track* trk)
{
  for (int i=0; i<m_numDivs; i++)
  {
    if (dist_PhiEta < coneSize[i])
    {
      //### NB: hasT, hasVelo etc only check that track passes thru the station ###//
      //### Does not check that their are any hits on the track                 ###//
      //### To do this, would need to start by extracting trk->lhcbIDs()        ###//
      //### Then loop over these IDs and check what type of hit each one is!    ###//

      if (msgLevel(MSG::DEBUG))
      {
        // Debug does not work yet
        debug() << "Distance in Phi-Eta space from this track to the Muon is " << dist_PhiEta << endmsg;
        debug() << "Track has hit in T station?\t "    << trk->hasT() << endmsg;
        debug() << "Track has hit in Velo?\t "         << trk->hasVelo() << endmsg;
        debug() << "Track has hit in TT?\t "           << trk->hasTT() << endmsg;
        debug() << "Track has likelihood of: "         << trk->likelihood() << endmsg;
        debug() << "Track has ghostProbability of: "   << trk->ghostProbability() << endmsg;
        debug() << "Distance in Phi-Eta space from this track to the Muon is " << dist_PhiEta << endmsg;
      }

      num_tot[i]++;
      num_T[i]           += trk->hasT();
      num_Velo[i]        += trk->hasVelo();
      num_TT[i]          += trk->hasTT();
      tot_meas[i]        += trk->nMeasurements();
      if (trk->ghostProbability()!=999)
        tot_ghost[i]     += trk->ghostProbability();
      tot_likelihood[i]  += trk->likelihood();

      ConeMomentum[i]    += trk->momentum();
    }
  }
}

//=============================================================================
// Update measurement for ProtoParticle
//=============================================================================
void IsolationMeasurement::UpdateMeasurement(double dist_PhiEta, LHCb::ProtoParticle* proto)
{
  const LHCb::Track* trk = proto->track();

  for (int i=0; i<m_numDivs; i++)
  {
    if (dist_PhiEta < coneSize[i])
    {
      //### NB: hasT, hasVelo etc only check that track passes thru the station ###//
      //### Does not check that their are any hits on the track                 ###//
      //### To do this, would need to start by extracting trk->lhcbIDs()        ###//
      //### Then loop over these IDs and check what type of hit each one is!    ###//
      if (msgLevel(MSG::DEBUG))
      {
        debug() << "Distance in Phi-Eta space from this track to the Muon is " << dist_PhiEta << endmsg;
        debug() << "Track has hit in T station?\t "       << trk->hasT() << endmsg;
        debug() << "Track has hit in Velo?\t "            << trk->hasVelo() << endmsg;
        debug() << "Track has hit in TT?\t "              << trk->hasTT() << endmsg;
        debug() << "Track has likelihood of: "            << trk->likelihood() << endmsg;
        debug() << "Track has ghostProbability of: "      << trk->ghostProbability() << endmsg;
        debug() << "Ecal Energy for this protoParticle: " << proto->info( proto->CaloEcalE, 0) << endmsg;
        debug() << "Hcal Energy for this protoParticle: " << proto->info( proto->CaloHcalE, 0) << endmsg;
      }
        
      num_tot[i]++;
      num_T[i]          += trk->hasT();
      num_Velo[i]       += trk->hasVelo();
      num_TT[i]         += trk->hasTT();
      tot_meas[i]       += trk->nMeasurements();
      //if (trk->ghostProbability()!=999) // Applied to tracks but not protos. Prob a mistake, maybe fix after comparison
      tot_ghost[i]      += trk->ghostProbability();
      tot_likelihood[i] += trk->likelihood();
      tot_EcalE[i]      += proto->info( proto->CaloEcalE, 0);
      tot_HcalE[i]      += proto->info( proto->CaloHcalE, 0);

      ConeMomentum[i]   += trk->momentum();
    }
  }
}

//=============================================================================
// Plot results
//=============================================================================
void IsolationMeasurement::PlotMeasurement(Tuples::Tuple tuple, const std::string& head, const std::string& type)
{
  //### Plot data for tracks inside cones ###//
  char num [] = "123456789";

  for (int i=0; i<m_numDivs; i++)
  {
    tuple->column("rec" + head + type + "_" + num[i] + "_ConeSize",   coneSize[i]);
    tuple->column("rec" + head + type + "_" + num[i] + "_Num",        num_tot[i]);
    tuple->column("rec" + head + type + "_" + num[i] + "_hasT",       num_T[i]);
    tuple->column("rec" + head + type + "_" + num[i] + "_hasVelo",    num_Velo[i]);
    tuple->column("rec" + head + type + "_" + num[i] + "_hasTT",      num_TT[i]);
    tuple->column("rec" + head + type + "_" + num[i] + "_NumMeasure", tot_meas[i]);
    tuple->column("rec" + head + type + "_" + num[i] + "_Likelihood", tot_likelihood[i]);
    tuple->column("rec" + head + type + "_" + num[i] + "_GhostProb",  tot_ghost[i]);
    tuple->column("rec" + head + type + "_" + num[i] + "_P",          ConeMomentum[i].R()); //*MeV??
    tuple->column("rec" + head + type + "_" + num[i] + "_Pt",         ConeMomentum[i].R()*sin(ConeMomentum[i].theta())); //*MeV??
    tuple->column("rec" + head + type + "_" + num[i] + "_Phi",        ConeMomentum[i].phi());
    tuple->column("rec" + head + type + "_" + num[i] + "_Eta",        ConeMomentum[i].Eta());
    
    if (type == "Protos")
    {
      tuple->column("rec" + head + "Protos_" + num[i] + "_EcalE",   tot_EcalE[i]);
      tuple->column("rec" + head + "Protos_" + num[i] + "_HcalE",   tot_HcalE[i]);
    }
  }
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode IsolationMeasurement::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  //return DVAlgorithm::finalize();
  return DaVinciTupleAlgorithm::finalize();
  //return DaVinciAlgorithm::finalize();
}

//=============================================================================
