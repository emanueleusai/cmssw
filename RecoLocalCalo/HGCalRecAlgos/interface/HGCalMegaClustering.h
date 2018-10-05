#ifndef RecoLocalCalo_HGCalRecAlgos_HGCalMegaClustering_h
#define RecoLocalCalo_HGCalRecAlgos_HGCalMegaClustering_h

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

#include "FastSimulation/Event/interface/FSimEvent.h"
#include "FastSimulation/Event/interface/FSimTrack.h"
#include "FastSimulation/Event/interface/FSimVertex.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "TrackPropagation/RungeKutta/interface/defaultRKPropagator.h"
#include "DataFormats/GeometrySurface/interface/PlaneBuilder.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeGeometry/interface/MagVolumeOutsideValidity.h"

// C/C++ headers
#include <string>
#include <vector>





class HGCalMegaClustering
{


public:


HGCalMegaClustering(const edm::ParameterSet& conf, std::string gun_type, float GEN_engpt, int pidSelected, int energyRadius=6, int frontRadius=3, int backRadius=8, bool doPileupSubtraction=true);
virtual ~HGCalMegaClustering(){}

void getEventSetup(const edm::EventSetup& es);

void getMegaClusters(
	const std::vector<SimTrack> & simTracks_,
  	const std::vector<SimVertex> & simVertices_,
	const std::vector<reco::HGCalMultiCluster> & multiClusters_,
	const HGCRecHitCollection & recHitsEE_,
	const HGCRecHitCollection & recHitsFH_,
	const HGCRecHitCollection & recHitsBH_,
	const reco::CaloClusterCollection & clusters_);

// void populate();
// void makeClusters();
// void getClusters();
// void getEvent();
// void getEventSetup();
// void reset();






private:

static constexpr unsigned int maxlayer = 52;
static constexpr float hgcalOuterRadius_ = 160.;
static constexpr float hgcalInnerRadius_ = 25.;
const std::vector<float> energyWeights ={ 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 
                                          1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 
                                          0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 
                                          1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12};

float getConeRadius(float frontRadius, float  backRadius, float  z, float maxval=9999.);
void retrieveLayerPositions(unsigned layers);
int pidSelected_;
std::string gun_type_;
float GEN_engpt_;
int energyRadius_;
int frontRadius_;
int backRadius_;
bool doPileupSubtraction_;
FSimEvent *mySimEvent_;
hgcal::RecHitTools rhtools_;
std::vector<float> layerPositions_;
MagneticField const *aField_;

};

#endif
