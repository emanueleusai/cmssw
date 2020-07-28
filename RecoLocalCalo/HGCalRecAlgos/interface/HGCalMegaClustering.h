#ifndef RecoLocalCalo_HGCalRecAlgos_HGCalMegaClustering_h
#define RecoLocalCalo_HGCalRecAlgos_HGCalMegaClustering_h

//Geometry includes
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/GeometrySurface/interface/PlaneBuilder.h"

//Calorimetry Data formats
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

//Simulation and event generator data formats
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FastSimulation/Event/interface/FSimEvent.h"
#include "FastSimulation/Event/interface/FSimTrack.h"
#include "FastSimulation/Event/interface/FSimVertex.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

//magnetic field includes
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeGeometry/interface/MagVolumeOutsideValidity.h"

//core CMSSW includes
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"

//math includes
#include "TVector3.h"
#include "DataFormats/Math/interface/Point3D.h"

//reconstruction tools
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "TrackPropagation/RungeKutta/interface/defaultRKPropagator.h"

// C/C++ headers
#include <string>
#include <vector>


class HGCalMegaClustering
{

public:


    HGCalMegaClustering(const edm::ParameterSet& conf, std::string gun_type, float GEN_engpt, int pidSelected, int energyRadius=6, int frontRadius=3, int backRadius=8, bool doPileupSubtraction=true);
    virtual ~HGCalMegaClustering() {}

//retrieve configuration
    void getEventSetup(const edm::EventSetup& es);

//main method, produces pairs (megacluster, particle momentum)
    std::vector<std::pair<reco::BasicCluster,float> > getMegaClusters(
        const std::vector<SimTrack> & simTracks_,
        const std::vector<SimVertex> & simVertices_,
        const std::vector<reco::HGCalMultiCluster> & multiClusters_,
        const HGCRecHitCollection & recHitsEE_,
        const HGCRecHitCollection & recHitsFH_,
        const HGCRecHitCollection & recHitsBH_,
        const reco::CaloClusterCollection & layerClusters_,
        const edm::HepMCProduct & hepmc_);

private:

//helper function to subtract energy from pileup interactions
    std::pair<float,float> pileupSubtraction(
        const std::pair<const FSimTrack*,const reco::HGCalMultiCluster*>& matchedMultiCluster,
        const std::vector<const reco::CaloCluster *>& selectedLayerClusters,
        int layer);

//detector parameters
    static constexpr unsigned int maxlayer = 52;
    static constexpr float hgcalOuterRadius_ = 160.;
    static constexpr float hgcalInnerRadius_ = 25.;
    const std::vector<float> energyWeights = { 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02,
                                               1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02,
                                               0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86,
                                               1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12, 1.12
                                             };

//cone-building helper functions
    float getConeRadius(float frontRadius, float  backRadius, float  z, float maxval=9999.);
    void retrieveLayerPositions(unsigned layers);

//configuration variables
    int pidSelected_;
    std::string gun_type_;
    float GEN_engpt_;
    float vz_;
    float energyRadius_;
    float frontRadius_;
    float backRadius_;
    bool doPileupSubtraction_;
    FSimEvent *mySimEvent_;

//variables for the helper tools
    hgcal::RecHitTools rhtools_;
    std::vector<float> layerPositions_;
    MagneticField const *aField_;
    std::map<DetId, const HGCRecHit *> hitmap_;


};

#endif