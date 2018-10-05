#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalMegaClustering.h"

// Geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
//
#include "DataFormats/CaloRecHit/interface/CaloID.h"


namespace HGCal_helpers {

class coordinates {
 public:
  coordinates() : x(0), y(0), z(0), eta(100), phi(0) {}
  float x, y, z, eta, phi;
  inline math::XYZTLorentzVectorD toVector() { return math::XYZTLorentzVectorD(x, y, z, 0); }
};
class simpleTrackPropagator {
 public:
  simpleTrackPropagator(MagneticField const *f)
      : field_(f), prod_(field_, alongMomentum, 5.e-5), absz_target_(0) {
    ROOT::Math::SMatrixIdentity id;
    AlgebraicSymMatrix55 C(id);
    C *= 0.001;
    err_ = CurvilinearTrajectoryError(C);
  }
  void setPropagationTargetZ(const float &z);

  bool propagate(const double px, const double py, const double pz, const double x, const double y,
                 const double z, const float charge, coordinates &coords) const;

  bool propagate(const math::XYZTLorentzVectorD &momentum, const math::XYZTLorentzVectorD &position,
                 const float charge, coordinates &coords) const;

 private:
  simpleTrackPropagator() : field_(0), prod_(field_, alongMomentum, 5.e-5), absz_target_(0) {}
  const RKPropagatorInS &RKProp() const { return prod_.propagator; }
  Plane::PlanePointer targetPlaneForward_, targetPlaneBackward_;
  MagneticField const *field_;
  CurvilinearTrajectoryError err_;
  defaultRKPropagator::Product prod_;
  float absz_target_;
};

void simpleTrackPropagator::setPropagationTargetZ(const float &z) {
  targetPlaneForward_ = Plane::build(Plane::PositionType(0, 0, std::abs(z)), Plane::RotationType());
  targetPlaneBackward_ =
      Plane::build(Plane::PositionType(0, 0, -std::abs(z)), Plane::RotationType());
  absz_target_ = std::abs(z);
}
bool simpleTrackPropagator::propagate(const double px, const double py, const double pz,
                                      const double x, const double y, const double z,
                                      const float charge, coordinates &output) const {
  output = coordinates();

  typedef TrajectoryStateOnSurface TSOS;
  GlobalPoint startingPosition(x, y, z);
  GlobalVector startingMomentum(px, py, pz);
  Plane::PlanePointer startingPlane =
      Plane::build(Plane::PositionType(x, y, z), Plane::RotationType());
  TSOS startingStateP(
      GlobalTrajectoryParameters(startingPosition, startingMomentum, charge, field_), err_,
      *startingPlane);

  TSOS trackStateP;
  if (pz > 0) {
    trackStateP = RKProp().propagate(startingStateP, *targetPlaneForward_);
  } else {
    trackStateP = RKProp().propagate(startingStateP, *targetPlaneBackward_);
  }
  if (trackStateP.isValid()) {
    output.x = trackStateP.globalPosition().x();
    output.y = trackStateP.globalPosition().y();
    output.z = trackStateP.globalPosition().z();
    output.phi = trackStateP.globalPosition().phi();
    output.eta = trackStateP.globalPosition().eta();
    return true;
  }
  return false;
}

bool simpleTrackPropagator::propagate(const math::XYZTLorentzVectorD &momentum,
                                      const math::XYZTLorentzVectorD &position, const float charge,
                                      coordinates &output) const {
  return propagate(momentum.px(), momentum.py(), momentum.pz(), position.x(), position.y(),
                   position.z(), charge, output);
}

}  // HGCal_helpers


HGCalMegaClustering::HGCalMegaClustering(const edm::ParameterSet& conf, std::string gun_type, float GEN_engpt, int pidSelected, int energyRadius, int frontRadius, int backRadius, bool doPileupSubtraction)
{
	gun_type_ = gun_type; 
	GEN_engpt_=GEN_engpt; 
	pidSelected_=pidSelected; 
	energyRadius_=energyRadius; 
	frontRadius_=frontRadius; 
	backRadius_=backRadius;
	doPileupSubtraction_=doPileupSubtraction;
	mySimEvent_ = new FSimEvent(conf);
}

void HGCalMegaClustering::retrieveLayerPositions(unsigned layers) {

  for (unsigned ilayer = 1; ilayer <= layers; ++ilayer) {
    const GlobalPoint pos = rhtools_.getPositionLayer(ilayer);
    layerPositions_.push_back(pos.z());
  }
}

float HGCalMegaClustering::getConeRadius(float frontRadius, float  backRadius, float  z, float maxval)
{
    float depthTerm = backRadius * (fabs(z)-320.7)/(407.8-320.7);
    float val = frontRadius + depthTerm;
    if (val > maxval)
    {
        return maxval;
    }
    else
    {
        return val;
    }

}

void HGCalMegaClustering::getEventSetup(const edm::EventSetup& es)
{
	edm::ESHandle<HepPDT::ParticleDataTable> pdt;
  	es.getData(pdt);
  	mySimEvent_->initializePdt(&(*pdt));
  	rhtools_.getEventSetup(es);
  	retrieveLayerPositions(52);
}

void HGCalMegaClustering::getMegaClusters(

  	const std::vector<SimTrack> & simTracks_,
  	const std::vector<SimVertex> & simVertices_,
	//const std::vector<reco::GenParticle> & genParticles_,
	const std::vector<reco::HGCalMultiCluster> & multiClusters_,
	const HGCRecHitCollection & recHitsEE_,
	const HGCRecHitCollection & recHitsFH_,
	const HGCRecHitCollection & recHitsBH_,
	const reco::CaloClusterCollection & clusters_)
{


  	mySimEvent_->fill(simTracks_, simVertices_);

	std::vector<FSimTrack *> selectedGen;
	for (unsigned int i = 0; i < mySimEvent_->nTracks()	; ++i) 
	{
		int reachedEE = 0;
		FSimTrack &myTrack(mySimEvent_->track(i));
		if (std::abs(myTrack.vertex().position().z()) >= layerPositions_[0]) continue;
		unsigned nlayers = rhtools_.lastLayerFH();
		if (myTrack.noEndVertex())  // || myTrack.genpartIndex()>=0)
		{



	}




		
  //     HGCal_helpers::coordinates propcoords;
  //     bool reachesHGCal = toHGCalPropagator.propagate(
  //         myTrack.momentum(), myTrack.vertex().position(), myTrack.charge(), propcoords);
  //     vtx = propcoords.toVector();

  //     if (reachesHGCal && vtx.Rho() < hgcalOuterRadius_ && vtx.Rho() > hgcalInnerRadius_) {
  //       reachedEE = 2;
  //       double dpt = 0;

  //       for (int i = 0; i < myTrack.nDaughters(); ++i) dpt += myTrack.daughter(i).momentum().pt();
  //       if (abs(myTrack.type()) == 11) fbrem = dpt / myTrack.momentum().pt();
  //     } else if (reachesHGCal && vtx.Rho() > hgcalOuterRadius_)
  //       reachedEE = 1;

  //     HGCal_helpers::simpleTrackPropagator indiv_particleProp(aField_);
  //     for (unsigned il = 0; il < nlayers; ++il) {
  //       const float charge = myTrack.charge();
  //       indiv_particleProp.setPropagationTargetZ(layerPositions_[il]);
  //       HGCal_helpers::coordinates propCoords;
  //       indiv_particleProp.propagate(myTrack.momentum(), myTrack.vertex().position(), charge,
  //                                    propCoords);

  //       xp.push_back(propCoords.x);
  //       yp.push_back(propCoords.y);
  //       zp.push_back(propCoords.z);
  //     }
  //   } else {
  //     vtx = myTrack.endVertex().position();
  //   }












	// }

	// unsigned int npart = mySimEvent_->nTracks();
 //  for (unsigned int i = 0; i < npart; ++i) {


	// selectedGen = genParticles[(abs(genParticles.pid) == pidSelected) & (genParticles.reachedEE > 0)]
 //    if gun_type == "pt":
 //        selectedGen = selectedGen[(selectedGen.pt >= GEN_engpt*.999)]
 //    else:
 //        selectedGen = selectedGen[(selectedGen.energy >= GEN_engpt*.999)]
}

// void HGCalMegaClustering::populate()
// {


// }

// void HGCalMegaClustering::makeClusters()
// {


// }

// void HGCalMegaClustering::getClusters()
// {


// }

// void HGCalMegaClustering::getEvent()
// {


// }

// void HGCalMegaClustering::getEventSetup()
// {


// }

// void HGCalMegaClustering::reset()
// {


// }
