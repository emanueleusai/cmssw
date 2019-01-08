#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalMegaClustering.h"

// Geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
//
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/Math/interface/deltaR.h"


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
                                      const math::XYZTLorentzVectorD &position, const float &charge,
                                      coordinates &output) const {
  return propagate(momentum.px(), momentum.py(), momentum.pz(), position.x(), position.y(),
                   position.z(), charge, output);
}


std::vector<std::pair<FSimTrack*,reco::HGCalMultiCluster*> > getHighestEnergyObjectIndex(const std::vector<FSimTrack *> & reference, const std::vector<reco::HGCalMultiCluster> & object, const float & deltaR=0.1)
{
  std::vector<std::pair<FSimTrack*,reco::HGCalMultiCluster*> > pairings;
  for (int i=0;i<reference.size(),i++)
  {
    int index = -1;
    float maxenergy = -1.0;
    for(int j=0;j<object.size(),j++)
    {
      if (object[j]->energy()>maxenergy &&
          reco::deltaR(reference[i],object[j])<deltaR)
      {
        index=j;
        maxenergy=object[j]->energy();
      }
    }
    if (index>-1)
    {
      pairings.push_back(std::pair<FSimTrack*,reco::HGCalMultiCluster*>(reference[i],object[index]));
    }
  }
  return pairings;
}

TVector2 convertToXY(const float & eta, const float & phi, const float & z)
{
    t = exp(-1. * eta);
    x = z * 2. * t * cos(phi) / (1. - t*t);
    y = z * 2. * t * sin(phi) / (1. - t*t);
    return TVector2(x,y);
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

  layerPositions_.clear();
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
  	edm::ESHandle<MagneticField> magfield;
  	es.get<IdealMagneticFieldRecord>().get(magfield);
  	aField_ = &(*magfield);
}


std::pair<float,float> HGCalMegaClustering::pileupSubtraction(
  const std::pair<FSimTrack*,reco::HGCalMultiCluster*>& matchedMultiCluster, 
  const std::vector<reco::CaloCluster *>& selectedLayerClusters, 
  int layer)
{
  float pTSum = 0;
  float energySum = 0;
  auto layer_z = layerPositions_[layer-1];
  //get multi cluster x and y coordinates
  float matchedMultiCluster_phi = matchedMultiCluster.second.phi() - M_PI;
  if (matchedMultiCluster_phi < -M_PI)
    matchedMultiCluster_phi += 2*M_PI;
  TVector2 multiClusPos(HGCal_helpers::convertToXY(matchedMultiCluster.second.eta(), matchedMultiCluster_phi, layer_z));
  //calculate radius based on current layer's z position
  auto coneRadius = getConeRadius(frontRadius_, backRadius_, layer_z);
  for (const auto & layerCluster: selectedLayerClusters)
  {
    TVector2 layerClusPos(layerCluster->x(),layerCluster->y());
    if((multiClusPos-layerClusPos).Mod()<coneRadius)
    {
      const std::vector<std::pair<DetId, float> > &hf = layerCluster.hitsAndFractions();
      DetId highest_energy_rh;
      float highest_energy=-1;
      //find most energetic cluster
      for (const auto & rh_pair: hf)
      {
        auto hit_energy = hitmap_[rh_pair.first]->energy();
        if (hit_energy>highest_energy)
        {
          highest_energy=hit_energy;
          highest_energy_rh = rh_pair.first;
        }
      }//rh loop
      if (highest_energy>0)
      {
        auto highest_energy_gp = rhtools_.getPosition(highest_energy_rh);
        TVector2 highest_energy_pos(position->x(),position->y());
        for (const auto & rh_pair: hf)
        {
          auto gp = rhtools_.getPosition(rh_pair.first);
          TVector2 rh_pos(gp.x(),gp.y());
          if((highest_energy_pos-rh_pos).Mod()<energyRadius_)
          {
            auto hit_energy = hitmap_[rh_pair.first]->energy();
            energySum+= hit_energy*energyWeights[layer-1]*1.38;
            pTSum+= rhtools_.getPt(gp, hit_energy, vz_)*energyWeights[layer-1]*1.38;
          }
        }//rh loop
      }
    }// good layer cluster if 
  }
  return std::pair<float,float>(energySum,pTSum);
}

std::vector<std::pair<reco::BasicCluster,float> > HGCalMegaClustering::getMegaClusters(

  const std::vector<SimTrack> & simTracks_,
  const std::vector<SimVertex> & simVertices_,
	//const std::vector<reco::GenParticle> & genParticles_,
	const std::vector<reco::HGCalMultiCluster> & multiClusters_,
	const HGCRecHitCollection & recHitsEE_,
	const HGCRecHitCollection & recHitsFH_,
	const HGCRecHitCollection & recHitsBH_,
	const reco::CaloClusterCollection & layerClusters_,
  const edm::HepMCProduct & hepmc_)
{

  HepMC::GenVertex *primaryVertex = *(hepmc_)->GetEvent()->vertices_begin();
  vz_ = primaryVertex->position().z() / 10.;

  std::vector<std::pair<reco::BasicCluster,float> > tmp_clusters;

  hitmap_.clear();
  for (unsigned int i = 0; i < rechitsEE.size(); ++i)
  {
    hitmap_[rechitsEE[i].detid()] = &rechitsEE[i];
  }
  for (unsigned int i = 0; i < rechitsFH.size(); ++i)
  {
    hitmap_[rechitsFH[i].detid()] = &rechitsFH[i];
  }
  for (unsigned int i = 0; i < rechitsBH.size(); ++i)
  {
    hitmap_[rechitsBH[i].detid()] = &rechitsBH[i];
  }

	std::vector<reco::BasicCluster> megaClusters;

  mySimEvent_->fill(simTracks_, simVertices_);
  HGCal_helpers::simpleTrackPropagator toHGCalPropagator(aField_);
  toHGCalPropagator.setPropagationTargetZ(layerPositions_[0]);
	std::vector<FSimTrack *> selectedGen;
	for (unsigned int i = 0; i < mySimEvent_->nTracks()	; ++i) 
	{
		//find if the particle has reached EE
		math::XYZTLorentzVectorD vtx(0, 0, 0, 0);
		int reachedEE = 0;
		FSimTrack &myTrack(mySimEvent_->track(i));
		if (std::abs(myTrack.vertex().position().z()) >= layerPositions_[0]) continue;
		//unsigned nlayers = rhtools_.lastLayerFH();
		if (myTrack.noEndVertex())  // || myTrack.genpartIndex()>=0)
		{
			HGCal_helpers::coordinates propcoords;
			auto reachesHGCal = toHGCalPropagator.propagate(
			myTrack.momentum(), myTrack.vertex().position(), myTrack.charge(), propcoords);
			vtx = propcoords.toVector();

			if (reachesHGCal && vtx.Rho() < hgcalOuterRadius_ && vtx.Rho() > hgcalInnerRadius_)
			{
        		reachedEE = 2;
      }
      else if (reachesHGCal && vtx.Rho() > hgcalOuterRadius_)
        		reachedEE = 1;
		}
		//now skim the particles
		if (abs(myTrack.type()) == pidSelected_ && reachedEE > 0 &&
			((gun_type_ == "pt" && myTrack.momentum().pt() >= GEN_engpt_*.999) ||
			 (gun_type_ != "pt" && myTrack.momentum().energy() >= GEN_engpt_*.999)))
		{
			selectedGen.push_back(&myTrack);
		}
  }// loop over sim particles

		//return empty collection of there are no multiclusters
		if (multiClusters_.size()==0)
			return megaClusters;

		//select multiclusters
		std::vector<std::pair<FSimTrack*,reco::HGCalMultiCluster*> > bestMultiClusters;
    if (GEN_engpt_<=7.5)
      bestMultiClusters=HGCal_helpers::getHighestEnergyObjectIndex(selectedGen, multiClusters_, 0.2)
    else
      bestMultiClusters=HGCal_helpers::getHighestEnergyObjectIndex(selectedGen, multiClusters_, 0.1)

    for (const auto & pair: bestMultiClusters)
    {
      float max_rechit_energy=-1;
      float energySum = 0;
      float pTSum = 0;
      for (int layer=1; layer<maxlayer+1;layer++)
      {
        std::vector<reco::CaloCluster *> selectedLayerClusters;
        auto layer_z = layerPositions_[layer-1];
        auto coneRadius = getConeRadius(frontRadius_, backRadius_, layer_z);
        TVector2 multiClusPos(HGCal_helpers::convertToXY(pair.second->eta(), pair.second->phi(), layer_z));
        for (const auto & layerCluster: layerClusters_)
        {
          const std::vector<std::pair<DetId, float> > &hf = layerCluster.hitsAndFractions();
          auto hfsize = hf.size();
          layer_of_the_cluster=0;
          if (hfsize>0)
          {
            layer_of_the_cluster = rhtools_.getLayerWithOffset(hf[0].first);
          }
          TVector2 layerClusPos(layerCluster->x(),layerCluster->y());
          if(layer==layer_of_the_cluster && layerCluster.eta()*pair.first->eta()>=0)
          {
            selectedLayerClusters.push_back(&layerCluster);
          }
          else continue;
          if((multiClusPos-layerClusPos).Mod()<coneRadius)
          {
            DetId highest_energy_rh;
            float highest_energy=-1;
            //find most energetic cluster
            for (const auto & rh_pair: hf)
            {
              auto hit_energy = hitmap_[rh_pair.first]->energy();
              if (hit_energy>highest_energy)
              {
                highest_energy=hit_energy;
                highest_energy_rh = rh_pair.first;
              }
            }//rh loop
            if (highest_energy>0)
            {
              auto highest_energy_gp = rhtools_.getPosition(highest_energy_rh);
              TVector2 highest_energy_pos(position->x(),position->y());
              for (const auto & rh_pair: hf)
              {
                auto gp = rhtools_.getPosition(rh_pair.first);
                TVector2 rh_pos(gp.x(),gp.y());
                if((highest_energy_pos-rh_pos).Mod()<energyRadius_)
                {
                  auto hit_energy = hitmap_[rh_pair.first]->energy();
                  energySum+= hit_energy*energyWeights[layer-1]*1.38;
                  pTSum+= rhtools_.getPt(gp, hit_energy, vz_)*energyWeights[layer-1]*1.38;
                }
              }//rh loop
            }
          }// good layer cluster if       
        }//layer cluster loop
        if(doPileupSubtraction_)
        {
          auto pileup_quantities = pileupSubtraction(pair, selectedLayerClusters, layer);
          energySum -= pileup_quantities.first;
          pTSum -= pileup_quantities.second;
        }
      }//loop over layers
      reco::BasicCluster tmp_bc(energySum,pair.second->position(),reco::CaloID::DET_HGCAL_ENDCAP);
      tmp_clusters.emplace_back(tmp_bc,pTSum);
    }//loop over multiclusters -- sim particle pairs
    return tmp_clusters;
	}//getMegaClusters


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
