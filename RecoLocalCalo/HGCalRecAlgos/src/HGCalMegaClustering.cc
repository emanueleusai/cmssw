#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalMegaClustering.h"

//Access the geometry of the HGCAL detector
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

//Access data structures of calorimeter-type detectors
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"

//General helper functions
#include "DataFormats/Math/interface/deltaR.h"

//HGCAL-specific helper functions
namespace HGCal_helpers {


//defines a set of coordinated in space
class coordinates {
public:
    coordinates(): x(0), y(0), z(0), eta(100), phi(0) {}
    float x, y, z, eta, phi;
    inline math::XYZTLorentzVectorD toVector() {
        return math::XYZTLorentzVectorD(x, y, z, 0);
    }
};


//propagate the fitted track of a (charged) particle in space, given a specific magnetic field configuration
class simpleTrackPropagator {
public:
    //simple constructor requesting the magnetic field configuration
    simpleTrackPropagator(MagneticField const * f):
        field_(f), prod_(field_, alongMomentum, 5.e-5), absz_target_(0) {
        ROOT::Math::SMatrixIdentity id;
        AlgebraicSymMatrix55 C(id);
        C *= 0.001;
        err_ = CurvilinearTrajectoryError(C);
    }

    //set the Z coordinate for target plane of the propagated track
    void setPropagationTargetZ(const float & z);

    //propagate the particle to the defined Z=const plane
    //px,py,pz: momentum of the particle at position (x,y,z)
    //the position of the propagated particle is stored in "coordinates"
    //returns true if propagation succeeded
    bool propagate(const double px,
                   const double py,
                   const double pz,
                   const double x,
                   const double y,
                   const double z,
                   const float charge, coordinates & output) const;

    //overloading of the propagation function using the ROOT Lorentz vector class
    bool propagate(const math::XYZTLorentzVectorD & momentum,
                   const math::XYZTLorentzVectorD & position,
                   const float charge, coordinates & output) const;

private:
    //simpleTrackPropagator(): field_(0), prod_(field_, alongMomentum, 5.e-5), absz_target_(0) {}
    //internal function, returns the Runge-Kutta track propagator object
    const RKPropagatorInS & RKProp() const {
        return prod_.propagator;
    }
    //internal variables of the propagator
    Plane::PlanePointer targetPlaneForward_, targetPlaneBackward_;
    MagneticField const * field_;
    CurvilinearTrajectoryError err_;
    defaultRKPropagator::Product prod_;
    float absz_target_;
};


//builds a Z=const plane
void simpleTrackPropagator::setPropagationTargetZ(const float & z) {
    targetPlaneForward_ = Plane::build(Plane::PositionType(0, 0, std::abs(z)), Plane::RotationType());
    targetPlaneBackward_ =
        Plane::build(Plane::PositionType(0, 0, -std::abs(z)), Plane::RotationType());
    absz_target_ = std::abs(z);
}


//main propagator method
bool simpleTrackPropagator::propagate(const double px,
                                      const double py,
                                      const double pz,
                                      const double x,
                                      const double y,
                                      const double z,
                                      const float charge, coordinates & output) const {

    //initialize the input objects for the propagator and the propagator itself
    output = coordinates();
    typedef TrajectoryStateOnSurface TSOS;
    GlobalPoint startingPosition(x, y, z);
    GlobalVector startingMomentum(px, py, pz);
    Plane::PlanePointer startingPlane =
        Plane::build(Plane::PositionType(x, y, z), Plane::RotationType());
    TSOS startingStateP(
        GlobalTrajectoryParameters(startingPosition, startingMomentum, charge, field_), err_,
        * startingPlane);

    //propagate the status of the particles, use fwd or bkwd plane depending on the z component of the momentum of the particle
    TSOS trackStateP;
    if (pz > 0) {
        trackStateP = RKProp().propagate(startingStateP, * targetPlaneForward_);
    } else {
        trackStateP = RKProp().propagate(startingStateP, * targetPlaneBackward_);
    }
    if (trackStateP.isValid()) {
        //fill out the output coordinate object
        output.x = trackStateP.globalPosition().x();
        output.y = trackStateP.globalPosition().y();
        output.z = trackStateP.globalPosition().z();
        output.phi = trackStateP.globalPosition().phi();
        output.eta = trackStateP.globalPosition().eta();
        return true;
    }
    return false;
}


//overloading of the propagator
bool simpleTrackPropagator::propagate(const math::XYZTLorentzVectorD & momentum,
                                      const math::XYZTLorentzVectorD & position,
                                      const float charge,
                                      coordinates & output) const {

    return propagate(momentum.px(), momentum.py(), momentum.pz(), position.x(), position.y(),
                     position.z(), charge, output);
}


//matches simulated tracks (particles) with the highest energy multicluster within
//dR = sqrt( dEta**2 + dPhi**2) < deltaR
std::vector < std::pair <const FSimTrack *, const reco::HGCalMultiCluster * > >
getHighestEnergyObjectIndex(const std::vector < FSimTrack * > & reference,
                            const std::vector < reco::HGCalMultiCluster > & object,
                            const float & deltaR = 0.1) {
    //define output object
    std::vector < std::pair <
    const FSimTrack *,
          const reco::HGCalMultiCluster * > > pairings;

    //loop over simulated tracks
    for (unsigned int i = 0; i < reference.size(); i++) {
        int index = -1;
        float maxenergy = -1.0;

        //loop over multiclusters
        for (unsigned int j = 0; j < object.size(); j++) {
            if (object[j].energy() > maxenergy && //energy threshold for the multicluster
                    //proximity condition between the track and the multicluster
                    reco::deltaR < math::XYZTLorentzVectorD, reco::HGCalMultiCluster > (reference[i] - > momentum(), object[j]) < deltaR) {
                //find highest energy multicluster satisfying the conditions
                index = j;
                maxenergy = object[j].energy();
            }
        }//multicluster loop

        if (index > -1) {//found a match
            pairings.emplace_back(reference[i], & object[index]);//fill output vector
        }
    }//simtrack loop
    return pairings;
}

//convert a point given in (eta,phi,z) coordinates to x,y(,z)
TVector2 convertToXY(const float & eta,
                     const float & phi,
                     const float & z) {
    auto t = exp(-1. * eta);
    auto x = z * 2. * t * cos(phi) / (1. - t * t);
    auto y = z * 2. * t * sin(phi) / (1. - t * t);
    return TVector2(x, y);
}

} // HGCal_helpers


//constructor for the main class. Read parameters from the edm::ParameterSet (usually provided through a python configuration file)
HGCalMegaClustering::HGCalMegaClustering(const edm::ParameterSet & conf, std::string gun_type, float GEN_engpt, int pidSelected, float energyRadius, float frontRadius, float backRadius, bool doPileupSubtraction) {
    gun_type_ = gun_type;//what kind of "particle gun" has been used to generate the data
    GEN_engpt_ = GEN_engpt;//energy of the generated particle
    pidSelected_ = pidSelected;//particleID of the particle we want to reconstruct
    energyRadius_ = energyRadius;//clustering radius for energy deposits on the same layer
    frontRadius_ = frontRadius;
    backRadius_ = backRadius;//front and back radii of the truncated cone
    doPileupSubtraction_ = doPileupSubtraction;//subtract energy from pile-up interactions
    mySimEvent_ = new FSimEvent(conf);//full configuration of the simulated event
}

//support function computing the Z position of each HGCAL layer
//once HGCAL design is finalized this can be substituted by a set of constants read from HGCAL geometry record
void HGCalMegaClustering::retrieveLayerPositions(unsigned layers) {
    layerPositions_.clear();
    for (unsigned ilayer = 1; ilayer <= layers; ++ilayer) {
        const GlobalPoint pos = rhtools_.getPositionLayer(ilayer);
        layerPositions_.push_back(pos.z());
    }
}


//computes the radius size at a specific z given front radius, back radius and maximum radius
float HGCalMegaClustering::getConeRadius(float frontRadius, float backRadius, float z, float maxval) {
    float depthTerm = backRadius * (fabs(z) - 320.7) / (407.8 - 320.7);
    float val = frontRadius + depthTerm;
    //cone radius is capped at maxval, effectively becomes a cylinder beyond that.
    if (val > maxval) {
        return maxval;
    } else {
        return val;
    }

}


//reads event setup from the edm::EventSetup CMSSW class and initializes internal variables
void HGCalMegaClustering::getEventSetup(const edm::EventSetup & es) {

    //initialize the fastsim data structure
    edm::ESHandle < HepPDT::ParticleDataTable > pdt;
    es.getData(pdt);
    mySimEvent_ - > initializePdt( & ( * pdt));
    //initialize the RHtools support functions
    rhtools_.getEventSetup(es);

    //extracts information about the magnetic field configuration
    edm::ESHandle < MagneticField > magfield;
    es.get < IdealMagneticFieldRecord > ().get(magfield);
    aField_ = & ( * magfield);
}


//support function to compute the amount of pileup energy to be subtracted from the megacluster
//inputs: track-multicluster pair, collection of layer clusters
//output: energy-momentum pair, from pileup interactions
std::pair < float, float > HGCalMegaClustering::pileupSubtraction(
    const std::pair <
    const FSimTrack *,
    const reco::HGCalMultiCluster * > & matchedMultiCluster,
    const std::vector <
    const reco::CaloCluster * > & selectedLayerClusters,
    int layer) {

    float pTSum = 0;
    float energySum = 0;
    auto layer_z = layerPositions_[layer - 1];

    //get multicluster X,Y coordinates
    float matchedMultiCluster_phi = matchedMultiCluster.second - > phi() - M_PI;
    if (matchedMultiCluster_phi < -M_PI)
        matchedMultiCluster_phi += 2 * M_PI;
    TVector2 multiClusPos(HGCal_helpers::convertToXY(matchedMultiCluster.second - > eta(), matchedMultiCluster_phi, layer_z));
    //calculate radius based on current layer's z position
    auto coneRadius = getConeRadius(frontRadius_, backRadius_, layer_z);

    //loop over layer clusters
    for (const auto & layerCluster: selectedLayerClusters) {

        TVector2 layerClusPos(layerCluster - > x(), layerCluster - > y());
        //check distance between multicluster center and the layer cluster
        if ((multiClusPos - layerClusPos).Mod() < coneRadius) {

            //access the collection of detector hits (rechits) associated to a layer cluster
            const std::vector < std::pair < DetId, float > > & hf = layerCluster - > hitsAndFractions();
            DetId highest_energy_rh;
            float highest_energy = -1;
            //find most energetic cluster
            //loop over rechit-energy fraction pairs
            for (const auto & rh_pair: hf) {

                //retrieve the hit energy and find the highest energy rechit
                auto hit_energy = hitmap_[rh_pair.first] - > energy();
                if (hit_energy > highest_energy) {
                    highest_energy = hit_energy;
                    highest_energy_rh = rh_pair.first;
                }
            } //rechit loop

            if (highest_energy > 0) {//found a matching rechit

                //get coordinates of the highest energy rechit and convert to X,Y
                auto highest_energy_gp = rhtools_.getPosition(highest_energy_rh);
                TVector2 highest_energy_pos(highest_energy_gp.x(), highest_energy_gp.y());
                //loop over rechits andsum rechit energies
                //within energyRadius from the most energetic rechit
                for (const auto & rh_pair: hf) {

                    auto gp = rhtools_.getPosition(rh_pair.first);
                    TVector2 rh_pos(gp.x(), gp.y());
                    //check distance condition
                    if ((highest_energy_pos - rh_pos).Mod() < energyRadius_) {

                        auto hit_energy = hitmap_[rh_pair.first] - > energy();
                        //sum energies and momenta weighed by layer
                        //the weights correct for EM--hadronic energy fractions
                        energySum += hit_energy * energyWeights[layer - 1] * 1.38;
                        pTSum += rhtools_.getPt(gp, hit_energy, vz_) * energyWeights[layer - 1] * 1.38;
                    }
                } //rh loop
            }
        } // good layer cluster condition
    }
    return std::pair < float, float > (energySum, pTSum);
}


//main method producing the megaclusters in a collision event
//inputs: simulated tracks and vertices, multicluster collection,
//rechits, layer clusters, MC generator information
//outputs: vector of pairs megacluster-transverse momentum
std::vector < std::pair < reco::BasicCluster, float > > HGCalMegaClustering::getMegaClusters(
    const std::vector < SimTrack > & simTracks_,
    const std::vector < SimVertex > & simVertices_,
    const std::vector < reco::HGCalMultiCluster > & multiClusters_,
    const HGCRecHitCollection & recHitsEE_,
    const HGCRecHitCollection & recHitsFH_,
    const HGCRecHitCollection & recHitsBH_,
    const reco::CaloClusterCollection & layerClusters_,
    const edm::HepMCProduct & hepmc_) {

    //access primary vertices, needed by RHTools
    auto primaryVertex = hepmc_.GetEvent() - > vertices_begin();
    vz_ = ( * primaryVertex) - > position().z() / 10.;
    //output collection
    std::vector < std::pair < reco::BasicCluster, float > > tmp_clusters;

    //put together all the reconstructed detector hits (rechits)
    //from the different HGCAL sections
    //(EE: electormagnetic section; FH,BH: hadronic sections)
    hitmap_.clear();
    for (unsigned int i = 0; i < recHitsEE_.size(); ++i) {
        hitmap_[recHitsEE_[i].detid()] = & recHitsEE_[i];
    }
    for (unsigned int i = 0; i < recHitsFH_.size(); ++i) {
        hitmap_[recHitsFH_[i].detid()] = & recHitsFH_[i];
    }
    for (unsigned int i = 0; i < recHitsBH_.size(); ++i) {
        hitmap_[recHitsBH_[i].detid()] = & recHitsBH_[i];
    }

    //access the simulated event structure and prepare the track propagator
    mySimEvent_ - > fill(simTracks_, simVertices_);
    HGCal_helpers::simpleTrackPropagator toHGCalPropagator(aField_);
    toHGCalPropagator.setPropagationTargetZ(layerPositions_[0]);
    std::vector < FSimTrack * > selectedGen;
    //loop over all simulated tracks
    for (unsigned int i = 0; i < mySimEvent_ - > nTracks(); ++i) {

        //find if the particle has reached EE
        math::XYZTLorentzVectorD vtx(0, 0, 0, 0);
        int reachedEE = 0;
        FSimTrack & myTrack(mySimEvent_ - > track(i));
        //if the skip track's vertex is after the first layer of EE, skip particle
        //it would not be reconstructed by the tracker
        if (std::abs(myTrack.vertex().position().z()) >= layerPositions_[0]) continue;
        //check particle does not naturally decay (and therefore cannot be propagated)
        if (myTrack.noEndVertex())
        {

            //make sure the particle can be propagated to HGCAL entry surface
            HGCal_helpers::coordinates propcoords;
            auto reachesHGCal = toHGCalPropagator.propagate(
                                    myTrack.momentum(), myTrack.vertex().position(), myTrack.charge(), propcoords);
            vtx = propcoords.toVector();

            if (reachesHGCal && vtx.Rho() < hgcalOuterRadius_ && vtx.Rho() > hgcalInnerRadius_) {
                //the particle reaches HGCAL EE surface and is inside the HGCAL disk
                reachedEE = 2;
            } else if (reachesHGCal && vtx.Rho() > hgcalOuterRadius_)
                //the particle reaches EE but is outside of the disk
                reachedEE = 1;
        }

        //now select charged pions with energy smaller than the generated energy
        //this part of code assumes the algorithm is fed a type of
        //simulation called "particle gun"
        if (abs(myTrack.type()) == pidSelected_ && reachedEE > 0 &&
                ((gun_type_ == "pt" && myTrack.momentum().pt() >= GEN_engpt_ ) ||
                 (gun_type_ != "pt" && myTrack.momentum().energy() >= GEN_engpt_ ))) {
            selectedGen.push_back( & myTrack);
        }
    } // loop over sim particles

    //return empty collection if there are no multiclusters
    if (multiClusters_.size() == 0)
        return tmp_clusters;

    //match multiclusters with the selected sim tracks and store
    //the pairs in bestMultiClusters
    std::vector < std::pair <
    const FSimTrack *,
          const reco::HGCalMultiCluster * > > bestMultiClusters;
    //use momentum-dependent deltaR selection criteria
    if (GEN_engpt_ <= 7.5)
        bestMultiClusters = HGCal_helpers::getHighestEnergyObjectIndex(selectedGen, multiClusters_, 0.2);
    else
        bestMultiClusters = HGCal_helpers::getHighestEnergyObjectIndex(selectedGen, multiClusters_, 0.1);

    //loop over track-multicluster pairs
    for (const auto & pair: bestMultiClusters) {

        float energySum = 0;
        float pTSum = 0;

        //loop over HGCAL layers
        for (unsigned int layer = 1; layer < maxlayer + 1; layer++) {

            std::vector <
            const reco::CaloCluster * > selectedLayerClusters;

            //compute cone-related variables
            auto layer_z = layerPositions_[layer - 1];
            auto coneRadius = getConeRadius(frontRadius_, backRadius_, layer_z);
            TVector2 multiClusPos(HGCal_helpers::convertToXY(pair.second - > eta(), pair.second - > phi(), layer_z));

            //loop over all layer clusters on a single HGCAL layer
            for (const auto & layerCluster: layerClusters_) {

                //access layer cluster position and associated hits
                const std::vector < std::pair < DetId, float > > & hf = layerCluster.hitsAndFractions();
                auto hfsize = hf.size();
                unsigned int layer_of_the_cluster = 0;
                if (hfsize > 0) {
                    layer_of_the_cluster = rhtools_.getLayerWithOffset(hf[0].first);
                }
                TVector2 layerClusPos(layerCluster.x(), layerCluster.y());

                //make sure the cluster and the particle are on the same side of the detector
                if (layer == layer_of_the_cluster && layerCluster.eta() * pair.first - > momentum().eta() >= 0) {
                    selectedLayerClusters.push_back( & layerCluster);
                } else continue;

                //consider only layer clusters within the cone
                if ((multiClusPos - layerClusPos).Mod() < coneRadius) {

                    DetId highest_energy_rh;
                    float highest_energy = -1;
                    //find the most energetic rechit
                    //loop over rechits
                    for (const auto & rh_pair: hf) {

                        auto hit_energy = hitmap_[rh_pair.first] - > energy();
                        if (hit_energy > highest_energy) {
                            highest_energy = hit_energy;
                            highest_energy_rh = rh_pair.first;
                        }
                    } //rechit loop

                    //found suitable rechit
                    if (highest_energy > 0) {
                        //extract position of the highest energy rechit
                        auto highest_energy_gp = rhtools_.getPosition(highest_energy_rh);
                        TVector2 highest_energy_pos(highest_energy_gp.x(), highest_energy_gp.y());
                        //loop over rechits
                        for (const auto & rh_pair: hf) {

                            //sum the energy of the rechits closest to the highest energy RH
                            auto gp = rhtools_.getPosition(rh_pair.first);
                            TVector2 rh_pos(gp.x(), gp.y());
                            if ((highest_energy_pos - rh_pos).Mod() < energyRadius_) {//adjacency condition

                                //sum energies and momenta with layer-dependent weights
                                //and EM--hadronic energy fraction corrections
                                auto hit_energy = hitmap_[rh_pair.first] - > energy();
                                energySum += hit_energy * energyWeights[layer - 1] * 1.38;
                                pTSum += rhtools_.getPt(gp, hit_energy, vz_) * energyWeights[layer - 1] * 1.38;
                            }
                        } //rechit loop
                    }
                } // good layer cluster condition
            } //layer cluster loop

            if (doPileupSubtraction_) {//perform pileup subtraction is configured to do so

                auto pileup_quantities = pileupSubtraction(pair, selectedLayerClusters, layer);
                energySum -= pileup_quantities.first;
                pTSum -= pileup_quantities.second;
            }
        } //loop over layers

        // fill CMSSW cluster data format with the mega cluster info and add it to the output vector
        reco::BasicCluster tmp_bc(energySum, pair.second - > position(), reco::CaloID::DET_HGCAL_ENDCAP);
        tmp_clusters.emplace_back(tmp_bc, pTSum);
    } //loop over multiclusters -- sim particle pairs
    return tmp_clusters;
} //getMegaClusters

