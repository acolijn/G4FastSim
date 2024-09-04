//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1/src/EventAction.cc
/// \brief Implementation of the G4FastSim::EventAction class

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4EventManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Navigator.hh"
#include "G4SDManager.hh"
#include "Hit.hh"
#include "GammaRayHelper.hh"

///namespace G4FastSim
///{
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//using namespace G4FastSim;

namespace G4FastSim {

//std::mutex EventAction::mtx;

EventAction::EventAction() : G4UserEventAction(), fGammaRayHelper(&GammaRayHelper::Instance())
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  //fHitsCollectionNames.push_back("LXeCollection");
  fHitsCollectionNames.push_back("LXeFiducialCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{
  if (verbosityLevel > 0)
    G4cout << "EventAction::BeginOfEventAction..... NEXT" << G4endl;	

  // Reset variables
  ResetVariables();

  //G4cout<<"EventAction::BeginOfEventAction next event...."<<G4endl;
  G4PrimaryVertex* primaryVertex = event->GetPrimaryVertex();
  fXp = primaryVertex->GetPosition().x();
  fYp = primaryVertex->GetPosition().y();
  fZp = primaryVertex->GetPosition().z();

  //auto def =  event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition();
  //G4cout << " def = " << def->GetParticleName() << G4endl;
  if (verbosityLevel > 0){
    G4cout << "EventAction::BeginOfEventAction Primary vertex: x = "<< primaryVertex->GetPosition() / cm << " (cm)"<<G4endl;
    G4cout << "                                                p = "<< primaryVertex->GetPrimary()->GetMomentumDirection() << G4endl;
    G4cout << "                                                E = "<< primaryVertex->GetPrimary()->GetKineticEnergy() / keV << " keV"<< G4endl;
  }

  if(IsFastSimulation()) {
    fGammaRayHelper->InitializeCDFs(primaryVertex->GetPrimary()->GetKineticEnergy());  // only when here for first time we do the initialization
  
    // Kill the event if the particle does not point to the fiducial volume
    //if (!IsWithinFiducialVolume(primaryVertex->GetPosition(), primaryVertex->GetPrimary()->GetMomentumDirection())){
    //  G4cout<<"EventAction::BeginOfEventAction: Killing event..."<<G4endl;
    //  G4RunManager::GetRunManager()->AbortEvent();
    //}
  }
  if(!fInitializedGraphs) {
    // Get the RunAction instance
    const RunAction* runAction = static_cast<const RunAction*>(G4RunManager::GetRunManager()->GetUserRunAction());
    G4double e0 = primaryVertex->GetPrimary()->GetKineticEnergy();
    runAction->DefineDifferentialCrossSectionNtuple(e0);
    fInitializedGraphs = true;
  }
}



bool EventAction::IsWithinFiducialVolume(const G4ThreeVector& position, const G4ThreeVector& direction) {
    G4Navigator* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

    // Define a point far along the direction to ensure it goes through the volume if it intersects
    G4ThreeVector farPoint = position + 1000 * cm * direction;

    G4ThreeVector intersectionPoint;
    G4double stepLength;

    // Check for intersection with the fiducial volume
    if (navigator->LocateGlobalPointAndSetup(position, 0, false)->GetLogicalVolume() == fFiducialVolume) {
        return true;
    }

    // Check if the far point intersects the fiducial volume
    if (navigator->ComputeStep(position, direction, 1000 * cm, stepLength) > 0) {
        intersectionPoint = position + stepLength * direction;
        if (navigator->LocateGlobalPointAndSetup(intersectionPoint, 0, false)->GetLogicalVolume() == fFiducialVolume) {
            return true;
        }
    }

    return false;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::ResetVariables() {
  fNumberOfScatters = 0;
  // this is used in the standard MC to see if the gamma ray has already been in xenon or not
  fHasBeenInXenon = false;
  // set the type of event: "direct_gamma (=0)",  "scattered_gamma (=1)",
  // this is set to "direct_gamma" in the beginning and changed to "scattered_gamma" if a scatter is made	prior to entering the xenon target	
  fEventType = DIRECT_GAMMA;
  // the avalaible energy is the maximum energy that can be deposited in the event
  // it will be reduced after every energy deposit

  fAvailableEnergy = fMaxEnergy;
  //G4cout << "EventAction::ResetVariables: fAvailableEnergy = " << fAvailableEnergy << G4endl;
  //G4cout << "EventAction::ResetVariables: fMaxEnergy = " << fMaxEnergy << G4endl;

  fEdep = 0.0;
  fLogWeight = 0.0;
  fNclusters = 0;
  fNphot = 0;
  fNcomp = 0;
  fXp = 0.0;
  fYp = 0.0;
  fZp = 0.0;
  fEd.clear();
  fX.clear();
  fY.clear();
  fZ.clear();
  fW.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  if(verbosityLevel>0) G4cout << "EventAction::EndOfEventAction..... Analyze hits and cluster...." << G4endl;
  AnalyzeHits(event);

  // if no scatters were made we are dealing with an event that did nothing inside the fiducial volume
  // such an event should have a weight=1.0
  if(fNumberOfScatters == 0) fLogWeight = 0.0;
  const G4Event* currentEvent = G4EventManager::GetEventManager()->GetConstCurrentEvent();
  fEventID = currentEvent->GetEventID();
  
  if(verbosityLevel<0) G4cout << "EventAction::EndOfEventAction..... Fill ntuple...." << G4endl;
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  //std::lock_guard<std::mutex> lock(mtx);
  
  // get the energy depositis in keV
  analysisManager->FillNtupleDColumn(0, 0, fEventID);
  analysisManager->FillNtupleDColumn(0, 1, fNclusters);
  analysisManager->FillNtupleDColumn(0, 2, fNphot);
  analysisManager->FillNtupleDColumn(0, 3, fNcomp);
  analysisManager->FillNtupleDColumn(0, 4, fEdep);
  analysisManager->FillNtupleDColumn(0, 5, fLogWeight);
  analysisManager->FillNtupleDColumn(0, 6, fEventType);
  analysisManager->FillNtupleDColumn(0, 7, fXp);
  analysisManager->FillNtupleDColumn(0, 8, fYp);
  analysisManager->FillNtupleDColumn(0, 9, fZp);
  analysisManager->AddNtupleRow(0);

  if(verbosityLevel>0) G4cout << "EventAction::EndOfEventAction: Done...." << G4endl;	

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * Analyzes the hits in the event.
 *
 * @param event The G4Event object representing the current event.
 */
void EventAction::AnalyzeHits(const G4Event* event) {
  G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  if (!HCE) {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found." << G4endl;
      G4Exception("EventAction::EndOfEventAction()", "MyCode0001", JustWarning, msg);
      return;
  }

  std::vector<Hit*> allHits;

  for (size_t i = 0; i < fHitsCollectionNames.size(); ++i) {
      G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollectionNames[i]);
      auto *fHitsCollection = static_cast<HitsCollection*>(HCE->GetHC(hcID));
      
      if (!fHitsCollection) continue;

      G4int n_hit = fHitsCollection->entries();
      if(verbosityLevel>0) G4cout << "Hits Collection: " << fHitsCollectionNames[i] << " has " << n_hit << " hits." << G4endl;
      for (G4int j = 0; j < n_hit; ++j) {
          Hit* hit = (*fHitsCollection)[j];
          //if(verbosityLevel>0) hit->Print();
          allHits.push_back(hit);
      }
  }

  // cluster hits based on spatial and time thresholds
  std::vector<Cluster> fClusters;
  G4double spatialThreshold = 10 * cm;
  G4double timeThreshold = 2.0 * ns;
  ClusterHits(allHits, spatialThreshold, timeThreshold, fClusters); 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EventAction::CalculateDistance(const G4ThreeVector& pos1, const G4ThreeVector& pos2) {
    return (pos1 - pos2).mag();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EventAction::CalculateTimeDifference(G4double time1, G4double time2) {
    return std::fabs(time1 - time2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


/**
 * @brief Clusters hits based on spatial and time thresholds.
 * 
 * This function takes a vector of hits and clusters them based on their spatial and time proximity.
 * Hits that are within the specified spatial and time thresholds are added to the same cluster.
 * If a hit does not belong to any existing cluster, a new cluster is created for that hit.
 * 
 * @param hits The vector of hits to be clustered.
 * @param spatialThreshold The maximum spatial distance for hits to be considered part of the same cluster.
 * @param timeThreshold The maximum time difference for hits to be considered part of the same cluster.
 * @param clusters The vector of clusters to store the clustered hits.
 */
void EventAction::ClusterHits(std::vector<Hit*>& hits, G4double spatialThreshold, G4double timeThreshold, std::vector<Cluster>& clusters) {

    if (hits.empty()) return;  // No hits, nothing to do.

    // Find the earliest hit time to normalize times relative to the start of the event.
    G4double startTime = hits[0]->time;
    for (const auto& hit : hits) {
        if (hit->time < startTime) {
            startTime = hit->time;
        }
    }

    // Normalize hit times to the start of the event.
    for (auto& hit : hits) {
        hit->time -= startTime;
    }

    //G4cout<<G4endl;
    //G4cout<<"ClusterHits: hits.size() = "<<hits.size()<<G4endl;
    // 
    // Finding the cluster seeds. A cluster seed is a hit from a primary track doing a Compton or photoelectric interaction.
    // For the fast simulation the primary track can have another trackID, so we do not need to check the ID.
    //
    for (auto& hit : hits){      
      G4String process = hit->processType;
      G4int trackID = hit->trackID;

      G4bool isRelevantProcess = (process == "compt" || process == "phot");
      G4bool isPrimaryTrack = (trackID == 1);

      if (IsFastSimulation() || isPrimaryTrack) {
          if (process == "compt") fNcomp++;
          if (process == "phot") fNphot++;
          if (isRelevantProcess) {
              clusters.push_back(Cluster{hit->position, hit->energyDeposit, hit->time, {hit}});
              hit->used = true;
          }
      }
    }

    //
    // Clustering the hits, with the seeds already in the cluster vector
    //
    for (auto& hit : hits) {
        if(verbosityLevel>0 && (fNcomp+fNphot == 1)) hit->Print();
        if (hit->used) continue;

        bool addedToCluster = false;
        for (auto& cluster : clusters) {
            if (CalculateDistance(hit->position, cluster.position) < spatialThreshold &&
                CalculateTimeDifference(hit->time, cluster.time) < timeThreshold) {
                
                G4double energyDeposit = hit->energyDeposit;
                if(energyDeposit > 0*eV) {
                  G4int clusterSize = cluster.hits.size();
                  cluster.position = (cluster.position * clusterSize + hit->position) / (clusterSize + 1);
                  cluster.energyDeposit += energyDeposit; // Sum energy deposits
                  cluster.time = (cluster.time * clusterSize + hit->time) / (clusterSize + 1); // Update average time
                  cluster.hits.push_back(hit);
                  addedToCluster = true;
                  break;
                }
            }
        }
        // If the hit was not added to any existing cluster, create a new cluster
        if (!addedToCluster) {
            clusters.push_back(Cluster{hit->position, hit->energyDeposit, hit->time, {hit}});
        }
    }

    // Merge close clusters

    for (size_t i = 0; i < clusters.size(); ++i) {
        for (size_t j = i + 1; j < clusters.size(); ) {
            if (CalculateDistance(clusters[i].position, clusters[j].position) < spatialThreshold &&
                CalculateTimeDifference(clusters[i].time, clusters[j].time) < timeThreshold) {

                // Merge cluster j into cluster i
                G4int totalHits = clusters[i].hits.size() + clusters[j].hits.size();
                clusters[i].position = (clusters[i].position * clusters[i].hits.size() + clusters[j].position * clusters[j].hits.size()) / totalHits;
                clusters[i].energyDeposit += clusters[j].energyDeposit;
                clusters[i].time = (clusters[i].time * clusters[i].hits.size() + clusters[j].time * clusters[j].hits.size()) / totalHits;
                clusters[i].hits.insert(clusters[i].hits.end(), clusters[j].hits.begin(), clusters[j].hits.end());

                // Remove cluster j
                clusters.erase(clusters.begin() + j);
            } else {
                ++j; // Only increment if no merge
            }
        }
    }

    //G4cout << "Number of clusters: " << clusters.size() << G4endl;
    // Calculate cluster positions and store into ntuple variables
    fNclusters = 0;

    for (auto& cluster : clusters) {
      if(verbosityLevel>1  && (fNcomp+fNphot == 1)){
        G4cout <<"Cluster: " << G4endl;
        G4cout <<"                 position: " << cluster.position << G4endl;
        G4cout <<"                 energyDeposit: " << cluster.energyDeposit << G4endl;
        G4cout <<"                 weight: " << fLogWeight << G4endl;
      }
      //cluster.position /= cluster.hits.size();
      if (cluster.energyDeposit > 0*eV) {
        fNclusters++;
        fEdep += cluster.energyDeposit / keV;
        fEd.push_back(cluster.energyDeposit / keV);
        fX.push_back(cluster.position.x());
        fY.push_back(cluster.position.y());
        fZ.push_back(cluster.position.z());
        fW.push_back(fLogWeight);
      }
    }

    // these are the weird events that I do not understand.......
    if ((fEdep < 800.0) && (fEdep > 0.0) && (fNphot == 1 ) && (fEventType == DIRECT_GAMMA) && (fNcomp == 0)) {
      G4cout << G4endl;
      G4cout << "EventAction::ClusterHits: fEdep = " << fEdep << G4endl;
      G4cout << "EventAction::ClusterHits: fNphot = " << fNphot << G4endl;
      G4cout << "EventAction::ClusterHits: fNcomp = " << fNcomp << G4endl;
      G4cout << "EventAction::ClusterHits: fEventType = " << fEventType << G4endl;
      G4cout << "EventAction::ClusterHits: fNclusters = " << fNclusters << G4endl;
      for (auto& cluster : clusters) {
        G4cout << "Cluster: " << G4endl;
        G4cout << "                 position: " << cluster.position << G4endl;
        G4cout << "                 energyDeposit: " << cluster.energyDeposit << G4endl;
        G4cout << "                 weight: " << fLogWeight << G4endl;
      }
      G4cout << G4endl;
      for (auto& hit : hits) {
        hit->Print();
      }
    }

}


} // namespace G4FastSim