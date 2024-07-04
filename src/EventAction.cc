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
#include "G4SDManager.hh"
#include "Hit.hh"
#include "GammaRayHelper.hh"

///namespace G4FastSim
///{
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//using namespace G4FastSim;

namespace G4FastSim {

//std::mutex EventAction::mtx;

EventAction::EventAction()
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  fHitsCollectionNames.push_back("HitsCollection1");
  fHitsCollectionNames.push_back("HitsCollection2");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{
  G4cout << "EventAction::BeginOfEventAction..... NEXT" << G4endl;	
  fEdep = 1.2345;
  fEd.clear();
  fX.clear();
  fY.clear();
  fZ.clear();

  //G4cout<<"EventAction::BeginOfEventAction next event...."<<G4endl;
  G4PrimaryVertex* primaryVertex = event->GetPrimaryVertex();
  fXp = primaryVertex->GetPosition().x();
  fYp = primaryVertex->GetPosition().y();
  fZp = primaryVertex->GetPosition().z();

  //auto def =  event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition();
  //G4cout << " def = " << def->GetParticleName() << G4endl;
  //G4cout << " p = "<<event->GetPrimaryVertex()->GetPrimary()->GetMomentumDirection() << G4endl;
  // Kill the event if the particle does not point to the fiducial volume
  //if (TBD) {
  //  G4cout<<"EventAction::BeginOfEventAction: Killing event with fZp = "<<fZp<<G4endl;
  //  G4RunManager::GetRunManager()->AbortEvent();
  //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  G4cout << "EventAction::EndOfEventAction..... Analyze hits and cluster...." << G4endl;
  AnalyzeHits(event);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  //std::lock_guard<std::mutex> lock(mtx);
  analysisManager->FillNtupleDColumn(0, 0, fEdep);
  analysisManager->FillNtupleDColumn(0, 1, fXp);
  analysisManager->FillNtupleDColumn(0, 2, fYp);
  analysisManager->FillNtupleDColumn(0, 3, fZp);
  analysisManager->AddNtupleRow(0);

}

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
      G4cout << "Hits Collection: " << fHitsCollectionNames[i] << " has " << n_hit << " hits." << G4endl;
      for (G4int j = 0; j < n_hit; ++j) {
          Hit* hit = (*fHitsCollection)[j];
          hit->Print();
          allHits.push_back(hit);
      }
  }

  // cluster hits based on spatial and time thresholds
  std::vector<Cluster> fClusters;
  G4double spatialThreshold = 50.0 * mm;
  G4double timeThreshold = 10.0 * ns;
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

    for (auto& hit : hits) {
        bool addedToCluster = false;
        for (auto& cluster : clusters) {
            if (CalculateDistance(hit->position, cluster.position) < spatialThreshold &&
                CalculateTimeDifference(hit->time, cluster.time) < timeThreshold) {
                  
                G4int clusterSize = cluster.hits.size();
                cluster.position = (cluster.position * clusterSize + hit->position) / (clusterSize + 1);
                cluster.energyDeposit += hit->energyDeposit; // Sum energy deposits
                cluster.time = (cluster.time * clusterSize + hit->time) / (clusterSize + 1); // Update average time
                cluster.hits.push_back(hit);
                addedToCluster = true;
                break;
            }
        }
        if (!addedToCluster) {
            clusters.push_back(Cluster{hit->position, hit->energyDeposit, hit->time, {hit}});
        }
    }

    G4cout << "Number of clusters: " << clusters.size() << G4endl;
    // Calculate cluster positions and store into ntuple variables
    for (auto& cluster : clusters) {
        //cluster.position /= cluster.hits.size();
        fEd.push_back(cluster.energyDeposit);
        fX.push_back(cluster.position.x());
        fY.push_back(cluster.position.y());
        fZ.push_back(cluster.position.z());
    }
}


} // namespace G4FastSim