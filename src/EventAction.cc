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
#include "SensitiveDetector.hh"
#include "G4LogicalVolumeStore.hh"

///namespace G4Sim
///{
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//using namespace G4Sim;

namespace G4Sim {

//std::mutex EventAction::mtx;

EventAction::EventAction() : G4UserEventAction(), fGammaRayHelper(&GammaRayHelper::Instance())
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);
}

/**
 * @brief Adds a hits collection name to the EventAction.
 * 
 * This function adds the specified hits collection name to the EventAction's list of hits collection names.
 * 
 * @param name The name of the hits collection to be added.
 */
void EventAction::AddHitsCollectionName(const G4String& name) {
    G4cout << "EventAction::AddHitsCollectionName: " << name << G4endl;
    fHitsCollectionNames.push_back(name);
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
    // Only when here for first time we do the initialization of teh CDFs (first time is taken care of inside function).
    fGammaRayHelper->InitializeCDFs(primaryVertex->GetPrimary()->GetKineticEnergy());  
  }

  if(!fInitializedGraphs) {
    // Get the RunAction instance
    const RunAction* runAction = static_cast<const RunAction*>(G4RunManager::GetRunManager()->GetUserRunAction());
    G4double e0 = primaryVertex->GetPrimary()->GetKineticEnergy();
    runAction->DefineDifferentialCrossSectionNtuple(e0);
    fInitializedGraphs = true;
  }
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

  fLogWeight = 0.0;
  fEventType = 0;
  fXp = 0.0;
  fYp = 0.0;
  fZp = 0.0;

  // cluster information
  fE.clear();
  fX.clear();
  fY.clear();
  fZ.clear();
  fW.clear();
  fID.clear();
  // detector information
  fEdet.clear();
  fNdet.clear();
  fNphot.clear();
  fNcomp.clear();
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
  // get the energy depositis in keV
  analysisManager->FillNtupleDColumn(0, 0, fEventID);
  analysisManager->FillNtupleDColumn(0, 1, fLogWeight);
  analysisManager->FillNtupleDColumn(0, 2, fEventType);
  analysisManager->FillNtupleDColumn(0, 3, fXp);
  analysisManager->FillNtupleDColumn(0, 4, fYp);
  analysisManager->FillNtupleDColumn(0, 5, fZp);
  analysisManager->AddNtupleRow(0);

  if(verbosityLevel>0) G4cout << "EventAction::EndOfEventAction: Done...." << G4endl;	

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Renormalizes the hit times in the hits collections.
 *
 * This function iterates over all hits collections specified in `fHitsCollectionNames`,
 * retrieves the hits, and renormalizes their times such that the smallest hit time
 * in the collection is set to zero. The renormalized time for each hit is calculated
 * by subtracting the smallest hit time from each hit's time.
 *
 * The function performs the following steps:
 * 1. Iterates over all hits collections specified in `fHitsCollectionNames`.
 * 2. Retrieves the hits from each collection and stores them in a vector.
 * 3. Finds the smallest hit time in the vector of all hits.
 * 4. Renormalizes the time of each hit by subtracting the smallest hit time.
 *
 * @note The function assumes that there is at least one hit in the collections.
 */
void EventAction::RenormalizeHitTimes(G4HCofThisEvent* HCE) {
    // renormalize the Hit times...
    std::vector<Hit*> allHits;
    G4double minTime = 0.0;

    for (size_t i = 0; i < fHitsCollectionNames.size(); ++i) {
        G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollectionNames[i]);
        auto* fHitsCollection = static_cast<HitsCollection*>(HCE->GetHC(hcID));
        if (!fHitsCollection) continue;

        G4int n_hit = fHitsCollection->entries();
        if (verbosityLevel > 0) G4cout << "Hits Collection: " << fHitsCollectionNames[i] << " has " << n_hit << " hits." << G4endl;

        for (G4int j = 0; j < n_hit; ++j) {
            Hit* hit = (*fHitsCollection)[j];
            allHits.push_back(hit);
            if (j==0) minTime = hit->time;
            else if (hit->time < minTime) minTime = hit->time;
        }
    }
    // renomalize the hit times to the smallest time in the collection
    for (auto& hit : allHits) {
        hit->time -= minTime;
    }

}


/**
 * @brief Counts the number of Compton and photoelectric interactions in a list of hits.
 *
 * This function iterates through a vector of Hit pointers and increments the counters
 * for Compton interactions (`ncomp`) and photoelectric interactions (`nphot`) based on
 * the `processType` of each hit. Only interactions of the primary gamma ray are counted.
 *
 * @param hits A vector of pointers to Hit objects to be analyzed.
 * @param ncomp An integer reference to store the count of Compton interactions.
 * @param nphot An integer reference to store the count of photoelectric interactions.
 */
void EventAction::CountInteractions(std::vector<Hit*>& hits, int& ncomp, int& nphot) {


    for (const auto& hit : hits) {
        if (hit->trackID != 1) continue;  // Only primary track hits
        if (hit->processType == "compt") ncomp++;
        if (hit->processType == "phot") nphot++;
    }
}

/**
 * @brief Prints the names of all sensitive detectors and their associated hit collections.
 *
 * This function retrieves the sensitive detector manager and hit collection table,
 * then iterates through all hit collections to print their names along with the
 * names of the associated sensitive detectors.
 *
 * The output is printed to the standard Geant4 output stream (G4cout).
 */
void PrintSensitiveDetectorNames() {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    G4HCtable* hcTable = sdManager->GetHCtable();
    
    G4int numberOfCollections = hcTable->entries();
    G4cout << "Number of hit collections: " << numberOfCollections << G4endl;

    for (G4int i = 0; i < numberOfCollections; ++i) {
        G4String collectionName = hcTable->GetHCname(i);
        G4String sdName = hcTable->GetSDname(i);

        G4cout << "Hit Collection: " << collectionName << ", Sensitive Detector: " << sdName << G4endl;
    }
}


/**
 * @brief Analyzes the hits in the given event.
 *
 * This function processes the hits collected during the event, renormalizes hit times,
 * and clusters the hits based on spatial and temporal thresholds. It also counts the
 * number of Compton and photoelectric interactions.
 *
 * @param event Pointer to the G4Event object containing the event data.
 *
 * The function performs the following steps:
 * 1. Retrieves the hits collection of the event.
 * 2. Renormalizes hit times to avoid numerical issues.
 * 3. Prints the names of sensitive detectors if verbosity level is greater than 0.
 * 4. Iterates over all sensitive detectors and retrieves their hit collections.
 * 5. For each hit collection, it retrieves the hits and adds them to a list.
 * 6. Counts the number of Compton and photoelectric interactions.
 * 7. Clusters the hits based on spatial and temporal thresholds.
 * 8. Stores the clustered data along with interaction counts.
 *
 * If any sensitive detector is not found or cannot be cast to the custom SensitiveDetector class,
 * the function prints an error message and exits the program.
 */
void EventAction::AnalyzeHits(const G4Event* event) {
    G4HCofThisEvent* HCE = event->GetHCofThisEvent();
    if (!HCE) {
        G4ExceptionDescription msg;
        msg << "No hits collection of this event found." << G4endl;
        G4Exception("EventAction::EndOfEventAction()", "MyCode0001", JustWarning, msg);
        return;
    }

    // Renormalize hit times to avoid numerical issues
    RenormalizeHitTimes(HCE);

    if(verbosityLevel>0) PrintSensitiveDetectorNames();

    // Get the list of sensitive detector names
    std::vector<G4String> sensitiveDetectorNames = GetSensitiveDetectorNames();
    // Loop over all sensitive detectors
    int collectionId = 0;
    for (const G4String& sdName : sensitiveDetectorNames) {
        G4VSensitiveDetector* sd = G4SDManager::GetSDMpointer()->FindSensitiveDetector(sdName);
        if (!sd) {
            G4cerr << "Sensitive detector " << sdName << " not found!" << G4endl;
            exit(-1);
        }

        // Cast to your custom sensitive detector class
        SensitiveDetector* mySD = dynamic_cast<SensitiveDetector*>(sd);
        if (!mySD) {
            G4cerr << "Error casting to SensitiveDetector for " << sdName << G4endl;
            exit(-1);
        }

        // Now iterate over all hit collections in this sensitive detector
        std::vector<Hit*> hitList;
        G4double spatialThreshold = 100.0 * mm;
        G4double timeThreshold = 10.0 * ns;
        for (auto* hitsCollection : mySD->GetHitsCollections()) {
            G4int n_hit = hitsCollection->entries();
            if(verbosityLevel>0)
                G4cout << "Hits Collection in " << sdName << " has " << n_hit << " hits." << G4endl;
            
            G4String collectionName = hitsCollection->GetName();
            spatialThreshold = GetSpatialThreshold(collectionName) * mm;
            timeThreshold = GetTimeThreshold(collectionName) * ns;

            // Find the hits in this collection
            for (G4int j = 0; j < n_hit; ++j) {
                Hit* hit = (*hitsCollection)[j];
                hitList.push_back(hit);  // Add hit to the list
            }


        }    
        // Call the function that processes the hits: all hist in teh same senstive detector are processed together
        // Count Compton and photoelectric interactions
        G4int ncomp = 0;
        G4int nphot = 0;
        CountInteractions(hitList, ncomp, nphot);
        // Cluster hits and store the data
        std::vector<Cluster> clusters;
        ClusterHits(hitList, spatialThreshold, timeThreshold, clusters, static_cast<int>(collectionId));  
        // Store the data for each collection
        StorePerCollectionData(clusters, ncomp, nphot);
        // Increment the collection ID
        collectionId++;
    }

    return;
}

/**
 * @brief Retrieves the names of all sensitive detectors in the simulation.
 *
 * This method accesses the G4SDManager to get the list of all hits collections
 * and extracts the names of the sensitive detectors associated with them.
 * It ensures that the returned list contains unique names only.
 *
 * @return A vector of G4String containing the names of all sensitive detectors.
 */
std::vector<G4String> EventAction::GetSensitiveDetectorNames() {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    std::vector<G4String> sdNames;

    // Retrieve the list of all hits collections and extract the sensitive detector names from them
    const G4HCtable* hctable = sdManager->GetHCtable();
    for (G4int i = 0; i < hctable->entries(); ++i) {
        G4String hcName = hctable->GetHCname(i);
        G4String sdName = hctable->GetSDname(i);  // Get the corresponding sensitive detector name
        if (std::find(sdNames.begin(), sdNames.end(), sdName) == sdNames.end()) {
            sdNames.push_back(sdName);  // Avoid duplicates
        }
    }

    return sdNames;
}

/**
 * Stores the cluster data for each tag group.
 *
 * @param clusters The clusters for the given tag group.
 * @param ncomp The number of Compton interactions.
 * @param nphot The number of photoelectric interactions.
 */
void EventAction::StorePerCollectionData(const std::vector<Cluster>& clusters, G4int ncomp, G4int nphot) {
    G4double edet = 0.0;
    G4int nclus = 0;

    for (const auto& cluster : clusters) {
        if (cluster.energyDeposit > 0 * eV) {
            nclus++;
            edet += cluster.energyDeposit / keV;

            fE.push_back(cluster.energyDeposit / keV);
            fX.push_back(cluster.position.x());
            fY.push_back(cluster.position.y());
            fZ.push_back(cluster.position.z());
            fID.push_back(cluster.collectionID);
            fW.push_back(fLogWeight);
        }
    }

    // Store the total energy deposit and number of clusters per collection (tag)
    fEdet.push_back(edet);
    fNdet.push_back(nclus);
    fNphot.push_back(nphot);
    fNcomp.push_back(ncomp);
}

/**
 * Clusters hits based on spatial and time thresholds.
 *
 * @param hits The vector of hits to be clustered.
 * @param spatialThreshold The maximum spatial distance for hits to be considered part of the same cluster.
 * @param timeThreshold The maximum time difference for hits to be considered part of the same cluster.
 * @param clusters The vector of clusters to store the clustered hits.
 * @param collectionID The collection ID for the hits.
 */
void EventAction::ClusterHits(std::vector<Hit*>& hits, G4double spatialThreshold, G4double timeThreshold, std::vector<Cluster>& clusters, int collectionID) {
    if (hits.empty()) return;

    // Clustering the hits
    for (auto& hit : hits) {
        if (verbosityLevel > 0) hit->Print();
        if (hit->used) continue;

        bool addedToCluster = false;
        for (auto& cluster : clusters) {
            if (CalculateDistance(hit->position, cluster.position) < spatialThreshold &&
                CalculateTimeDifference(hit->time, cluster.time) < timeThreshold) {

                G4double energyDeposit = hit->energyDeposit;
                if (energyDeposit > 0 * eV) {
                    G4int clusterSize = cluster.hits.size();
                    cluster.position = (cluster.position * clusterSize + hit->position) / (clusterSize + 1);
                    cluster.energyDeposit += energyDeposit;
                    cluster.time = (cluster.time * clusterSize + hit->time) / (clusterSize + 1);
                    cluster.hits.push_back(hit);
                    addedToCluster = true;
                    break;
                }
            }
        }

        // If the hit was not added to any existing cluster, create a new cluster
        if (!addedToCluster) {
            clusters.push_back(Cluster{hit->position, hit->energyDeposit, hit->time, {hit}, collectionID});
        }
    }

    // Merge close clusters
    MergeClusters(clusters, spatialThreshold, timeThreshold);
}

/**
 * Merges clusters that are close in space and time.
 *
 * @param clusters The vector of clusters to be merged.
 * @param spatialThreshold The maximum spatial distance for merging clusters.
 * @param timeThreshold The maximum time difference for merging clusters.
 */
void EventAction::MergeClusters(std::vector<Cluster>& clusters, G4double spatialThreshold, G4double timeThreshold) {
    for (size_t i = 0; i < clusters.size(); ++i) {
        for (size_t j = i + 1; j < clusters.size();) {
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
                ++j;  // Increment only if no merge
            }
        }
    }
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

G4double EventAction::GetSpatialThreshold(const G4String& collectionName) {
    if (fClusteringParameters.find(collectionName) != fClusteringParameters.end()) {
        return fClusteringParameters[collectionName].first;
    }
    return 10.0 * mm;  // default value if not found
}

G4double EventAction::GetTimeThreshold(const G4String& collectionName) {
    if (fClusteringParameters.find(collectionName) != fClusteringParameters.end()) {
        return fClusteringParameters[collectionName].second;
    }
    return 100.0 * ns;  // default value if not found
}

std::map<G4String, std::pair<G4double, G4double>> EventAction::fClusteringParameters;

void EventAction::SetClusteringParameters(const std::map<G4String, std::pair<G4double, G4double>>& params) {
    fClusteringParameters = params;
}

} // namespace G4Sim