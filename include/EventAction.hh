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
/// \file B1/include/EventAction.hh
/// \brief Definition of the B1::EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"	
#include "G4String.hh"
#include "G4Types.hh"
#include "Cluster.hh"
#include "globals.hh"

#include "GammaRayHelper.hh"
#include <vector>
#include <mutex>

///
/// Event action class
///

class G4Event;
class G4HCofThisEvent;
class G4VHitsCollection;

/**
 * @namespace G4Sim
 * @brief A namespace for fast simulation classes in Geant4.
/*/
namespace G4Sim
{

// Use enum to define the constants
enum EventType {
    DIRECT_GAMMA = 0,
    SCATTERED_GAMMA = 1
};

/**
 * @class EventAction
 * @brief Class responsible for defining actions to be taken at the beginning and end of each event.
 *
 * This class inherits from G4UserEventAction and is used to define the actions to be taken at the beginning and end of each event.
 * It provides methods to analyze the hits in the event, reset the variables for each event, and access the variables stored in the ntuple tree.
 * It also provides methods to set and get the number of scatters, maximum energy, available energy, weight, and event type.
 */
class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event);// override;
    void EndOfEventAction(const G4Event* event);// override;

    std::vector<G4double>& GetX(){return fX;};
    std::vector<G4double>& GetY(){return fY;};
    std::vector<G4double>& GetZ(){return fZ;};
    std::vector<G4double>& GetE(){return fE;};
    std::vector<G4double>& GetW(){return fW;};
    std::vector<G4int>& GetID(){return fID;};

    // per detector information
    std::vector<G4double>& GetEdet(){return fEdet;};
    std::vector<G4int>& GetNdet(){return fNdet;};
    std::vector<G4int>& GetNphot(){return fNphot;};
    std::vector<G4int>& GetNcomp(){return fNcomp;};

    void AddWeight(G4double weight) { fLogWeight += weight; }
    void AnalyzeHits(const G4Event* event);
    void ResetVariables();

    G4bool IsFastSimulation() { return fFastSimulation; }
    G4int GetNumberOfScatters() { return fNumberOfScatters;}
    G4double GetMaxEnergy() { return fMaxEnergy; }
    G4int GetNumberOfScattersMax() { return fNumberOfScattersMax; }
    G4double GetAvailableEnergy() { return fAvailableEnergy; }
    G4double GetWeight() { return fLogWeight; }
    G4bool HasBeenInXenon() { return fHasBeenInXenon; }
    G4int GetEventType() { return fEventType; }

    void SetNumberOfScattersMax(G4int n) { fNumberOfScattersMax = n; }
    void SetFastSimulation(G4bool fast) { fFastSimulation = fast; }
    void SetNumberOfScatters(G4int n) { fNumberOfScatters = n; }
    void SetMaxEnergy(G4double e) { fMaxEnergy = e; }
    void SetAvailableEnergy(G4double e) { fAvailableEnergy = e; }
    void ReduceAvailableEnergy(G4double e) { fAvailableEnergy -= e; } 
    void SetHasBeenInXenon(G4bool b) { fHasBeenInXenon = b; }
    void SetEventType(G4int type) { fEventType = type; }

    G4double GetSpatialThreshold(const G4String& collectionName);
    G4double GetTimeThreshold(const G4String& collectionName);
    static void SetClusteringParameters(const std::map<G4String, std::pair<G4double, G4double>>& params);
    void AddHitsCollectionName(const G4String& name);
    void RenormalizeHitTimes(G4HCofThisEvent* HCE);
    void CountInteractions(std::vector<Hit*>& hits, G4int& ncomp, G4int& nphot);

  private:
    //
    // functions for hit clustering
    //
    G4double fSpatialThreshold;  // Spatial threshold for clustering
    G4double fTimeThreshold;     // Time threshold for clustering

    //void ClusterHits(std::vector<G4Sim::Hit*>& hits, G4double spatialThreshold, G4double timeThreshold, std::vector<Cluster>& clusters);
    void ClusterHits(std::vector<G4Sim::Hit*>& hits, G4double spatialThreshold, G4double timeThreshold, std::vector<Cluster>& clusters, int collectionID);
    G4double CalculateDistance(const G4ThreeVector& pos1, const G4ThreeVector& pos2);
    G4double CalculateTimeDifference(G4double time1, G4double time2);

    // define here all the variables that you want to store for each event in the 
    // ntuple tree  
    G4double fLogWeight;
    G4int fEventID;
    G4int fEventType;
    G4double fXp;
    G4double fYp;
    G4double fZp;

    std::vector<G4double> fE;
    std::vector<G4double> fX;
    std::vector<G4double> fY;
    std::vector<G4double> fZ;
    std::vector<G4double> fW;
    std::vector<G4int> fID;
    
    std::vector<G4double> fEdet;
    std::vector<G4int> fNdet;
    std::vector<G4int> fNphot;
    std::vector<G4int> fNcomp;

    G4bool fFastSimulation = false;
    G4bool fInitializedGraphs = false;
    G4int fNumberOfScattersMax = 0;
    G4int fNumberOfScatters = 0;
    G4bool fHasBeenInXenon = false;

    // maximum energy deposit that is allowed: /run/setMaxEnergy
    G4double fMaxEnergy = 0.0;
    // avaliable energy for scattering energy deposits (relevant for muliple scattering events)
    G4double fAvailableEnergy = 0.0;

    //static std::mutex mtx; // Mutex for thread safety

    std::vector<G4String> fHitsCollectionNames;

    GammaRayHelper* fGammaRayHelper;

    G4int verbosityLevel=0;
    static std::map<G4String, std::pair<G4double, G4double>> fClusteringParameters;

};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


