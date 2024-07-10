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

/// Event action class
///

class G4Event;
class G4HCofThisEvent;
class G4VHitsCollection;

namespace G4FastSim
{



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
    std::vector<G4double>& GetE(){return fEd;};

    void AddEdep(G4double edep) { fEdep += edep; }
    void AddWeight(G4double weight) { fWeight *= weight; }
    void AnalyzeHits(const G4Event* event);
    void ResetVariables();

    G4bool IsFastSimulation() { return fFastSimulation; }
    G4int GetNumberOfScatters() { return fEd.size();}
    G4double GetMaxEnergy() { return fMaxEnergy; }
    G4int GetNumberOfScattersMax() { return fNumberOfScattersMax; }

    void SetNumberOfScattersMax(G4int n) { fNumberOfScattersMax = n; }
    void SetFastSimulation(G4bool fast) { fFastSimulation = fast; }
    void SetNumberOfScatters(G4int n) { fNumberOfScatters = n; }
    void SetMaxEnergy(G4double e) { fMaxEnergy = e; }



  private:
    //
    // functions for hit clustering
    //
    void ClusterHits(std::vector<G4FastSim::Hit*>& hits, G4double spatialThreshold, G4double timeThreshold, std::vector<Cluster>& clusters);
    G4double CalculateDistance(const G4ThreeVector& pos1, const G4ThreeVector& pos2);
    G4double CalculateTimeDifference(G4double time1, G4double time2);

    // define here all the variables that you want to store for each event in the 
    // ntuple tree  
    G4double fWeight;
    G4double fEdep;
    G4double fXp;
    G4double fYp;
    G4double fZp;

    std::vector<G4double> fEd;
    std::vector<G4double> fX;
    std::vector<G4double> fY;
    std::vector<G4double> fZ;

    G4bool fFastSimulation = false;
    G4bool fInitializedGraphs = false;
    G4int fNumberOfScattersMax = 0;
    G4int fNumberOfScatters = 0;
    G4double fMaxEnergy = 0.0;

    //static std::mutex mtx; // Mutex for thread safety

    std::vector<G4String> fHitsCollectionNames;

    GammaRayHelper* fGammaRayHelper;


};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


