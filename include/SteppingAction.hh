#ifndef _SteppingAction_h
#define _SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"
#include "GammaRayHelper.hh"
#include "Hit.hh"
#include "SensitiveDetector.hh"
#include <map>
#include <vector>
#include <utility>

class G4LogicalVolume;

/**
 * @namespace G4FastSim
 * @brief A namespace for fast simulation classes in Geant4.
/*/
namespace G4FastSim
{

class EventAction;

/**
 * @class SteppingAction
 * @brief Class responsible for handling the stepping action of particles in the simulation.
 *
 * This class inherits from the G4UserSteppingAction class and implements the UserSteppingAction method.
 * It also provides methods for setting the verbosity level, printing step information, adding hits to collections,
 * performing scattering calculations, and analyzing standard steps.
 *
 * @note The SteppingAction class requires an EventAction object and a GammaRayHelper object to be passed to its constructor.
 * These objects are used for event-level actions and gamma ray calculations, respectively.
 */
class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(EventAction* eventAction, GammaRayHelper* helper);
    ~SteppingAction();// override = default;

    // method from the base class
    void UserSteppingAction(const G4Step*) override;
    void SetVerbosity(G4int level) { verbosityLevel = level; }  
    void Print(const G4Step* step);

  private:
    void AddHitToCollection(Hit *newhit, G4String collectionName);
    G4double DoScatter(const G4Step* step, G4ThreeVector x0);
    void AnalyzeStandardStep(const G4Step* step);

    EventAction* fEventAction;
    G4LogicalVolume* fScoringVolume;
    G4ParticleTable* particleTable;
    GammaRayHelper* fGammaRayHelper;

    std::map<G4String, HitsCollection*> fHitsCollections;
    std::map<G4String, G4int> fHitsCollectionIDs;
    G4bool fHitsCollectionInitialized;

    G4int verbosityLevel=0;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
