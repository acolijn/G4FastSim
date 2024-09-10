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
 * @namespace G4Sim
 * @brief A namespace for fast simulation classes in Geant4.
/*/
namespace G4Sim
{

class EventAction;

/**
 * @class SteppingAction
 * @brief Class responsible for handling the stepping action of particles in the simulation.
 *
 * The SteppingAction class handles particle stepping in the simulation, recording hits,
 * calculating scattering, and managing verbosity for diagnostic purposes.
 *
 * @details
 * This class inherits from G4UserSteppingAction and overrides the UserSteppingAction method.
 * It also provides additional private methods for adding hits, performing scatter calculations,
 * and analyzing the steps.
 *
 * @note The class requires an EventAction object and a GammaRayHelper object in its constructor.
 */
class SteppingAction : public G4UserSteppingAction
{
  public:
    /**
     * @brief Constructor for the SteppingAction class.
     * 
     * @param eventAction Pointer to the EventAction object for event-level actions.
     * @param helper Pointer to the GammaRayHelper object for gamma ray-related calculations.
     */
    SteppingAction(EventAction* eventAction, GammaRayHelper* helper);

    /**
     * @brief Destructor for the SteppingAction class.
     */
    ~SteppingAction(); // override = default;

    /**
     * @brief Method that is called at each step of a particle's trajectory.
     * 
     * @param step Pointer to the current G4Step object.
     */
    void UserSteppingAction(const G4Step*) override;

    /**
     * @brief Sets the verbosity level for debugging output.
     * 
     * @param level Integer representing the verbosity level.
     */
    void SetVerbosity(G4int level) { verbosityLevel = level; }

    /**
     * @brief Prints step information for diagnostic purposes.
     * 
     * @param step Pointer to the current G4Step object.
     */
    void Print(const G4Step* step);

  private:
    /* 
     * @brief Adds a hit to the specified hit collection.
     * 
     * @param newhit Pointer to the new hit object to be added.
     * @param collectionName Name of the hit collection where the hit will be added.
     */
    void AddHitToCollection(Hit *newhit, G4String collectionName);

    /**
     * @brief Performs scattering calculations for a given step.
     * 
     * @param step Pointer to the current G4Step object.
     * @param x0 Initial position of the particle.
     * @return G4double The result of the scatter calculation.
     */
    G4double DoScatter(const G4Step* step, G4ThreeVector x0);

    /**
     * @brief Analyzes the standard step for particle interaction.
     * 
     * @param step Pointer to the current G4Step object.
     */
    void AnalyzeStandardStep(const G4Step* step);

    /// Pointer to the EventAction object.
    EventAction* fEventAction;

    /// Pointer to the logical volume used for scoring.
    G4LogicalVolume* fScoringVolume;

    /// Pointer to the G4ParticleTable, used for particle information.
    G4ParticleTable* particleTable;

    /// Pointer to the GammaRayHelper object for gamma ray calculations.
    GammaRayHelper* fGammaRayHelper;

    /// Map of hit collections by name.
    std::map<G4String, HitsCollection*> fHitsCollections;

    /// Map of hit collection IDs by name.
    std::map<G4String, G4int> fHitsCollectionIDs;

    /// Boolean flag indicating if the hit collections have been initialized.
    G4bool fHitsCollectionInitialized;

    /// Verbosity level for debug output.
    G4int verbosityLevel = 0;
};

}

#endif
