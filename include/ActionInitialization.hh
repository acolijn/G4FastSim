#ifndef B1ActionInitialization_h
#define B1ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "GammaRayHelper.hh"
#include "EventAction.hh"

/**
 * @namespace G4FastSim
 * @brief A namespace for fast simulation classes in Geant4.
/*/
namespace G4FastSim
{

/**
 * @class ActionInitialization
 * @brief Initializes the actions for the simulation.
 *
 * This class is responsible for initializing the actions required for the simulation.
 * It provides methods to build the actions for the master thread and for worker threads.
 * The actions are built based on a GammaRayHelper object.
 */
class ActionInitialization : public G4VUserActionInitialization
{
  public:
    ActionInitialization(GammaRayHelper* helper);
    ~ActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;
  private:
    GammaRayHelper* fGammaRayHelper;
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
