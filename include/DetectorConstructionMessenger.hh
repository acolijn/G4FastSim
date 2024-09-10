#ifndef _DETECTORCONSTRUCTIONMESSENGER_HH
#define _DETECTORCONSTRUCTIONMESSENGER_HH

#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

/**
 * @namespace G4Sim
 * @brief A namespace for fast simulation classes in Geant4.
/*/
namespace G4Sim 
{

class DetectorConstruction;

/**
 * @class DetectorConstructionMessenger
 * @brief A class responsible for handling user commands related to the detector construction.
 *
 * This class inherits from G4UImessenger and provides methods to set new values for various parameters of the detector construction.
 * The class is used to interact with the user interface and update the corresponding values in the DetectorConstruction class.
 */
class DetectorConstructionMessenger : public G4UImessenger {
    public:
        DetectorConstructionMessenger(DetectorConstruction* detector);
        ~DetectorConstructionMessenger();

        void SetNewValue(G4UIcommand* command, G4String newValue) override;

    private:
        DetectorConstruction* fDetectorConstruction;

        G4UIcmdWithADoubleAndUnit* fOuterCryostatRadiusCmd;
        G4UIcmdWithADoubleAndUnit* fOuterCryostatHeightCmd;
        G4UIcmdWithADoubleAndUnit* fOuterCryostatWallThicknessCmd;
        G4UIcmdWithADoubleAndUnit* fInnerCryostatRadiusCmd;
        G4UIcmdWithADoubleAndUnit* fInnerCryostatHeightCmd;
        G4UIcmdWithADoubleAndUnit* fInnerCryostatWallThicknessCmd;
        G4UIcmdWithADoubleAndUnit* fFiducialRadiusCmd;
        G4UIcmdWithADoubleAndUnit* fFiducialHeightCmd;
};

}
#endif // DETECTORCONSTRUCTIONMESSENGER_HH
