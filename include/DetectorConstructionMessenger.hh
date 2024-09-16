#ifndef _DETECTORCONSTRUCTIONMESSENGER_HH
#define _DETECTORCONSTRUCTIONMESSENGER_HH

#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

/**
 * @namespace G4Sim
 * @brief Namespace for the G4Sim library.
/*/
namespace G4Sim {

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
        G4UIcmdWithAString* fGeometryFileNameCmd;  // New command to set geometry file name
        G4UIcmdWithAString* fMaterialFileNameCmd;  // New command to set material file name

};

}
#endif // DETECTORCONSTRUCTIONMESSENGER_HH
