#include "DetectorConstructionMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

namespace G4FastSim{

DetectorConstructionMessenger::DetectorConstructionMessenger(DetectorConstruction* detector)
    : fDetectorConstruction(detector) {

    G4UIdirectory* dir = new G4UIdirectory("/detector/");
    dir->SetGuidance("UI commands for the detector setup");

    fOuterCryostatRadiusCmd = new G4UIcmdWithADoubleAndUnit("/detector/setOuterCryostatRadius", this);
    fOuterCryostatRadiusCmd->SetGuidance("Set the radius of the outer cryostat");
    fOuterCryostatRadiusCmd->SetParameterName("Radius", false);
    fOuterCryostatRadiusCmd->SetUnitCategory("Length");
    fOuterCryostatRadiusCmd->SetDefaultValue(1.3 * m);

    fOuterCryostatHeightCmd = new G4UIcmdWithADoubleAndUnit("/detector/setOuterCryostatHeight", this);
    fOuterCryostatHeightCmd->SetGuidance("Set the height of the outer cryostat");
    fOuterCryostatHeightCmd->SetParameterName("Height", false);
    fOuterCryostatHeightCmd->SetUnitCategory("Length");
    fOuterCryostatHeightCmd->SetDefaultValue(2.0 * m);

    fOuterCryostatWallThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/setOuterCryostatWallThickness", this);
    fOuterCryostatWallThicknessCmd->SetGuidance("Set the wall thickness of the outer cryostat");
    fOuterCryostatWallThicknessCmd->SetParameterName("Thickness", false);
    fOuterCryostatWallThicknessCmd->SetUnitCategory("Length");
    fOuterCryostatWallThicknessCmd->SetDefaultValue(1 * cm);

    fInnerCryostatRadiusCmd = new G4UIcmdWithADoubleAndUnit("/detector/setInnerCryostatRadius", this);
    fInnerCryostatRadiusCmd->SetGuidance("Set the radius of the inner cryostat");
    fInnerCryostatRadiusCmd->SetParameterName("Radius", false);
    fInnerCryostatRadiusCmd->SetUnitCategory("Length");
    fInnerCryostatRadiusCmd->SetDefaultValue(1.2 * m);

    fInnerCryostatHeightCmd = new G4UIcmdWithADoubleAndUnit("/detector/setInnerCryostatHeight", this);
    fInnerCryostatHeightCmd->SetGuidance("Set the height of the inner cryostat");
    fInnerCryostatHeightCmd->SetParameterName("Height", false);
    fInnerCryostatHeightCmd->SetUnitCategory("Length");
    fInnerCryostatHeightCmd->SetDefaultValue(1.8 * m);

    fInnerCryostatWallThicknessCmd = new G4UIcmdWithADoubleAndUnit("/detector/setInnerCryostatWallThickness", this);
    fInnerCryostatWallThicknessCmd->SetGuidance("Set the wall thickness of the inner cryostat");
    fInnerCryostatWallThicknessCmd->SetParameterName("Thickness", false);
    fInnerCryostatWallThicknessCmd->SetUnitCategory("Length");
    fInnerCryostatWallThicknessCmd->SetDefaultValue(0.5 * cm);

    fFiducialRadiusCmd = new G4UIcmdWithADoubleAndUnit("/detector/setFiducialRadius", this);
    fFiducialRadiusCmd->SetGuidance("Set the radius of the fiducial volume");
    fFiducialRadiusCmd->SetParameterName("Radius", false);
    fFiducialRadiusCmd->SetUnitCategory("Length");
    fFiducialRadiusCmd->SetDefaultValue(0.8 * m);

    fFiducialHeightCmd = new G4UIcmdWithADoubleAndUnit("/detector/setFiducialHeight", this);
    fFiducialHeightCmd->SetGuidance("Set the height of the fiducial volume");
    fFiducialHeightCmd->SetParameterName("Height", false);
    fFiducialHeightCmd->SetUnitCategory("Length");
    fFiducialHeightCmd->SetDefaultValue(1.0 * m);
}

DetectorConstructionMessenger::~DetectorConstructionMessenger() {
    delete fOuterCryostatRadiusCmd;
    delete fOuterCryostatHeightCmd;
    delete fOuterCryostatWallThicknessCmd;
    delete fInnerCryostatRadiusCmd;
    delete fInnerCryostatHeightCmd;
    delete fInnerCryostatWallThicknessCmd;
    delete fFiducialRadiusCmd;
    delete fFiducialHeightCmd;
}

/**
 * @brief Sets the new value for a given command.
 *
 * This function is called when a new value is assigned to a command in the user interface.
 * It checks which command is being modified and calls the corresponding setter function in the `DetectorConstruction` class.
 *
 * @param command The command being modified.
 * @param newValue The new value assigned to the command.
 */
void DetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {

    if (command == fOuterCryostatRadiusCmd) {
        fDetectorConstruction->SetOuterCryostatRadius(fOuterCryostatRadiusCmd->GetNewDoubleValue(newValue));
    } else if (command == fOuterCryostatHeightCmd) {
        fDetectorConstruction->SetOuterCryostatHeight(fOuterCryostatHeightCmd->GetNewDoubleValue(newValue));
    } else if (command == fOuterCryostatWallThicknessCmd) {
        fDetectorConstruction->SetOuterCryostatWallThickness(fOuterCryostatWallThicknessCmd->GetNewDoubleValue(newValue));
    } else if (command == fInnerCryostatRadiusCmd) {
        fDetectorConstruction->SetInnerCryostatRadius(fInnerCryostatRadiusCmd->GetNewDoubleValue(newValue));
    } else if (command == fInnerCryostatHeightCmd) {
        fDetectorConstruction->SetInnerCryostatHeight(fInnerCryostatHeightCmd->GetNewDoubleValue(newValue));
    } else if (command == fInnerCryostatWallThicknessCmd) {
        fDetectorConstruction->SetInnerCryostatWallThickness(fInnerCryostatWallThicknessCmd->GetNewDoubleValue(newValue));
    } else if (command == fFiducialRadiusCmd) {
        fDetectorConstruction->SetFiducialRadius(fFiducialRadiusCmd->GetNewDoubleValue(newValue));
    } else if (command == fFiducialHeightCmd) {
        fDetectorConstruction->SetFiducialHeight(fFiducialHeightCmd->GetNewDoubleValue(newValue));
    }
}

}