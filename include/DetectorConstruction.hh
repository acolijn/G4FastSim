#ifndef __DetectorConstruction__
#define __DetectorConstruction__ 1

#include "G4VUserDetectorConstruction.hh"
#include "DetectorConstructionMessenger.hh"
#include "GammaRayHelper.hh"
#include "Materials.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class DetectorConstructionMessenger;

/**
 * @namespace G4FastSim
 * @brief A namespace for fast simulation classes in Geant4.
/*/
namespace G4FastSim
{

/**
 * @class DetectorConstruction
 * @brief Class responsible for constructing the detector geometry.
 *
 * This class inherits from G4VUserDetectorConstruction and is used to define the geometry of the detector.
 * It constructs various volumes such as the world volume, water tank, outer cryostat, inner cryostat, LXe, and fiducial volume.
 * It also defines a sensitive detector and provides methods to access the scoring volume and the GammaRayHelper object.
 * The dimensions of the various volumes can be set using the provided setter methods.
 * The class also handles checking for overlaps between volumes.
 */
class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction(G4FastSim::GammaRayHelper *gammaRayHelper);
    ~DetectorConstruction();//override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructWorld();
    void ConstructWaterTank();
    void ConstructOuterCryostat();
    void ConstructInnerCryostat();    
    void ConstructLXe();
    void ConstructFiducialVolume();
    void DefineSensitiveDetector();

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    GammaRayHelper* GetGammaRayHelper() const { return fGammaRayHelper; }

    void SetOuterCryostatRadius(G4double radius) { outer_cryostat_radius = radius; }
    void SetOuterCryostatHeight(G4double height) { outer_cryostat_height = height; }
    void SetOuterCryostatWallThickness(G4double thickness) { outer_cryostat_wall_thickness = thickness; }
    void SetInnerCryostatRadius(G4double radius) { inner_cryostat_radius = radius; }
    void SetInnerCryostatHeight(G4double height) { inner_cryostat_height = height; }
    void SetInnerCryostatWallThickness(G4double thickness) { inner_cryostat_wall_thickness = thickness; }
    void SetFiducialRadius(G4double radius) { fiducial_radius = radius; }
    void SetFiducialHeight(G4double height) { fiducial_height = height; }

  private:
    Materials *fMaterials = nullptr;

    GammaRayHelper *fGammaRayHelper;

    // check for overlaps
    G4bool fCheckOverlaps = true;

    //dimensions
    G4double outer_cryostat_radius;
    G4double outer_cryostat_height;
    G4double outer_cryostat_wall_thickness;

    G4double inner_cryostat_radius;
    G4double inner_cryostat_height;
    G4double inner_cryostat_wall_thickness;

    G4double fiducial_radius;
    G4double fiducial_height;

    //logical volumes
    G4LogicalVolume* fWorldLogical = nullptr;
    G4LogicalVolume* fWaterTankLogical = nullptr;
    G4LogicalVolume* fOuterCryostatLogical = nullptr;
    G4LogicalVolume* fVacuumLogical = nullptr;
    G4LogicalVolume* fInnerCryostatLogical = nullptr;
    G4LogicalVolume* fLXeLogical = nullptr;
    G4LogicalVolume* fPTFELogical = nullptr;
    G4LogicalVolume* fLXeFiducialLogical = nullptr;

    //physical volumes
    G4VPhysicalVolume* fWorldPhysical = nullptr;
    G4VPhysicalVolume* fWaterTankPhysical = nullptr;
    G4VPhysicalVolume* fOuterCryostatPhysical = nullptr;
    G4VPhysicalVolume* fVacuumPhysical = nullptr;
    G4VPhysicalVolume* fInnerCryostatPhysical = nullptr;
    G4VPhysicalVolume* fLXePhysical = nullptr;
    G4VPhysicalVolume* fPTFEPhysical = nullptr;
    G4VPhysicalVolume* fLXeFiducialPhysical = nullptr;

  protected:
    G4LogicalVolume* fScoringVolume = nullptr;

    DetectorConstructionMessenger* fMessenger;

};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
