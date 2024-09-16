#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSolid.hh"
#include "DetectorConstructionMessenger.hh"
#include "Materials.hh"

#include "nlohmann/json.hpp"

#include <map>
#include <string>

/**
 * @namespace G4Sim
 * @brief Namespace for the G4Sim library.
/*/
namespace G4Sim {

/**
 * @class DetectorConstruction
 * @brief Class for constructing the detector geometry.
 *
 * This class inherits from G4VUserDetectorConstruction and is responsible for constructing the detector geometry.
 * It defines methods for constructing the physical volume, setting the JSON geometry file name, loading the geometry from a JSON file,
 * making a volume sensitive, constructing a logical volume, creating a solid, and getting a logical volume.
 * The class also contains maps to store logical and physical volumes for easy lookup, as well as the world physical and logical volumes.
 * It has a flag to check for overlaps and a pointer to the Materials class.
 * The class also has a member variable to store the JSON file name and a pointer to the DetectorConstructionMessenger class.
 */
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
        
    // set the JSON geometry file name
    void SetGeometryFileName(const std::string& fileName);
    void SetMaterialFileName(const std::string& fileName);

private:
    G4VPhysicalVolume* PlaceVolume(const nlohmann::json& volumeDef, G4LogicalVolume* logicalVolume);
    G4RotationMatrix* GetRotationMatrix(const nlohmann::json& volumeDef);
    void SetAttributes(const nlohmann::json& volumeDef, G4LogicalVolume* logicalVolume);
    void LoadGeometryFromJson(const std::string& jsonFileName);
    void MakeVolumeSensitive(const G4String& volumeName, const G4String& collectionName);
    G4LogicalVolume* ConstructVolume(const nlohmann::json& volumeDef);
    G4VSolid* CreateSolid(const nlohmann::json& solidDef);
    G4LogicalVolume* GetLogicalVolume(const G4String& name);

    // Maps to store logical and physical volumes for easy lookup
    std::map<G4String, G4LogicalVolume*> logicalVolumeMap;
    std::map<G4String, G4VPhysicalVolume*> physicalVolumeMap;

    // World physical and logical volumes
    G4LogicalVolume* fWorldLogical = nullptr;
    G4VPhysicalVolume* fWorldPhysical = nullptr;

    // Flag to check overlaps
    G4bool fCheckOverlaps = true;
    Materials *fMaterials = nullptr;

    std::string geoFileName;
    std::string matFileName;
    
    std::map<G4String, std::pair<G4double, G4double>> fClusteringParameters;

    DetectorConstructionMessenger* fMessenger;
};

}  // namespace G4XamsSim

#endif
