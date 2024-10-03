#include "DetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"

#include "DetectorConstructionMessenger.hh"
#include "Materials.hh"
#include "SensitiveDetector.hh"
#include "EventAction.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"


#include "nlohmann/json.hpp"
#include <fstream>
#include <iostream>

using namespace G4Sim;
using json = nlohmann::json;

/**
 * @namespace G4Sim
 * @brief Namespace for the G4Sim library..
/*/
namespace G4Sim {

/**
 * @brief Constructor for the DetectorConstruction class.
 */
DetectorConstruction::DetectorConstruction() 
  : G4VUserDetectorConstruction() {
    // Load geometry from JSON
    fMessenger = new DetectorConstructionMessenger(this);
}

/**
 * @brief Destructor for the DetectorConstruction class.
 */
DetectorConstruction::~DetectorConstruction() {}

/**
 * @brief Constructs the physical volume of the detector.
 * 
 * This function loads the geometry from a JSON file and constructs the physical volume of the detector.
 * 
 * @return The constructed physical volume of the detector.
 */
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // construct materials
    fMaterials = new Materials(matFileName);
    fMaterials->DefineMaterials();

    // construct geometry
    G4cout << "DetectorConstruction::Construct: Loading geometry from JSON file: " << geoFileName << G4endl;
    LoadGeometryFromJson(geoFileName);
    G4cout << "DetectorConstruction::Construct: Geometry loaded successfully!" << G4endl;

    return fWorldPhysical;
}

/**
 * @brief Sets the geometry file name.
 * 
 * This function sets the name of the JSON file that contains the geometry information.
 * 
 * @param fileName The name of the JSON file.
 */
void DetectorConstruction::SetGeometryFileName(const std::string& fileName) {
    geoFileName = fileName;
    G4cout << "DetectorConstruction::SetGeometryFilename: JSON file name set to: " << geoFileName << G4endl;
}

/**
 * @brief Sets the material file name for the detector construction.
 *
 * This function assigns the provided file name to the material file name
 * used in the detector construction process. It also outputs a message
 * to the console indicating the new file name.
 *
 * @param fileName The name of the JSON file containing material information.
 */
void DetectorConstruction::SetMaterialFileName(const std::string& fileName) {
    matFileName = fileName;
    G4cout << "DetectorConstruction::SetMaterialFilename: JSON file name set to: " << matFileName << G4endl;
}   

/**
 * Loads the geometry from a JSON file.
 * 
 * @param jsonFileName The path to the JSON file containing the geometry information.
 */
void DetectorConstruction::LoadGeometryFromJson(const std::string& geoFileName) {
    std::ifstream inputFile(geoFileName);
    if (!inputFile.is_open()) {
        G4cerr << "DetectorConstruction::LoadGeometryFromJson: Error: Could not open geometry JSON file: " << geoFileName << G4endl;
        exit(-1);    
    }

    json geometryJson;
    inputFile >> geometryJson;

    // First construct the world volume
    G4Material* worldMaterial = G4Material::GetMaterial("G4_AIR");
    G4double worldSize = geometryJson["world"]["size"].get<double>() * m;
    G4Box* worldBox = new G4Box("World", worldSize / 2, worldSize / 2, worldSize / 2);
    fWorldLogical = new G4LogicalVolume(worldBox, worldMaterial, "World");
    fWorldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), fWorldLogical, "World", nullptr, false, 0, fCheckOverlaps);

    logicalVolumeMap["World"] = fWorldLogical;

    // Loop through volumes and set up the sensitive detectors
    std::map<G4String, std::vector<G4String>> detectorsToVolumes; // Map of detector names to volume names
    for (const auto& volume : geometryJson["volumes"]) {
        G4LogicalVolume* logVol = ConstructVolume(volume);
        if (logVol) {
            G4String name = volume["name"].get<std::string>();
            logicalVolumeMap[name] = logVol;

            // If the volume is active, associate it with a sensitive detector
            if (volume.contains("active") && volume["active"].get<bool>()) {
                G4String detectorName = volume.contains("detectorName") 
                                         ? volume["detectorName"].get<std::string>()
                                         : name;  // Use the name of the volume if detectorName is not specified
                
                // Associate volume name with detector name. This may seems a bit weird, but in this way
                // we can have multiple logical volumes (i.e. hit collections) associated with the same detector.
                // This is crucial for the clustering algorithm, when we have an artificial fidcutial volume
                // inside the LXe. This is needed for the accelerated MC only.
                detectorsToVolumes[detectorName].push_back(name); 

                // Define the clusternig parameters......
                G4double spatialThreshold = 10.0 * mm;  // default
                G4double timeThreshold = 100.0 * ns;    // default

                if (volume["clustering"].contains("spatialThreshold")) {
                    spatialThreshold = volume["clustering"]["spatialThreshold"].get<double>() * mm;
                }

                if (volume["clustering"].contains("timeThreshold")) {
                    timeThreshold = volume["clustering"]["timeThreshold"].get<double>() * ns;
                }

                // Store the thresholds in the map
                fClusteringParameters[name] = std::make_pair(spatialThreshold, timeThreshold);
            }

            // Check if multiple placements are needed (via 'repetitions' field)
            if (volume.contains("repetitions")) {
                // Place the volume multiple times
                PlaceMultipleVolumes(volume, logVol);
            } else {
                // Place the volume once
                PlaceSingleVolume(volume, logVol, -1);
            }
        }
    }

    // Now create sensitive detectors for each detector name and assign associated volumes
    for (const auto& [detectorName, volumeNames] : detectorsToVolumes) {
        MakeVolumeSensitive(detectorName, volumeNames);  // Pass detector name and associated volumes
    }

    EventAction::SetClusteringParameters(fClusteringParameters);
}

/**
 * Makes a volume sensitive by assigning a sensitive detector to it.
 * 
 * @param volumeName The name of the volume to make sensitive.
 * @param collectionName The name of the collection associated with the sensitive detector.
 */
void DetectorConstruction::MakeVolumeSensitive(const G4String& detectorName, const std::vector<G4String>& volumeNames) {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();

    // Create a sensitive detector with the detector name
    auto* sensitiveDetector = new SensitiveDetector(detectorName, volumeNames);
    sdManager->AddNewDetector(sensitiveDetector);

    for (const auto& volumeName : volumeNames) {
        G4LogicalVolume* logicalVolume = GetLogicalVolume(volumeName);
        if (logicalVolume) {
            logicalVolume->SetSensitiveDetector(sensitiveDetector);
            G4cout << "Assigned sensitive detector " << detectorName << " to volume: " << volumeName << G4endl;

            // Register the hits collection name with EventAction
            auto* eventAction = const_cast<EventAction*>(dynamic_cast<const EventAction*>(G4RunManager::GetRunManager()->GetUserEventAction()));
            if (eventAction) {
                eventAction->AddHitsCollectionName(volumeName + "Collection");
            }
        } else {
            G4cerr << "Error: Logical volume " << volumeName << " not found!" << G4endl;
        }
    }
}



/**
 * Constructs a G4LogicalVolume based on the provided volume definition.
 *
 * @param volumeDef The JSON object containing the volume definition.
 * @return The constructed G4LogicalVolume, or nullptr if an error occurred.
 */
G4LogicalVolume* DetectorConstruction::ConstructVolume(const json& volumeDef) {
    G4String name = volumeDef["name"].get<std::string>();
    G4Material* material = G4Material::GetMaterial(volumeDef["material"].get<std::string>());
    G4LogicalVolume* logicalVolume = nullptr;

    G4cout << "DetectorConstruction::ConstructVolume: Constructing volume: " << name << G4endl;
    G4cout << "DetectorConstruction::ConstructVolume: Material: " << material->GetName() << G4endl;  

    // if shape is 'union'
    if (volumeDef["shape"].get<std::string>() == "union"){
        G4VSolid* unionSolid = nullptr;
        for (const auto& component : volumeDef["components"]) {
            // create the solid
            G4VSolid* solid = CreateSolid(component);
            // its relative position
            G4ThreeVector position(component["placement"]["x"].get<double>() * mm,
                                   component["placement"]["y"].get<double>() * mm,
                                   component["placement"]["z"].get<double>() * mm);
            // its rotation (if exists)
            G4RotationMatrix* rotation = GetRotationMatrix(component);
            // Call the function that handles rotation
            if (!unionSolid) {
                unionSolid = solid;  // First component
            } else {
                unionSolid = new G4UnionSolid(name, unionSolid, solid, rotation, position);
            }
        }
        logicalVolume = new G4LogicalVolume(unionSolid, material, name);
    // if shape is 'subtraction'
    } else if (volumeDef["shape"].get<std::string>() == "subtraction"){
        G4VSolid* subtractionSolid = nullptr;
        for (const auto& component : volumeDef["components"]) {
            // create the solid
            G4VSolid* solid = CreateSolid(component);
            // its relative position
            G4ThreeVector position(component["placement"]["x"].get<double>() * mm,
                                   component["placement"]["y"].get<double>() * mm,
                                   component["placement"]["z"].get<double>() * mm);
            // its rotation (if exists)
            G4RotationMatrix* rotation = GetRotationMatrix(component);
            // Call the function that handles rotation
            if (!subtractionSolid) {
                subtractionSolid = solid;  // First component
            } else {
                subtractionSolid = new G4SubtractionSolid(name, subtractionSolid, solid, rotation, position);
            }
        }
        logicalVolume = new G4LogicalVolume(subtractionSolid, material, name);	
    // if shape is 'tubs' or 'box'
    } else {
        G4VSolid* solid = CreateSolid(volumeDef);
        logicalVolume = new G4LogicalVolume(solid, material, name);
    }

    // Set the attributes of the logical volume, like visibility, color, transparency, etc.
    SetAttributes(volumeDef, logicalVolume);
 
    return logicalVolume;
}

/**
 * @brief Sets the attributes of a given logical volume based on the provided JSON definition.
 *
 * This function configures the attributes of the specified logical volume using the
 * parameters defined in the JSON object. The attributes may include properties such as
 * material, color, visibility, and other visualization settings.
 *
 * @param volumeDef A JSON object containing the definition of the volume attributes.
 * @param logicalVolume A pointer to the G4LogicalVolume whose attributes are to be set.
 */
void DetectorConstruction::SetAttributes(const json& volumeDef, G4LogicalVolume* logicalVolume) {
    G4VisAttributes* visAttributes = new G4VisAttributes();

    if (volumeDef.contains("visible")) {
        visAttributes->SetVisibility(volumeDef["visible"].get<bool>());
    } else {
        visAttributes->SetVisibility(true);
    }

    if (volumeDef.contains("color")) {
        std::vector<double> color = volumeDef["color"].get<std::vector<double>>();
        visAttributes->SetColour(G4Colour(color[0], color[1], color[2], color[3]));
    } else {
        visAttributes->SetColour(G4Colour(0.5, 0.5, 0.5, 0.5));
    }

    logicalVolume->SetVisAttributes(visAttributes);
}

/**
 * @brief Places a volume inside its parent volume based on the provided JSON definition.
 *
 * This function reads the volume definition from a JSON object and places the volume
 * inside its parent volume. If a parent volume is specified in the JSON definition, 
 * it retrieves the parent volume by name. If no parent is specified, it defaults to 
 * placing the volume inside the world volume. The function also handles the position 
 * and optional rotation of the volume.
 *
 * @param volumeDef A JSON object containing the volume definition, including its name, 
 *                  parent volume (optional), and placement information (position and rotation).
 * @param logicalVolume A pointer to the logical volume to be placed.
 * @return A pointer to the placed physical volume, or nullptr if an error occurs.
 */
G4VPhysicalVolume* DetectorConstruction::PlaceVolume(const json& volumeDef, G4LogicalVolume* logicalVolume) {

    G4String name = volumeDef["name"].get<std::string>();
    G4VPhysicalVolume* physicalVolume = nullptr;

   // Place volume inside its parent
    if (logicalVolume) {
        G4LogicalVolume* parentVolume = fWorldLogical;  // Default to world

        if (volumeDef.contains("parent")) {
            G4String parentName = volumeDef["parent"].get<std::string>();
            parentVolume = GetLogicalVolume(parentName);
            if (!parentVolume) {
                G4cerr << "Error: Parent volume " << parentName << " not found!" << G4endl;
                exit(-1);
            }
        }

        G4ThreeVector position(volumeDef["placement"]["x"].get<double>() * mm, 
                               volumeDef["placement"]["y"].get<double>() * mm, 
                               volumeDef["placement"]["z"].get<double>() * mm);

        // Extract rotation (if exists)
        G4RotationMatrix* rotation = GetRotationMatrix(volumeDef);

        // place the volume
        physicalVolume = new G4PVPlacement(
            rotation,  // Rotation matrix (can be nullptr)
            position,        // Position vector
            logicalVolume,   // Logical volume to place
            name,      // Name of the volume
            parentVolume,    // Parent logical volume
            false,           // No boolean operation
            0,               // Copy number
            fCheckOverlaps   // Overlap checking
        );
    }

    return physicalVolume;
}

/**
 * @brief Creates a G4RotationMatrix based on the rotation parameters provided in the JSON object.
 *
 * This function extracts the rotation parameters (if they exist) from the given JSON object and 
 * creates a G4RotationMatrix using the specified Euler angles (x, y, z). The angles are expected 
 * to be in degrees and will be converted to radians internally.
 *
 * @param volumeDef A JSON object containing the volume definition, including placement and rotation parameters.
 * @return A pointer to a G4RotationMatrix object initialized with the specified rotation angles, or nullptr if no rotation is specified.
 */
G4RotationMatrix* DetectorConstruction::GetRotationMatrix(const json& volumeDef){

    G4RotationMatrix* rotationMatrix = nullptr;

    // Extract rotation (if exists)
    json rotationJson;
    if (volumeDef["placement"].contains("rotation")) {
        rotationJson = volumeDef["placement"]["rotation"];
    }

    if (rotationJson.contains("x") || rotationJson.contains("y") || rotationJson.contains("z")) {
        G4double rotationX = 0.0;
        G4double rotationY = 0.0;
        G4double rotationZ = 0.0;

        if (rotationJson.contains("x")) {
            rotationX = rotationJson["x"].get<double>() * deg;
        }
        if (rotationJson.contains("y")) {
            rotationY = rotationJson["y"].get<double>() * deg;
        }
        if (rotationJson.contains("z")) {
            rotationZ = rotationJson["z"].get<double>() * deg;
        }

        // Create the rotation matrix with the Euler angles
        rotationMatrix = new G4RotationMatrix();
        rotationMatrix->rotateX(rotationX);
        rotationMatrix->rotateY(rotationY);
        rotationMatrix->rotateZ(rotationZ);
    }

    return rotationMatrix;
}

/**
 * @brief Creates a G4VSolid object based on the provided solid definition.
 * 
 * @param solidDef The JSON object containing the solid definition.
 * @return A pointer to the created G4VSolid object. Returns nullptr if the shape is not supported.
 */
G4VSolid* DetectorConstruction::CreateSolid(const json& solidDef) {

    G4String shape = solidDef["shape"].get<std::string>();
    if (shape == "tubs") {
        G4double rMin = solidDef["dimensions"]["rMin"].get<double>() * mm;
        G4double rMax = solidDef["dimensions"]["rMax"].get<double>() * mm;
        G4double z = solidDef["dimensions"]["z"].get<double>() * mm;
        G4double startAngle = solidDef["dimensions"]["startAngle"].get<double>() * deg;
        G4double spanningAngle = solidDef["dimensions"]["spanningAngle"].get<double>() * deg;
        return new G4Tubs("Tubs", rMin, rMax, z / 2, startAngle, spanningAngle);
    } else if (shape == "box"){
        G4double x = solidDef["dimensions"]["x"].get<double>() * mm;
        G4double y = solidDef["dimensions"]["y"].get<double>() * mm;
        G4double z = solidDef["dimensions"]["z"].get<double>() * mm;
        return new G4Box("Box", x / 2, y / 2, z / 2);
    } else if (shape == "sphere"){
        G4double rMin = solidDef["dimensions"]["rMin"].get<double>() * mm;
        G4double rMax = solidDef["dimensions"]["rMax"].get<double>() * mm;
        G4double startPhi = solidDef["dimensions"]["startPhi"].get<double>() * deg;
        G4double endPhi = solidDef["dimensions"]["endPhi"].get<double>() * deg;
        G4double startTheta = solidDef["dimensions"]["startTheta"].get<double>() * deg;
        G4double endTheta = solidDef["dimensions"]["endTheta"].get<double>() * deg;
        return new G4Sphere("Sphere", rMin, rMax, startPhi, endPhi, startTheta, endTheta);
    } else {
        G4cerr << "Error: Unsupported shape: " << shape << G4endl;
        exit(-1);
    }
    // Add more shapes as needed (G4Box, G4Sphere, etc.)

    return nullptr;
}

/**
 * @brief Retrieves the logical volume with the specified name.
 * 
 * @param name The name of the logical volume to retrieve.
 * @return A pointer to the logical volume if found, otherwise nullptr.
 */
G4LogicalVolume* DetectorConstruction::GetLogicalVolume(const G4String& name) {
    if (logicalVolumeMap.find(name) != logicalVolumeMap.end()) {
        return logicalVolumeMap[name];
    }
    return nullptr;
}

}  // namespace G4XamsSim

/**
 * @brief Places multiple instances of a volume within a logical volume based on the provided JSON definition.
 *
 * This function extracts the base position, rotation, and repetition information from the JSON definition
 * and places multiple copies of the specified volume within the given logical volume.
 *
 * @param volumeDef A JSON object containing the volume definition, including name, placement, rotation, and repetition details.
 * @param logicalVolume A pointer to the G4LogicalVolume where the volumes will be placed.
 *
 * The JSON structure for volumeDef should include:
 * - "name": The base name of the volume.
 * - "placement": An object with "x", "y", and "z" coordinates for the base position.
 * - "rotation": (Optional) An object defining the rotation matrix.
 * - "repetitions": An object with:
 *   - "count": The number of repetitions.
 *   - "dx", "dy", "dz": The displacements in x, y, and z directions for each repetition.
 * - "parent": The name of the parent logical volume.
 *
 * The function will create unique names for each instance by appending an index suffix to the base name.
 * It will also store each created physical volume in the physicalVolumeMap using the unique name as the key.
 */
void DetectorConstruction::PlaceMultipleVolumes(const json& volumeDef, G4LogicalVolume* logicalVolume) {
    G4String name = volumeDef["name"].get<std::string>();
    
    // Extract the base position
    G4double baseX = volumeDef["placement"]["x"].get<double>() * mm;
    G4double baseY = volumeDef["placement"]["y"].get<double>() * mm;
    G4double baseZ = volumeDef["placement"]["z"].get<double>() * mm;
    
    // Default rotation matrix (identity)
    G4RotationMatrix* baseRotation = GetRotationMatrix(volumeDef);
    if (!baseRotation) {
        baseRotation = new G4RotationMatrix();  // Identity matrix (no rotation)
    }

    // Get repetitions info
    int count = volumeDef["repetitions"]["count"].get<int>();
    G4double dx = volumeDef["repetitions"]["dx"].get<double>() * mm;
    G4double dy = volumeDef["repetitions"]["dy"].get<double>() * mm;
    G4double dz = volumeDef["repetitions"]["dz"].get<double>() * mm;

    // Place multiple copies
    for (int i = 0; i < count; ++i) {
        G4ThreeVector position(baseX + i * dx, baseY + i * dy, baseZ + i * dz);
        
        // Append a suffix to make the physical volume name unique
        G4String instanceName = name + "_" + std::to_string(i);

        // Place the volume
        G4VPhysicalVolume* physicalVolume = new G4PVPlacement(
            baseRotation,       // Apply the rotation matrix
            position,           // Position
            logicalVolume,      // Logical volume
            instanceName,       // Unique name for each instance
            GetLogicalVolume(volumeDef["parent"].get<std::string>()),  // Parent volume
            false,              // No boolean operation
            i,                  // Copy number
            fCheckOverlaps      // Check for overlaps
        );
        
        physicalVolumeMap[instanceName] = physicalVolume;
    }
}

/**
 * @brief Places a single volume inside its parent volume based on the provided JSON definition.
 *
 * This function creates and places a physical volume inside its parent logical volume.
 * The placement details such as position and rotation are extracted from the JSON definition.
 * If a parent volume is specified in the JSON, it is used; otherwise, the world volume is used as the parent.
 * The placed volume is stored in a map with a unique key.
 *
 * @param volumeDef A JSON object containing the definition of the volume to be placed.
 *                  Expected keys:
 *                  - "name": The name of the volume (string).
 *                  - "parent": The name of the parent volume (optional, string).
 *                  - "placement": An object containing the position coordinates (x, y, z) in millimeters.
 * @param logicalVolume The logical volume to be placed.
 * @param copyNumber An integer representing the unique copy number of the volume instance (default is 0).
 * @return A pointer to the placed G4VPhysicalVolume.
 */
G4VPhysicalVolume* DetectorConstruction::PlaceSingleVolume(const json& volumeDef, G4LogicalVolume* logicalVolume, int copyNumber = 0) {
    G4String name = volumeDef["name"].get<std::string>();
    G4VPhysicalVolume* physicalVolume = nullptr;

    // Place volume inside its parent
    if (logicalVolume) {
        G4LogicalVolume* parentVolume = fWorldLogical;  // Default to world

        if (volumeDef.contains("parent")) {
            G4String parentName = volumeDef["parent"].get<std::string>();
            parentVolume = GetLogicalVolume(parentName);
            if (!parentVolume) {
                G4cerr << "Error: Parent volume " << parentName << " not found!" << G4endl;
                exit(-1);
            }
        }

        G4ThreeVector position(volumeDef["placement"]["x"].get<double>() * mm, 
                               volumeDef["placement"]["y"].get<double>() * mm, 
                               volumeDef["placement"]["z"].get<double>() * mm);

        // Extract rotation (if exists)
        G4RotationMatrix* rotation = GetRotationMatrix(volumeDef);

        // Place the volume

        G4String instanceName = name;
        if (copyNumber >= 0) instanceName = name+ "_" + std::to_string(copyNumber);
        
        physicalVolume = new G4PVPlacement(
            rotation,                    // Rotation matrix (can be nullptr)
            position,                    // Position vector
            logicalVolume,               // Logical volume to place
            instanceName,                // Unique name of the volume instance
            parentVolume,                // Parent logical volume
            false,                       // No boolean operation
            copyNumber,                  // Unique copy number
            fCheckOverlaps               // Overlap checking
        );

        // Store the placed volume in the physicalVolumeMap with unique key
        physicalVolumeMap[instanceName] = physicalVolume;
    }

    return physicalVolume;
}

