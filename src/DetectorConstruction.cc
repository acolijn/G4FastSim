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
    //std::map<G4String, std::vector<G4String>> detectorsToVolumes; // Map of detector names to volume names
    for (const auto& volume : geometryJson["volumes"]) {
        // Define the volumes. Assembly volumes area special case.....
        G4cout << "DetectorConstruction::LoadGeometryFromJson: Constructing volume: " << volume["name"].get<std::string>() << G4endl;

        if (volume["shape"] == "assembly") {
            // Handle assembly volume... (the logical volume is already made active inside the assembly)
            G4AssemblyVolume* assembly = ConstructAssembly(volume);
            // Place the assembly in the world or parent volume
            PlaceAssembly(volume, assembly);
        } else {
            // Construct the logical volume
            G4LogicalVolume* logVol = ConstructVolume(volume);
            // Place volumes 
            if (volume.contains("repetitions")) {
                PlaceRepeatedVolume(volume, logVol);
            } else {
                // place all incarnations of the volume
                PlaceVolume(volume, logVol);
            } // end if volume
        } // end if assembly
    } // end for volume

    // Set the clustering parameters in the EventAction
    EventAction::SetClusteringParameters(fClusteringParameters);
}

/**
 * @brief Places an assembly volume in the specified parent volume based on the provided JSON configuration.
 * 
 * This function iterates over the "placement" array in the JSON object to position and optionally rotate 
 * the assembly volume within the parent volume. The position is specified in millimeters, and the rotation 
 * is applied if provided.
 * 
 * @param volume A JSON object containing the configuration for the assembly placement. It includes:
 *               - "placement": An array of objects, each specifying:
 *                   - "x": The x-coordinate of the position in millimeters.
 *                   - "y": The y-coordinate of the position in millimeters.
 *                   - "z": The z-coordinate of the position in millimeters.
 *                   - "rotation" (optional): An object specifying the rotation matrix.
 *               - "parent": A string specifying the name of the parent logical volume.
 * @param assembly A pointer to the G4AssemblyVolume object that will be placed in the parent volume.
 */
void DetectorConstruction::PlaceAssembly(const json& volume, G4AssemblyVolume* assembly) {

    // Place the assembly in the world or parent volume
    for (const auto& placement : volume["placement"]) {
        G4ThreeVector position(placement["x"].get<double>() * mm, 
                            placement["y"].get<double>() * mm, 
                            placement["z"].get<double>() * mm);
        G4RotationMatrix* rotation = nullptr;
        if (placement.contains("rotation")) {
            rotation = GetRotationMatrix(placement["rotation"]);
        }

        G4LogicalVolume* parentVolume = GetLogicalVolume(volume["parent"].get<std::string>());
        assembly->MakeImprint(parentVolume, position, rotation);  // Place the assembly
    }
}

/**
 * @brief Constructs an assembly volume from a JSON definition.
 *
 * This function creates a new G4AssemblyVolume and iterates over the components
 * defined in the provided JSON object. For each component, it constructs a logical
 * volume and places it within the assembly at the specified position and rotation.
 *
 * @param volumeDef A JSON object defining the components of the assembly. Each component
 *                  should include placement information with x, y, z coordinates and
 *                  optionally a rotation matrix.
 * @return A pointer to the constructed G4AssemblyVolume.
 *
 * @note If a logical volume cannot be constructed for a component, the function will
 *       output an error message and terminate the program.
 */
G4AssemblyVolume* DetectorConstruction::ConstructAssembly(const json& volumeDef) {
    G4cout << "DetectorConstruction::ConstructAssembly: Constructing assembly volume: " << volumeDef["name"].get<std::string>() << G4endl;
    G4AssemblyVolume* assembly = new G4AssemblyVolume();

    // Loop over components and construct the logical volumes...
    // This is diffent from the union and subtraction volumes, since there we glue the components together.
    for (const auto& component : volumeDef["components"]) {
        G4LogicalVolume* logVol = ConstructVolume(component);
        
        if(logVol){
            G4String name = component["name"].get<std::string>();
            logicalVolumeMap[name] = logVol;
            // Create the position vector
            G4ThreeVector position(component["placement"][0]["x"].get<double>() * mm,
                                component["placement"][0]["y"].get<double>() * mm,
                                component["placement"][0]["z"].get<double>() * mm);
            // Create the rotation matrix
            G4RotationMatrix* rotation = nullptr;
            if (component["placement"][0].contains("rotation")) {
                rotation = GetRotationMatrix(component["placement"][0]["rotation"]);
            }
            // Add the logical volume to the assembly
            assembly->AddPlacedVolume(logVol, position, rotation);
        } else {
            G4cerr << "Error: Logical volume " << volumeDef["name"].get<std::string>() << " not found!" << G4endl;
            exit(-1);
        }
    }
    return assembly;
}

/**
 * Makes a volume sensitive by assigning a sensitive detector to it.
 * 
 * @param volumeName The name of the volume to make sensitive.
 * @param collectionName The name of the collection associated with the sensitive detector.
 */
void DetectorConstruction::MakeVolumeSensitive(const G4String& detectorName, const G4String& volumeName) {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    // Create a sensitive detector with the detector name
    // check if the sensitive detector already exists. the sensitive detector is created only once
    // if it exists get the pointer to it otherwise create a new one
    auto* sensitiveDetector = dynamic_cast<SensitiveDetector*>(sdManager->FindSensitiveDetector(detectorName));
    if (!sensitiveDetector) {
        G4cout << "DetectorConstruction::MakeVolumeSensitive: Creating new sensitive detector: " << detectorName << G4endl;
        sensitiveDetector = new SensitiveDetector(detectorName, volumeName);
        sdManager->AddNewDetector(sensitiveDetector);
    } else {
        G4cout << "DetectorConstruction::MakeVolumeSensitive: Using existing sensitive detector: " << detectorName << G4endl;
        // just add the collection....
        sensitiveDetector->AddCollection(volumeName);
    }

    // Assign the sensitive detector to the logical volume
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

    // Handle "union" or "subtraction"
    if (volumeDef["shape"] == "union" || volumeDef["shape"] == "subtraction") {
        G4VSolid* compoundSolid = nullptr;
        for (const auto& component : volumeDef["components"]) {
            G4VSolid* solid = CreateSolid(component);
            G4ThreeVector position(component["placement"][0]["x"].get<double>() * mm,
                                   component["placement"][0]["y"].get<double>() * mm,
                                   component["placement"][0]["z"].get<double>() * mm);
            G4RotationMatrix* rotation = nullptr;
            if (component["placement"][0].contains("rotation")) {
                rotation = GetRotationMatrix(component["placement"][0]["rotation"]);
            }

            if (!compoundSolid) {
                compoundSolid = solid;  // First component
            } else if (volumeDef["shape"] == "union") {
                compoundSolid = new G4UnionSolid(name, compoundSolid, solid, rotation, position);
            } else if (volumeDef["shape"] == "subtraction") {
                compoundSolid = new G4SubtractionSolid(name, compoundSolid, solid, rotation, position);
            }
        }
        logicalVolume = new G4LogicalVolume(compoundSolid, material, name);
    } else {
        // Simple solids like "tubs", "box", etc.
        G4VSolid* solid = CreateSolid(volumeDef);
        logicalVolume = new G4LogicalVolume(solid, material, name);
    }
    logicalVolumeMap[name] = logicalVolume;

    // Handle active volumes and detectors
    if (volumeDef.contains("active") && volumeDef["active"].get<bool>()) {
        G4String detectorName = volumeDef.contains("detectorName") 
                                ? volumeDef["detectorName"].get<std::string>()
                                : name;
        // make the volume sensitive
        MakeVolumeSensitive(detectorName, name);  // Pass detector name and associated volumes
        // get clustering parameters
        G4double spatialThreshold = volumeDef["clustering"].contains("spatialThreshold") ? 
                                    volumeDef["clustering"]["spatialThreshold"].get<double>() * mm : 10.0 * mm;
        G4double timeThreshold = volumeDef["clustering"].contains("timeThreshold") ? 
                                volumeDef["clustering"]["timeThreshold"].get<double>() * ns : 100.0 * ns;
        // store for later use....
        fClusteringParameters[name] = std::make_pair(spatialThreshold, timeThreshold);
    }

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
 * @brief Creates a G4RotationMatrix based on the rotation parameters provided in the JSON object.
 *
 * This function extracts the rotation parameters (if they exist) from the given JSON object and 
 * creates a G4RotationMatrix using the specified Euler angles (x, y, z). The angles are expected 
 * to be in degrees and will be converted to radians internally.
 *
 * @param volumeDef A JSON object containing the volume definition, including placement and rotation parameters.
 * @return A pointer to a G4RotationMatrix object initialized with the specified rotation angles, or nullptr if no rotation is specified.
 */
//G4RotationMatrix* DetectorConstruction::GetRotationMatrix(const json& volumeDef){
G4RotationMatrix* DetectorConstruction::GetRotationMatrix(const json& rotationJson) {

    G4RotationMatrix* rotationMatrix = nullptr;

    // Extract rotation (if exists)
    //json rotationJson;
    //if (placementDef["placement"].contains("rotation")) {
    //     rotationJson = placementDef["placement"]["rotation"];
    //}

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
void DetectorConstruction::PlaceRepeatedVolume(const json& volumeDef, G4LogicalVolume* logicalVolume) {
    G4String name = volumeDef["name"].get<std::string>();
    
    // Extract the base position
    G4double baseX = volumeDef["placement"]["x"].get<double>() * mm;
    G4double baseY = volumeDef["placement"]["y"].get<double>() * mm;
    G4double baseZ = volumeDef["placement"]["z"].get<double>() * mm;
    
    // Default rotation matrix (identity)
    G4RotationMatrix* baseRotation = nullptr;
    if (volumeDef["placement"].contains("rotation")) {
        baseRotation = GetRotationMatrix(volumeDef["placement"]["rotation"]);
    } else {
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
 * @brief Places a volume within a parent logical volume based on the provided JSON definition.
 *
 * This function reads the volume definition from a JSON object, retrieves the parent logical volume,
 * and places the volume at specified positions and rotations within the parent volume. It handles
 * multiple placements and assigns unique instance names to each placed volume.
 *
 * @param volumeDef A JSON object containing the volume definition, including name, parent, and placement details.
 * @param logicalVolume A pointer to the logical volume to be placed.
 *
 * The JSON object `volumeDef` should have the following structure:
 * {
 *   "name": "volumeName",
 *   "parent": "parentVolumeName",
 *   "placement": [
 *     {
 *       "x": <double>, "y": <double>, "z": <double>,
 *       "rotation": <optional rotation definition>
 *     },
 *     ...
 *   ]
 * }
 *
 * The function will place the volume at each specified position and rotation, and store the placed
 * volumes in the `physicalVolumeMap` with unique instance names.
 *
 * @throws std::runtime_error if the parent volume is not found.
 */
void DetectorConstruction::PlaceVolume(const json& volumeDef, G4LogicalVolume* logicalVolume){//}, const json& placement, int copyNumber) {
    G4String name = volumeDef["name"].get<std::string>();
    G4LogicalVolume* parentVolume = GetLogicalVolume(volumeDef["parent"].get<std::string>());

    // Check if the parent volume exists
    if (!parentVolume) {
        G4cerr << "Error: Parent volume " << volumeDef["parent"].get<std::string>() << " not found!" << G4endl;
        exit(-1);
    }
    // loop over placement and count copyNumber
    int copyNumber = 0;
    for (const auto& placement : volumeDef["placement"]) {
        // Extract position and rotation
        G4ThreeVector position(placement["x"].get<double>() * mm, placement["y"].get<double>() * mm, placement["z"].get<double>() * mm);
        G4RotationMatrix* rotation = nullptr;
        if (placement.contains("rotation")) {
            rotation = GetRotationMatrix(placement["rotation"]);
        }

        // Place the volume
        G4String instanceName = name;
        if (copyNumber == 0){
            instanceName = name;
        } else {
            instanceName = name+ "_" + std::to_string(copyNumber);
        }
        G4VPhysicalVolume* physicalVolume = new G4PVPlacement(
            rotation, 
            position, 
            logicalVolume, 
            instanceName, 
            parentVolume, 
            false, 
            copyNumber, 
            fCheckOverlaps);

        // Store the placed volume in the physicalVolumeMap with unique key
        physicalVolumeMap[instanceName] = physicalVolume;
        copyNumber++;
    } // end for placement

    //return physicalVolume;
}
