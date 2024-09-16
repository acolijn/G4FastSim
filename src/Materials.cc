#include "Materials.hh"

#include "Materials.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
//#include "G4State.hh"
#include "nlohmann/json.hpp"
#include <fstream>
#include <iostream>

//#include "G4Element.hh"
#include "G4SystemOfUnits.hh"
//#include <regex>

using json = nlohmann::json;


/**
 * @namespace G4Sim
 * @brief Namespace for the G4Sim library.
/*/
namespace G4Sim {

    // Function to convert string to G4State
G4State StringToState(const std::string& stateStr) {
    if (stateStr == "kStateSolid") return kStateSolid;
    if (stateStr == "kStateLiquid") return kStateLiquid;
    if (stateStr == "kStateGas") return kStateGas;
    return kStateUndefined; // Default if not specified
}

// Function to parse value with units (e.g., "2.862 g/cm3" -> 2.862 * g/cm3)
G4double ParseValueWithUnits(const std::string& valueWithUnits) {
    std::istringstream iss(valueWithUnits);
    G4double value;
    std::string unitStr;
    iss >> value >> unitStr;

    if (unitStr == "g/cm3") return value * g / cm3;
    else if (unitStr == "K") return value * kelvin;
    else if (unitStr == "atm") return value * atmosphere;
    else if (unitStr == "mm") return value * mm;
    else if (unitStr == "cm") return value * cm;
    else if (unitStr == "m") return value * m;
    else if (unitStr == "eV") return value * eV;
    // Add more unit handling as needed

    G4cerr << "Unknown unit: " << unitStr << G4endl;
    return 0.0;  // Default if no match
}

/**
 * @brief Constructs a Materials object with the specified file name.
 * 
 * This constructor initializes a Materials object and sets the material file name.
 * 
 * @param fileName The name of the file containing material definitions.
 */
Materials::Materials(G4String fileName)
    : matFileName(fileName)
{
}

Materials::~Materials()
{
}

/**
 * @brief Defines materials based on a JSON configuration file.
 *
 * This function reads a JSON file that contains material definitions and uses
 * the Geant4 NIST material manager to define these materials. It supports both
 * predefined NIST materials and custom materials with specified components.
 *
 * @param jsonFile The path to the JSON file containing material definitions.
 *
 * The JSON file should have the following structure:
 * {
 *   "materials": {
 *     "MaterialName1": {
 *       "nist": "NIST_Material_Name"
 *     },
 *     "MaterialName2": {
 *       "density": <density_in_g/cm3>,
 *       "state": "kStateSolid" | "kStateLiquid" | "kStateGas",
 *       "temperature": <temperature_in_kelvin>,
 *       "pressure": <pressure_in_atmosphere>,
 *       "components": {
 *         "ElementName1": <fraction>,
 *         "ElementName2": <fraction>,
 *         ...
 *       }
 *     },
 *     ...
 *   }
 * }
 *
 * - "nist": Specifies a predefined NIST material.
 * - "density": The density of the custom material in g/cm3.
 * - "state": The state of the custom material (solid, liquid, or gas).
 * - "temperature": The temperature of the custom material in kelvin.
 * - "pressure": The pressure of the custom material in atmosphere.
 * - "components": A dictionary of elements and their fractions in the custom material.
 *
 * If a material is defined using the "nist" key, it will be created using the
 * predefined NIST material. If a material is defined using the "components" key,
 * it will be created as a custom material with the specified components.
 *
 * Example usage:
 * @code
 * Materials materials;
 * materials.DefineMaterials("path/to/materials.json");
 * @endcode
 */
void Materials::DefineMaterials() {
    G4cout << "Materials::DefineMaterials" << G4endl;

    // Get NIST material manager
    G4NistManager* nist = G4NistManager::Instance();

    // Load the materials configuration from the JSON file
    G4cout << "Materials::DefineMaterials: Loading materials from JSON file: " << matFileName << G4endl;
    std::ifstream inputFile(matFileName);
    json materialJson;
    inputFile >> materialJson;

    for (const auto& [materialName, materialData] : materialJson["materials"].items()) {
        G4cout << "Defining material: " << materialName << G4endl;

        G4Material* material = nullptr;
        // Handle NIST predefined materials
        if (materialData.contains("nist")) {
            material = nist->FindOrBuildMaterial(materialData["nist"].get<std::string>());
        }
        // Handle custom materials
        else if (materialData.contains("components")) {
            // Parse density
            G4double density = materialData.contains("density") ? ParseValueWithUnits(materialData["density"].get<std::string>()) : 1.0 * g / cm3;
            // Parse temperature
            G4double temperature = materialData.contains("temperature") ? ParseValueWithUnits(materialData["temperature"].get<std::string>()) : STP_Temperature;
            // Parse pressure
            G4double pressure = materialData.contains("pressure") ? ParseValueWithUnits(materialData["pressure"].get<std::string>()) : STP_Pressure;
            // Convert state from string to G4State enum
            G4State state = materialData.contains("state") ? StringToState(materialData["state"].get<std::string>()) : kStateUndefined;
            // Create custom material
            material = new G4Material(materialName, density, materialData["components"].size(), state, temperature, pressure);

            // Define elements for custom materials
            for (const auto& [elementName, fraction] : materialData["components"].items()) {
                G4Element* element = nist->FindOrBuildElement(elementName);
                material->AddElement(element, fraction.get<double>());
            }

        } else {
            G4cerr << "Materials::DefineMaterials: Error: Unknown material type for " << materialName << G4endl;
            exit(-1);
        }
        Print(material);

    }

    G4cout << "Materials::DefineMaterials done" << G4endl;
}

void Materials::Print(G4Material* material) {
    G4cout << "   Material: " << material->GetName() << G4endl;
    G4cout << "   Density: " << material->GetDensity() / (g / cm3) << " g/cm3" << G4endl;
    G4cout << "   State: " << material->GetState() << G4endl;
    G4cout << "   Temperature: " << material->GetTemperature() / kelvin << " K" << G4endl;
    G4cout << "   Pressure: " << material->GetPressure() / atmosphere << " atm" << G4endl;
    G4cout << "   Number of elements: " << material->GetNumberOfElements() << G4endl;
    for (size_t i = 0; i < material->GetNumberOfElements(); i++) {
        const G4Element* element = material->GetElement(i);
        G4cout << "   Element: " << element->GetName() << " Fraction: " << material->GetFractionVector()[i] << G4endl;
    }
}

} // namespace G4Sim
