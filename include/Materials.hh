#ifndef __MATERIALS_HH__
#define __MATERIALS_HH__

#include "G4String.hh"
#include "G4Material.hh"

/**
 * @namespace G4Sim
 * @brief Namespace for the G4Sim library.
/*/
namespace G4Sim {

/**
 * @class Materials
 * @brief A class for defining materials using the G4NistManager.
 *
 * The Materials class provides a way to define materials using the G4NistManager, which is a G4MaterialDatabase.
 * It is a singleton class, meaning that there is only one instance of it in the program.
 * The G4NistManager can be accessed by calling the static method G4NistManager::Instance().
 * The class provides a method called DefineMaterials() to define materials.
 */
class Materials {
    //
    // Materials
    //
    // Materials are defined in the G4NistManager, which is a G4MaterialDatabase.
    // The G4NistManager provides a list of predefined materials from the NIST database.
    // The G4NistManager is a singleton class, which means that there is only one instance of it in the program.
    // The G4NistManager can be accessed by calling the static method G4NistManager::Instance().
    // The G4NistManager provides a method called FindOrBuildMaterial() that returns a pointer to a G4Material object.
    // The FindOrBuildMaterial() method takes a string as an argument, which is the name of the material.

    public:
        Materials(G4String fileName);
        ~Materials();
        void DefineMaterials();
        void Print(G4Material* material);
    private:

        G4String matFileName;

};

}
#endif