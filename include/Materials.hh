#ifndef __MATERIALS_HH__
#define __MATERIALS_HH__


namespace G4FastSim {

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
        Materials();
        ~Materials();
        void DefineMaterials();
    private:

};

}
#endif