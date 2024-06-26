//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "Materials.hh"
#include "GammaRayHelper.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"	
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4AnalysisManager.hh"



namespace G4FastSim
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(GammaRayHelper* helper)
  : G4VUserDetectorConstruction(),
  fGammaRayHelper(helper){

  }

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  //fGammaRayHelper = &GammaRayHelper::Instance();
  //
  // Verification of geometry
  //
  fCheckOverlaps = true;
  //
  // Define materials
  //
  fMaterials = new Materials();
  fMaterials->DefineMaterials();
  // 
  // Construct world volume
  //
  ConstructWorld();
  //
  // Construct water tank
  //
  ConstructWaterTank();
  //
  // Construct outer cryostat and fill with vacuum
  //
  ConstructOuterCryostat();
  //
  // Construct inner cryostat and fill with LXe
  //
  ConstructInnerCryostat();
  //
  // Initialize the gamma ray helper class with all the materials that have been defined
  //

  // WARNING: do not define ay new materials after this point, as the helper class will not be aware of them
  
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  for (size_t i = 0; i < materialTable->size(); ++i) {
      G4Material* material = (*materialTable)[i];
      fGammaRayHelper->Initialize(material);
      G4cout << "Material: " << material->GetName() << " Attenuation length:" <<  fGammaRayHelper->GetAttenuationLength(1.0 * MeV, material) /cm << " cm" << G4endl;

  }

  
  //
  //always return the physical World
  //
  return fWorldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructWorld(){
  //
  // Construct the world volume
  //
  G4double world_sizeXY = 20*m, world_sizeZ = 15*m;
  
  G4Material* world_mat = G4Material::GetMaterial("G4_AIR");
  auto solidWorld = new G4Box("World", 0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);
  fWorldLogical = new G4LogicalVolume(solidWorld, world_mat, "World");

  fWorldPhysical = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    fWorldLogical,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    fCheckOverlaps);                            // overlaps checking

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructWaterTank(){
  //
  // Water tank
  //
  G4double tank_radius = 5*m;
  G4double tank_height = 10*m;

  G4Material* water = G4Material::GetMaterial("G4_WATER");
  auto solidWaterTank = new G4Tubs("WaterTank",           // its name
    0., tank_radius, tank_height / 2., 0. * deg, 360. * deg);  // its size

  fWaterTankLogical = new G4LogicalVolume(solidWaterTank,  // its solid
    water,                                         // its material
    "WaterTank");                                    // its name

  fWaterTankPhysical = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    fWaterTankLogical,        // its logical volume
    "WaterTank",              // its name
    fWorldLogical,            // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps);           // overlaps checking
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructOuterCryostat(){
  //
  // Outer cryostat
  //
  G4double outer_cryostat_radius = 1.3*m;
  G4double outer_cryostat_height = 2.0*m;
  G4double wall_thickness = 1*cm;

  G4Material* stainless_steel = G4Material::GetMaterial("StainlessSteel");
  auto solidOuterCryostat = new G4Tubs("OuterCryostat",           // its name
    0, outer_cryostat_radius, outer_cryostat_height / 2., 0. * deg, 360. * deg);  // its size

  fOuterCryostatLogical = new G4LogicalVolume(solidOuterCryostat,  // its solid
    stainless_steel,                                         // its material
    "OuterCryostat");                                    // its name

  fOuterCryostatPhysical = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    fOuterCryostatLogical,     // its logical volume
    "OuterCryostat",          // its name
    fWaterTankLogical,        // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps);           // overlaps checking

  // Fill the outer cryostat with vacuum
  G4Material* vacuum = G4Material::GetMaterial("Vacuum");
  auto solidVacuum = new G4Tubs("Vacuum",           // its name
    0., outer_cryostat_radius - wall_thickness, outer_cryostat_height / 2. - wall_thickness, 0. * deg, 360. * deg);  // its size  

  fVacuumLogical = new G4LogicalVolume(solidVacuum,  // its solid
    vacuum,                                         // its material
    "Vacuum");                                    // its name

  fVacuumPhysical = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    fVacuumLogical,             // its logical volume
    "Vacuum",                 // its name
    fOuterCryostatLogical,     // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps);           // overlaps checking

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructInnerCryostat(){
  //
  // Construct inner cryostat and fill with LXe
  //
  G4double inner_cryostat_radius = 1.2*m;
  G4double inner_cryostat_height = 1.8*m;
  G4double wall_thickness = 0.5*cm;

  // the vessel
  G4Material* stainless_steel = G4Material::GetMaterial("StainlessSteel");
  auto solidInnerCryostat = new G4Tubs("InnerCryostat",           // its name
    0., inner_cryostat_radius, inner_cryostat_height / 2., 0. * deg, 360. * deg);  // its size

  fInnerCryostatLogical = new G4LogicalVolume(solidInnerCryostat,  // its solid
    stainless_steel,                                         // its material
    "InnerCryostat");                                    // its name

  fInnerCryostatPhysical = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    fInnerCryostatLogical,     // its logical volume
    "InnerCryostat",          // its name
    fVacuumLogical,            // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps);           // overlaps checking

  // the liquid xenon
  G4Material* lXe = G4Material::GetMaterial("LXe");
  auto solidLXe = new G4Tubs("LXe",           // its name
    0., inner_cryostat_radius - wall_thickness, inner_cryostat_height / 2. - wall_thickness, 0. * deg, 360. * deg);  // its size

  fLXeLogical = new G4LogicalVolume(solidLXe,  // its solid
    lXe,                                         // its material
    "LXe");                                    // its name

  fLXePhysical = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    fLXeLogical,               // its logical volume
    "LXe",                    // its name
    fInnerCryostatLogical,     // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps);           // overlaps checking
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
