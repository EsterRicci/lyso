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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fTargetMater(0), fLogicTarget(0),
 fDetectorMater(0), fLogicDetector(0), 
 fWorldMater(0), fPhysiWorld(0),
 fAuxiliaryDetectorMater(0),fLogicAuxiliaryDetector(0),
 fDetectorMessenger(0)
{
  fTargetLength      = 2*mm; 
  fTargetRadius      = 2*mm;
  fDetectorLength    = 2.0*cm; 
  fDetectorThickness = 2.5*cm;

  fAuxiliaryDetectorRadius=4*cm;
  fAuxiliaryDetectorLength=10*cm;
  
  
  fWorldLength = 30*cm;//std::max(fTargetLength,fDetectorLength);
  fWorldRadius = 30*cm;//fTargetRadius + fDetectorThickness;
      
  DefineMaterials();
    
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // build materials
  //
  G4NistManager* man = G4NistManager::Instance();

  G4Element* Lutetium = man->FindOrBuildElement("Lu");
  G4Element* Yttrium = man->FindOrBuildElement("Y");
  G4Element* Oxygen = man->FindOrBuildElement("O");
  G4Element* Cerium = man->FindOrBuildElement("Ce");
  G4Element* Silicon = man->FindOrBuildElement("Si");

  G4Material* LYSO = new G4Material("LYSO",7.1*g/cm3,(G4int) 5);
  LYSO->AddElement(Lutetium,71.43*perCent);
  LYSO->AddElement(Yttrium,4.03*perCent);
  LYSO->AddElement(Silicon,6.37*perCent);
  LYSO->AddElement(Oxygen,18.14*perCent);
  LYSO->AddElement(Cerium,0.02*perCent); 
  fDetectorMater = LYSO;
   // new G4Material("Germanium", 32, 72.61*g/mole, 5.323*g/cm3);
  

  G4Element* N  = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",   "O", 8, 16.00*g/mole);
  //
  G4int ncomponents; G4double fractionmass;      
  G4Material* Air20 = new G4Material("Air", 1.205e-12*mg/cm3, ncomponents=2,
                      kStateGas, 293.*kelvin, 1*atmosphere);
    Air20->AddElement(N, fractionmass=0.7);
    Air20->AddElement(O, fractionmass=0.3);

  G4Material* Vacuum = man->FindOrBuildMaterial("G4_Galactic");
  //G4Material* air =man->FindOrBuildMaterial("G4_AIR");

   fWorldMater = Vacuum;

  G4Material* NaI=man->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  fAuxiliaryDetectorMater=NaI;
 
  
  // or use G4 materials data base
  //
  // G4NistManager* man = G4NistManager::Instance();
  //  G4Element* Aluminium = man->FindOrBuildElement("Al");
  //  G4Material* Al1 = new G4Material("Aluminium",2.7*g/cm3,(G4int) 1);
  //  Al1->AddElement(Aluminium);
  //  fTargetMater = man->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  fTargetMater = man->FindOrBuildMaterial("G4_Al");
                   
 ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // World
  //
  // (re) compute World dimensions if necessary
  //fWorldLength = 2.*fTargetLength + fDetectorThickness;
  //fWorldRadius = fTargetRadius + fDetectorThickness;
    
  // G4Box*
  //sWorld = new G4Box("World",                                 //name
   //		     0.,fWorldRadius, fWorldLength, 0.,twopi); //dimensions    
  
  G4Box*
  sWorld = new G4Box("World",                                 //name
		     2.*fWorldRadius,2.*fWorldRadius, 2.*fWorldLength); //dimensions  
                   
  G4LogicalVolume*
  lWorld = new G4LogicalVolume(sWorld,                  //shape
                             fWorldMater,               //material
                             "World");                  //name

  fPhysiWorld = new G4PVPlacement(0,                    //no rotation
                            G4ThreeVector(),            //at (0,0,0)
                            lWorld,                     //logical volume
                            "World",                    //name
                            0,                          //mother volume
                            false,                      //no boolean operation
                            0);                         //copy number
                            
  // Target
  //
  //  G4Tubs* 
  //sTarget = new G4Tubs("Target",                                   //name
  //                0., fTargetRadius, 0.5*fTargetLength, 0.,twopi); //dimensions

  
  G4Box* 
  sTarget = new G4Box("Target",                                   //name
                    fTargetRadius, fTargetRadius, fTargetLength); //dimensions


  fLogicTarget = new G4LogicalVolume(sTarget,           //shape
                             fTargetMater,              //material
                             "Target");                 //name
                               
           new G4PVPlacement(0,                         //no rotation
			     G4ThreeVector(0, 0, (fDetectorThickness + 0.5*cm)),             //at (0,0,0)
                           fLogicTarget,                //logical volume
                           "Target",                    //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  // Detector
  //
  //G4Tubs* 
  //sDetector = new G4Tubs("Detector",  
  //            fTargetRadius, fWorldRadius, 0.5*fDetectorLength, 0.,twopi);
  
  G4Box* 
  sDetector = new G4Box("Detector",  
			fDetectorLength, fDetectorThickness, fDetectorThickness);


  fLogicDetector = new G4LogicalVolume(sDetector,       //shape
                             fDetectorMater,            //material
                             "Detector");               //name
                               
           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(),             //at (0,0,0)
                           fLogicDetector,              //logical volume
                           "Detector",                  //name
                           lWorld,                      //mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

//Sodium Iodide detectors
G4Tubs* sAuxiliaryDetector=new G4Tubs("sAuxiliaryDetector",
  0,fAuxiliaryDetectorRadius,0.5*fAuxiliaryDetectorLength,0.,twopi);

fLogicAuxiliaryDetector=new G4LogicalVolume(sAuxiliaryDetector,fAuxiliaryDetectorMater,"AuxiliaryDetector");

new G4PVPlacement(0,
  G4ThreeVector(0,0,fDetectorThickness + fAuxiliaryDetectorLength*0.5+fDetectorLength*0.5+0.5*cm),
  fLogicAuxiliaryDetector,
  "AuxiliaryDetector",
  lWorld,
  false,
  0);

//sodium iodide detector shield
//This object is made of two parts: the top and the part sourrounfing the detector. It is all in alluminium
G4Tubs* sAuxiliaryCover1=new G4Tubs("sAuxiliaryCover1",0,fAuxiliaryDetectorRadius+4*mm,1*mm*0.5,0,twopi);

G4Tubs* sAuxiliaryCover2=new G4Tubs("sAuxiliaryCover2",fAuxiliaryDetectorRadius+2*mm,fAuxiliaryDetectorRadius+4*mm,fAuxiliaryDetectorLength*0.5,0.,twopi);

G4UnionSolid* sAuxiliaryCover=new G4UnionSolid("sAuxiliaryCover",sAuxiliaryCover1,sAuxiliaryCover2,0,G4ThreeVector(0,0,fAuxiliaryDetectorLength*0.5));

G4LogicalVolume* lAuxiliaryCover=new G4LogicalVolume(sAuxiliaryCover,fTargetMater,"AuxiliaryCover");

new G4PVPlacement(0,
  G4ThreeVector(0,0,fDetectorThickness +fDetectorLength*0.5+0.4*cm),
  lAuxiliaryCover,
  "AuxiliaryCover",
  lWorld,
  false,
  0);

G4VisAttributes* visAuxiliaryCover = new G4VisAttributes();
  visAuxiliaryCover->SetVisibility(true);
  visAuxiliaryCover->SetForceWireframe(true);

//lAuxiliaryCover->SetVisAttributes(visAuxiliaryCover);

  PrintParameters();
  
  //always return the root volume
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n Target : Length = " << G4BestUnit(fTargetLength,"Length")
         << " Radius = " << G4BestUnit(fTargetRadius,"Length")  
         << " Material = " << fTargetMater->GetName();
  G4cout << "\n Detector : Length = " << G4BestUnit(fDetectorLength,"Length")
         << " Tickness = " << G4BestUnit(fDetectorThickness,"Length")  
         << " Material = " << fDetectorMater->GetName() << G4endl;          
  G4cout << "\n" << fTargetMater << "\n" << fDetectorMater << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fTargetMater = pttoMaterial;
    if(fLogicTarget) { fLogicTarget->SetMaterial(fTargetMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetTargetMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fDetectorMater = pttoMaterial;
    if(fLogicDetector) { fLogicDetector->SetMaterial(fDetectorMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetRadius(G4double value)
{
  fTargetRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetLength(G4double value)
{
  fTargetLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorThickness(G4double value)
{
  fDetectorThickness = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorLength(G4double value)
{
  fDetectorLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetLength()
{
  return fTargetLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetRadius()
{
  return fTargetRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetTargetMaterial()
{
  return fTargetMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicTarget()
{
  return fLogicTarget;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorLength()
{
  return fDetectorLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorThickness()
{
  return fDetectorThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetDetectorMaterial()
{
  return fDetectorMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicDetector()
{
  return fLogicDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
