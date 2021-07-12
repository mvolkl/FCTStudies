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
/// \file B2aDetectorConstruction.cc
/// \brief Implementation of the B2aDetectorConstruction class
 
#include "B2aDetectorConstruction.hh"
#include "B2aDetectorMessenger.hh"
#include "B2TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4ThreadLocal 
G4GlobalMagFieldMessenger* B2aDetectorConstruction::fMagFieldMessenger = 0;

B2aDetectorConstruction::B2aDetectorConstruction()
:G4VUserDetectorConstruction(), 
 fNbOfChambers(0),
 fLogicTarget(NULL), fLogicChamber(NULL), 
 fTargetMaterial(NULL), fChamberMaterial(NULL), 
 fStepLimit(NULL),
 fCheckOverlaps(true)
{
  fMessenger = new B2aDetectorMessenger(this);

  // Geometry setting A
  // fNbOfChambers = 10;

  // Geometry setting B
  // fNbOfChambers = 5;

  // Geomtry setting C
  fNbOfChambers = 7;
  
  fLogicChamber = new G4LogicalVolume*[fNbOfChambers];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
B2aDetectorConstruction::~B2aDetectorConstruction()
{
  delete [] fLogicChamber; 
  delete fStepLimit;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* B2aDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::DefineMaterials()
{
  // Material definition 

  G4NistManager* nistManager = G4NistManager::Instance();

  // Air defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  // Al defined using NIST Manager
  // fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_Al");

  // Lead defined using NIST Manager
  fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_Pb");

  // Xenon gas defined using NIST Manager
  // fChamberMaterial = nistManager->FindOrBuildMaterial("G4_Xe");

  // Silicon defined using NIST Manager
  fChamberMaterial = nistManager->FindOrBuildMaterial("G4_Si");

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B2aDetectorConstruction::DefineVolumes()
{
  G4Material* air  = G4Material::GetMaterial("G4_AIR");

  // Sizes of the principal geometrical components (solids)
  
  //
  // Geometry A (10 chambers, definend above 2 cm and 10 cm spacing)
  //
  // G4double chamberSpacingSector1 = 2*cm; // from chamber center to center!
  // G4double chamberSpacingSector2 = 10*cm; // from chamber center to center!
  // G4int fNbOfChambersSector1 = 5;
  // G4int fNbOfChambersSector2 = fNbOfChambers - fNbOfChambersSector1;  

  // geometry A2 
  // G4double chamberSpacingSector1 = 0.5*cm; // from chamber center to center!
  // G4double chamberSpacingSector2 = 5*cm; // from chamber center to center!
  // G4int fNbOfChambersSector1 = 8;
  // G4int fNbOfChambersSector2 = fNbOfChambers - fNbOfChambersSector1;  

  //
  // Geometry B (5 chambers)
  //

  // geometry B1
  // G4double chamberSpacingSector1 = 1*cm; // from chamber center to center!
  // G4double chamberSpacingSector2 = 5*cm; // from chamber center to center!
  // G4int fNbOfChambersSector1 = 3;
  // G4int fNbOfChambersSector2 = fNbOfChambers - fNbOfChambersSector1;  

  // geometry B2
  // G4double chamberSpacingSector1 = 0.5*cm; // from chamber center to center!
  // G4double chamberSpacingSector2 = 3*cm; // from chamber center to center!
  // G4int fNbOfChambersSector1 = 3;
  // G4int fNbOfChambersSector2 = fNbOfChambers - fNbOfChambersSector1;  

  // geometry B3
  // G4double chamberSpacingSector1 = 1*cm; // from chamber center to center!
  // G4double chamberSpacingSector2 = 2*cm; // from chamber center to center!
  // G4int fNbOfChambersSector1 = 3;
  // G4int fNbOfChambersSector2 = fNbOfChambers - fNbOfChambersSector1;  


  //
  // Geometry C (11 chambers, positions of the chambers given explicitely)
  //
  G4double fChamberZPosRel[7] = {0.77*m, 1.0*m, 1.22*m, 1.5*m, 1.8*m, 2.2*m, 2.79*m};
  
  // G4double trackerLength = (fNbOfChambers+1)*chamberSpacing;

  // Geomtry A, B
  //G4double trackerLength = fNbOfChambersSector1 * chamberSpacingSector1 + 
  // 	(fNbOfChambersSector2+1) * chamberSpacingSector2;

  // Geomtry C
  G4double trackerLength = 3.2*m;//650*mm;
  
  //G4double chamberWidth = 0.05*mm; // width of the chambers
  G4double chamberWidth = 0.937*mm; // width of the chambers
  // G4double chamberHeightHalf = 150*mm;
  G4double chamberHeightHalf = 0.5*m;
  // G4double targetLength =  0.1*mm; // full length of Target, Geomtry C1
  G4double targetLength =  0.01*mm; // Geomtry C2
  G4double targetHeightHalf = chamberHeightHalf*0.9;
  G4double targetWidthHalf = chamberHeightHalf*0.9;

  // G4double worldLength = 1.2 * (2*targetLength + trackerLength);
  G4double worldHeightHalf = 1.2 * chamberHeightHalf;
  G4double worldLength = 1.2 * (2*targetLength + trackerLength);

  // G4double targetRadius  = 0.5*targetLength;   // Radius of Target
  targetLength = 0.5*targetLength*0.01;             // Half length of the Target  
  G4double trackerSize   = 0.5*trackerLength;  // Half length of the Tracker

  // Definitions of Solids, Logical Volumes, Physical Volumes

  // World

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLength);

  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;

  // G4Box* worldS
  //  = new G4Box("world",                                    //its name
  //              worldLength/2,worldLength/2,worldLength/2); //its size
  G4Box* worldS
    = new G4Box("world",                                    //its name
                worldHeightHalf,worldHeightHalf,worldLength/2); //its size
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,   //its solid
                 air,      //its material
                 "World"); //its name
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(), // at (0,0,0)
                 worldLV,         // its logical volume
                 "World",         // its name
                 0,               // its mother  volume
                 false,           // no boolean operations
                 0,               // copy number
                 fCheckOverlaps); // checking overlaps 

  // Target
  
  G4ThreeVector positionTarget = G4ThreeVector(0,0,-(targetLength+trackerSize));


  G4double fTargetZmin = (-(targetLength+trackerSize) - targetLength) / cm;
  G4double fTargetZmax = (-(targetLength+trackerSize) + targetLength) / cm;

  G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! target z position (cm) = "
         << -(targetLength+trackerSize) / cm 
         << ", range [" << fTargetZmin << ", " << fTargetZmax << "]" << G4endl;

  // G4Tubs* targetS
  //  = new G4Tubs("target",0.,targetRadius,targetLength,0.*deg,360.*deg);

  G4Box* targetS
  	= new G4Box("target", targetWidthHalf, targetHeightHalf, targetLength);

  fLogicTarget
    //= new G4LogicalVolume(targetS, fTargetMaterial,"Target",0,0,0);
    = new G4LogicalVolume(targetS, air,"Target",0,0,0);
  /*new G4PVPlacement(0,               // no rotation
                    positionTarget,  // at (x,y,z)
                    fLogicTarget,    // its logical volume
                    "Target",        // its name
                    worldLV,         // its mother volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps */

  G4cout << "Target is " << 2*targetLength/cm << " cm of "
         << fTargetMaterial->GetName() << G4endl;

  // Tracker
 
  G4ThreeVector positionTracker = G4ThreeVector(0,0,0);

  // G4Tubs* trackerS
  //   = new G4Tubs("tracker",0,trackerSize,trackerSize, 0.*deg, 360.*deg);

  // G4Box* trackerS
  // 	= new G4Box("tracker", trackerSize, trackerSize, trackerSize);
  G4Box* trackerS
    = new G4Box("tracker", 1.1*chamberHeightHalf, 1.1*chamberHeightHalf, trackerSize);

  G4LogicalVolume* trackerLV
    = new G4LogicalVolume(trackerS, air, "Tracker",0,0,0);  
  new G4PVPlacement(0,               // no rotation
                    positionTracker, // at (x,y,z)
                    trackerLV,       // its logical volume
                    "Tracker",       // its name
                    worldLV,         // its mother  volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps 

  // Visualization attributes

  G4VisAttributes* boxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));

  worldLV      ->SetVisAttributes(boxVisAtt);
  fLogicTarget ->SetVisAttributes(boxVisAtt);
  trackerLV    ->SetVisAttributes(boxVisAtt);

  // Tracker segments

  // G4cout << "There are " << fNbOfChambers << " chambers in the tracker region. "
  //        << G4endl
  //        << "The chambers are " << chamberWidth/cm << " cm of "
  //        << fChamberMaterial->GetName() << G4endl
  //        << "The distance between chamber is " << chamberSpacing/cm << " cm" 
  //        << G4endl;
  
  // G4double firstPosition = -trackerSize + chamberSpacing;
  // G4double firstLength   = trackerLength/10;
  // G4double lastLength    = trackerLength;

  G4double halfWidth = 0.5*chamberWidth;
  // G4double rmaxFirst = 0.5 * firstLength;

  // G4double rmaxIncr = 0.0;
  // if( fNbOfChambers > 0 ){
  //   rmaxIncr =  0.5 * (lastLength-firstLength)/(fNbOfChambers-1);
  //   if (chamberSpacing  < chamberWidth) {
  //      G4Exception("B2aDetectorConstruction::DefineVolumes()",
  //                  "InvalidSetup", FatalException,
  //                  "Width>Spacing");
  //   }
  // }

  for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {

      G4double Zposition = -trackerSize;
      
      //
      // Geomtry A, B
      //
      // if(copyNo < fNbOfChambersSector1) Zposition += (copyNo+1) * chamberSpacingSector1;
      // else Zposition += fNbOfChambersSector1 * chamberSpacingSector1 
      // + (copyNo + 1 - fNbOfChambersSector1) * chamberSpacingSector2;

      //
      // Geometry C
      //
      Zposition += fChamberZPosRel[copyNo];

      // G4double Zposition = firstPosition + copyNo * chamberSpacing;
      // G4double rmax =  rmaxFirst + copyNo * rmaxIncr;

      // G4Tubs* chamberS
      //   = new G4Tubs("Chamber_solid", 0, rmax, halfWidth, 0.*deg, 360.*deg);

      // G4Box* chamberS
  		// = new G4Box("Chamber_solid", 0.5*trackerSize, 0.5*trackerSize, halfWidth);
      G4Box* chamberS
        = new G4Box("Chamber_solid", chamberHeightHalf, chamberHeightHalf, halfWidth);
      
      fLogicChamber[copyNo] =
              new G4LogicalVolume(chamberS,fChamberMaterial,"Chamber_LV",0,0,0);

      fLogicChamber[copyNo]->SetVisAttributes(chamberVisAtt);

      new G4PVPlacement(0,                            // no rotation
                        G4ThreeVector(0,0,Zposition), // at (x,y,z)
                        fLogicChamber[copyNo],        // its logical volume
                        "Chamber_PV",                 // its name
                        trackerLV,                    // its mother  volume
                        false,                        // no boolean operations
                        copyNo,                       // copy number
                        fCheckOverlaps);              // checking overlaps 

  }

  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  G4double maxStep = 0.5*chamberWidth;
  fStepLimit = new G4UserLimits(maxStep);
  trackerLV->SetUserLimits(fStepLimit);
 
  /// Set additional contraints on the track, with G4UserSpecialCuts
  ///
  /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
  ///                                           maxLength,
  ///                                           maxTime,
  ///                                           minEkin));

  // Always return the physical world

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2aDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "B2/TrackerChamberSD";
  B2TrackerSD* aTrackerSD = new B2TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name 
  // of "Chamber_LV".
  SetSensitiveDetector("Chamber_LV", aTrackerSD, true);

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector(0.,0.,0.5*tesla);
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2aDetectorConstruction::SetTargetMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName);

  if (fTargetMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fTargetMaterial = pttoMaterial;
        if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
        G4cout 
          << G4endl 
          << "----> The target is made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetTargetMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetChamberMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fChamberMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fChamberMaterial = pttoMaterial;
        for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {
            if (fLogicChamber[copyNo]) fLogicChamber[copyNo]->
                                               SetMaterial(fChamberMaterial);
        }
        G4cout 
          << G4endl 
          << "----> The chambers are made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetChamberMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}  
