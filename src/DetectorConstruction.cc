#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "SiliconPlateConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Torus.hh"
#include "G4Hype.hh"

#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4PSEnergyDeposit.hh"
#include <G4VPrimitiveScorer.hh>

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VoxelLimits.hh"

#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

#include "G4GDMLParser.hh"

#include <G4VisAttributes.hh>
#include <iostream>
#include <fstream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
Constructs DetectorConstruction, defines default values.
*/
DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),fPBox(nullptr), fLBox(nullptr),
  fBox(nullptr)
{
  fDetectorMessenger = new DetectorMessenger(this);
  fTargetMPT = new G4MaterialPropertiesTable();
  fExpHall_x = fExpHall_y = fExpHall_z = 0.15*m;
  fTargetName = "holder";
  fThickness = 1*mm;
  fTargetThickness = 3*mm;
  fDetectorType = 4;
  fABSL = 1;
  fRES=4.0;
  fLY=10500./MeV;
  fDetectorName = "6pmt_coverage_pe";
  fVolName = "World";
  fSigAlpha = 0.5;
  DefineMaterials();
//  SetTargetMaterial("Scint");
  SetWorldMaterial("Air");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
Sets thickness of target.
*/
void DetectorConstruction::SetSize(G4double value){
  fTargetThickness=value;
  if(fBox){
    fBox->SetZHalfLength(fTargetThickness/2);
  }
  UpdateGeometry();

  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetPosition(G4double value){
  fPosition=value;
  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetLY(G4double value){
  fLY=value;
  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetRes(G4double value){
  fRES=value;
  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets which detector geometry is used.
*/
void DetectorConstruction::SetDetectorType(G4int value){
  fDetectorType=value;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetABS(G4double value){
  fABSL=value;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  DefineMaterials();
}

void DetectorConstruction::SetSigAlpha(G4double value){
  fSigAlpha=value;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  DefineMaterials();
}

void DetectorConstruction::SetDetectorName(G4String name){
  fDetectorName=name;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets material of target.
*/
void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fTargetMaterial = pttoMaterial;
    fTargetName = fTargetMaterial->GetName();
    if ( fLBox ) { fLBox->SetMaterial(fTargetMaterial); }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets material of world volume.
*/
void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fWorldMaterial = pttoMaterial;
    if ( fWLBox ) { fWLBox->SetMaterial(fWorldMaterial); }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Defines materials used in simulation. Sets material properties for PEN and other optical components.
*/
void DetectorConstruction::DefineMaterials(){// ------------- Materials -------------
  G4double a, z, density;
  G4int nelements;

  // fAir
  //
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  fAir = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  fAir->AddElement(N, 70.*perCent);
  fAir->AddElement(O, 30.*perCent);

  G4NistManager* man = G4NistManager::Instance();
  // Water
  //
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);

  G4Element* C = new G4Element("Carbon", "C", z=12, a=12*g/mole);
  G4Element* Pb = new G4Element("Lead", "Pb", z=87, a=207*g/mole);

  // Scintillators

  fPEN = new G4Material("PEN", density= 1.3*g/cm3, nelements=3);
  G4int number_of_atoms;
  fPEN->AddElement(O, number_of_atoms=4);
  fPEN->AddElement(H, number_of_atoms=10);
  fPEN->AddElement(C, number_of_atoms=14);

  G4double wavelength;
  char filler;
  G4double varabsorlength;
  G4double ems;
  G4double rindex;

  G4double absEnergy[102]  = {0};
  G4double abs[102]={0};
  G4double emission[102]={0};
  G4double rIndex[102]={0};
  G4double rIndex_fAir[102]={0};
  G4double ems_abs[102]={0};

  G4int absEntries = 0;
  ifstream ReadAbs;

  G4String abs_file = "../input_files/ABSTest.csv";
  G4double emission_fibre[102]={0};
  ReadAbs.open(abs_file);
  G4double var = GetABS();
  if(ReadAbs.is_open())
  {
    while(!ReadAbs.eof())
    {
      ReadAbs>>wavelength>>filler>>varabsorlength>>filler>>ems>>filler>>rindex;
      if(ReadAbs.eof()){
        break;
      }
      absEnergy[absEntries] = (1240/wavelength)*eV;
      abs[absEntries] = varabsorlength*mm;
// G4cout << varabsorlength << G4endl;
      emission[absEntries] = ems;
      rIndex[absEntries] = 1.65;
      rIndex_fAir[absEntries]=1.0;
      ems_abs[absEntries]=0.0;
      emission_fibre[absEntries]=1.0;
      absEntries++;
    }
  }

  else G4cout<<"Error opening file: " <<abs_file<<G4endl;
  ReadAbs.close();
  absEntries--;

  const G4int nEntries1 = sizeof(absEnergy)/sizeof(G4double);
  assert(sizeof(rIndex) == sizeof(absEnergy));
  assert(sizeof(abs) == sizeof(absEnergy));
  assert(sizeof(emission) == sizeof(absEnergy));
  assert(sizeof(rIndex_fAir == sizeof(absEnergy)));

  fTargetMPT->AddProperty("RINDEX",       absEnergy, rIndex, nEntries1)->SetSpline(true);
  fTargetMPT->AddProperty("ABSLENGTH",    absEnergy, abs, nEntries1)->SetSpline(true); // *
  fTargetMPT->AddProperty("FASTCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);
  fTargetMPT->AddProperty("SLOWCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);

  fTargetMPT->AddConstProperty("SCINTILLATIONYIELD",10500./MeV); // * 2.5 * PEN = PS, 10*PEN=PS
  fTargetMPT->AddConstProperty("RESOLUTIONSCALE",4.0); // * 1, 4, 8
  fTargetMPT->AddConstProperty("FASTTIMECONSTANT", 5.198*ns);
  fTargetMPT->AddConstProperty("SLOWTIMECONSTANT",24.336*ns);
  fTargetMPT->AddConstProperty("YIELDRATIO",0.05);

  fPEN->SetMaterialPropertiesTable(fTargetMPT);

  density = universe_mean_density;    //from PhysicalConstants.h
  fVacuum = new G4Material("Galactic", z=1., a=1.008*g/mole, density,
                           kStateGas,2.73*kelvin,3.e-18*pascal);
  //
  // fAir
  G4MaterialPropertiesTable* worldMPT = new G4MaterialPropertiesTable();
  worldMPT->AddProperty("RINDEX", absEnergy, rIndex_fAir, nEntries1)->SetSpline(true);

  fAir->SetMaterialPropertiesTable(worldMPT);
  fVacuum->SetMaterialPropertiesTable(worldMPT);
  fSi = man->FindOrBuildMaterial("G4_Si");
  fAl = man->FindOrBuildMaterial("G4_Al");
}

void DetectorConstruction::SetVolName(G4ThreeVector thePoint){
  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4VPhysicalVolume* myVolume= theNavigator->LocateGlobalPointAndSetup(thePoint);
  fVolName =  myVolume->GetName();
}

void DetectorConstruction::SetPropertyTable(G4Material* material, G4MaterialPropertiesTable* table){
  material->SetMaterialPropertiesTable(table);
}

void DetectorConstruction::UpdateGeometry(){
  DefineMaterials();
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

/*
Clears stored geometry, then constructs all volumes that can be used in the simulation.

Builds and places volumes in world.

Defines detector sensitivities and properties.
*/
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4GDMLParser parser;
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

// ------------- Volumes --------------

// The experimental Hall
  fWorldBox = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);

  fWLBox = new G4LogicalVolume(fWorldBox,fWorldMaterial,"World",0,0,0);

  fWPBox = new G4PVPlacement(0,G4ThreeVector(),fWLBox,"World",0,false,0);

  double target_width = 15*mm;
  fBox = new G4Box("target", target_width, target_width, fTargetThickness);
  fLBox = new G4LogicalVolume(fBox,fPEN, "target",0,0,0);
  double position = 0;

  double rod_sides = 2.5*mm;
  double rod_length = 50*mm;

  double curve = 2.5*mm;

  G4Box* rod_stick = new G4Box("target",rod_sides,rod_sides,rod_length);

  G4RotationMatrix* rm = new G4RotationMatrix(0,-90*deg,0);
  G4Tubs* curve_edge = new G4Tubs("curve",0,curve,rod_sides,0,180.*deg);

  G4UnionSolid* rod = new G4UnionSolid("final_rod",rod_stick, curve_edge, rm, G4ThreeVector(0,0,-rod_length));
  G4LogicalVolume* rod_log = new G4LogicalVolume(rod,fPEN,"target",0,0,0);

  G4Box* collimator = new G4Box("collimator",15*mm, 5*mm, 10*mm);


  G4Tubs* col_hole = new G4Tubs("hole", 0, 1*mm, 10*mm, 0, 360.*deg);
  G4SubtractionSolid* col_comp = new G4SubtractionSolid("final_col", collimator, col_hole,rm,G4ThreeVector(0,0,0));
  G4LogicalVolume* col_log = new G4LogicalVolume(col_comp, fAl,"collimator",0,0,0);
  G4NistManager* man = G4NistManager::Instance();

  G4ThreeVector point = G4ThreeVector(0,0,5*cm);
  G4Navigator* pointNavigator = new G4Navigator();
  pointNavigator->SetWorldVolume(fWPBox);
  pointNavigator->LocateGlobalPointAndSetup(point);

  // --------------Detectors--------------

  char filler;
  G4double wavelength;
  G4double cath_eff;
  G4double photocath_energy[57];
  G4double photocath_EFF[57];
  G4double perfect_EFF[57];
  G4double perfect_REFL[57];
  G4double photocath_REFL[57]={0};
  G4String pmt_file = "../input_files/pmtQE.csv";

  ifstream ReadEff;
  G4int effCounter = 0;
  ReadEff.open(pmt_file);

  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>wavelength>>filler>>cath_eff;
      if(ReadEff.eof()){
        break;
      }
      photocath_energy[57-effCounter] = (1240/wavelength)*eV;
      photocath_EFF[57-effCounter] = cath_eff;
      perfect_EFF[57-effCounter] = 1;
      perfect_REFL[57-effCounter] = 0;
      effCounter++;
    }
  }

  else G4cout<<"Error opening file: " <<pmt_file<<G4endl;
  ReadEff.close();
  effCounter--;

  const G4int nPMT_EFF = sizeof(photocath_energy)/sizeof(G4double);

  G4OpticalSurface* perfect_optsurf = new G4OpticalSurface("perfect",glisur,polished, dielectric_metal);
  G4MaterialPropertiesTable* detector_MT = new G4MaterialPropertiesTable();
  detector_MT->AddProperty("EFFICIENCY", photocath_energy, perfect_EFF,nPMT_EFF);
  detector_MT->AddProperty("REFLECTIVITY", photocath_energy, perfect_REFL,nPMT_EFF);
  perfect_optsurf->SetMaterialPropertiesTable(detector_MT);

  // Spectrometer Sensor

  G4double spectrometer_side = 2.5*mm;
  G4double spectrometer_top = 2.5*mm;
  G4double spectrometer_thickness = 0.01*mm;

  G4Box* spectrometer_sensor = new G4Box("spec_sensor",spectrometer_top,spectrometer_side,spectrometer_thickness);
  G4LogicalVolume* spec_log = new G4LogicalVolume(spectrometer_sensor,fSi,"spec_sensor_log");
  G4LogicalVolume* abs_log = new G4LogicalVolume(spectrometer_sensor,fSi,"spec_abs_log");

  G4double spec_energy[64];
  G4double spec_EFF[64];
  G4double spec_REFL[64]={0};
  G4String spec_file = "../input_files/specQE.csv";

  effCounter = 0;
  ReadEff.open(spec_file);

  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>wavelength>>filler>>cath_eff;
      if(ReadEff.eof()){
        break;
      }
      spec_energy[64-effCounter] = (1240/wavelength)*eV;
      spec_EFF[64-effCounter] = cath_eff;
      effCounter++;
    }
  }

  else G4cout<<"Error opening file: " <<spec_file<<G4endl;
  ReadEff.close();
  effCounter--;

  const G4int nSPEC_EFF = sizeof(spec_energy)/sizeof(G4double);

  G4OpticalSurface* spec_optsurf = new G4OpticalSurface("photocath_opsurf",glisur,polished, dielectric_metal);
  G4MaterialPropertiesTable* spec_MT = new G4MaterialPropertiesTable();
  spec_MT->AddProperty("EFFICIENCY", spec_energy, spec_EFF,nSPEC_EFF);
  spec_MT->AddProperty("REFLECTIVITY", spec_energy, spec_REFL,nSPEC_EFF);
  spec_optsurf->SetMaterialPropertiesTable(spec_MT);
  new G4LogicalSkinSurface("spec_surf",spec_log,spec_optsurf);
  //
  G4OpticalSurface* AirPEN = new G4OpticalSurface("AirPEN",glisur, ground, dielectric_dielectric);
  AirPEN -> SetPolish(0.9);
  AirPEN -> SetMaterialPropertiesTable(fTargetMPT);

  G4VPhysicalVolume* siPM_placement;
  G4VPhysicalVolume* siPM_case_placement;
  G4VPhysicalVolume* siPM_board_placement;
  G4VPhysicalVolume* col_placement;

  G4VPhysicalVolume* vacPlacement;

  G4VisAttributes* tileAttr = new G4VisAttributes(G4Colour::Blue());
  tileAttr->SetVisibility(true);
  rod_log->SetVisAttributes(tileAttr);


  /*
  0 - PMT on base of tile, collimator included.
  */
  //fDetectorType = 0;
  fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,rod_length),rod_log,"rod",fWLBox,false,0,false);
  siPM_placement = new G4PVPlacement(0, G4ThreeVector(0,0,-2.*cm),spec_log,"spec",fWLBox,false,0,false);
  siPM_board_placement = new G4PVPlacement(0, G4ThreeVector(0,0,2*rod_length+spectrometer_thickness),abs_log,"abs",fWLBox,false,0,false);
  col_placement = new G4PVPlacement(0, G4ThreeVector(0,7.5*mm,fPosition),col_log,"col",fWLBox,false,0,false);

  G4LogicalBorderSurface* surfaceAirPEN = new G4LogicalBorderSurface("AirPEN",fWPBox,fPBox,AirPEN);
  G4LogicalBorderSurface* surfacePENAir = new G4LogicalBorderSurface("AirPEN",fPBox,fWPBox,AirPEN);

  return fWPBox;
}
