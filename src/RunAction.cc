// Make this appear first!
#include "G4Timer.hh"

#include "RunAction.hh"
#include "RunActionMessenger.hh"

#include "G4Run.hh"

#include "Analysis.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "EventAction.hh"

#include <string>

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
 : G4UserRunAction(),fDetector(det),fPrimary(kin),
   fTimer(0)
{
  fTimer = new G4Timer;
  fMan = G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fTimer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  fMan = G4AnalysisManager::Instance();

  G4String targetName = fDetector->GetTargetMaterialName();
  G4String targetThickness = G4BestUnit(fDetector->GetTargetSize(),"Length");
  G4int sourceName = fPrimary->GetSourceType();
  G4String sourceString;

  // if(sourceName == 0 || sourceName==7){
  //   G4String sourceEnergy = G4BestUnit(fPrimary->GetSourceEnergy(),"Energy");
  //   sourceString = "gamma-"+ sourceEnergy;
  // }
  // else{
  //   sourceString = fPrimary->GetParticleName();
  // }

  //G4int directory = fDetector->GetDetectorType();
  //G4String directorName;

  //G4String position = G4BestUnit((fDetector->GetTargetSize()-fPrimary->GetSourcePosition()),"Length");

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  fTimer->Start();
  G4String abs_string = std::to_string(1);
  replace(abs_string.begin(),abs_string.end(),'.',',');

  G4String alpha_string = std::to_string(0.9);
  replace(alpha_string.begin(),alpha_string.end(),'.',',');

  fFileName = "../Data/"+fDetector->GetDetectorName()+".csv";

  fMan->SetVerboseLevel(0);
  fMan->OpenFile(fFileName);

  fMan->CreateNtuple("DetectedPhotons","Escapes");
  fMan->CreateNtupleDColumn("ABS");
  fMan->CreateNtupleDColumn("Alpha");
  fMan->CreateNtupleDColumn("Wavelength");

  fMan->FinishNtuple();
  fMan->FillNtupleDColumn(0,0,0.1);
  fMan->FillNtupleDColumn(0,1,0.95);
}

void RunAction::SetFileName(G4String fileName)
{
  fFileName = fileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  fTimer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent()
         << " " << *fTimer << G4endl;
    G4cout << "End Run" << G4endl;

    fMan->Write();
    fMan->CloseFile();
    delete fMan;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
