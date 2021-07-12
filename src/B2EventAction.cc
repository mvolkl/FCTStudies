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
/// \file B2EventAction.cc
/// \brief Implementation of the B2EventAction class

#include "B2EventAction.hh"

#include "B2Analysis.hh"
#include "B2TrackerHit.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

// #include "G4Navigator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2EventAction::B2EventAction()
: G4UserEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2EventAction::~B2EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2EventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2EventAction::EndOfEventAction(const G4Event* event)
{
  // get number of stored trajectories

  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);

  // periodic printing
  G4int eventID = event->GetEventID();
  if ( eventID < 100 || eventID % 100 == 0) {
    G4cout << ">>> Event: " << eventID  << G4endl;
    if ( trajectoryContainer ) {
      G4cout << "    " << n_trajectories
             << " trajectories stored in this event." << G4endl;
    }    
    G4cout << "    "  
           << hc->GetSize() << " hits stored in this event" << G4endl;

    // hc->PrintAllHits();
  }

  auto eventManager = G4EventManager::GetEventManager();
  auto analysisManager = G4AnalysisManager::Instance();

  // if (hc->GetSize() > 2) eventManager->KeepTheCurrentEvent();

  // G4cout << G4endl;
  // G4cout << "Trajectories:"<<
  //   "--------------------------------------------------------------" << G4endl;
 
  //
  // loop over trajectories
  //
  const int idElectron = 11;
  const int idPositron = -11;

  std::map<G4int, bool> storeTrajectory; 
  // std::map<G4int, G4double> InitialMomentum;

  for(G4int i=0; i<n_trajectories; i++) {
  	G4int  trackID = (*trajectoryContainer)[i]->GetTrackID();

    double z = (*trajectoryContainer)[i]->GetPoint(0)->GetPosition().z() / cm;
    int pdg = (*trajectoryContainer)[i]->GetPDGEncoding();
    int parentID = (*trajectoryContainer)[i]->GetParentID();
  	double px = (*trajectoryContainer)[i]->GetInitialMomentum().x();
  	double py = (*trajectoryContainer)[i]->GetInitialMomentum().y();
  	double pz = (*trajectoryContainer)[i]->GetInitialMomentum().z();
  	double InitialMomentum = sqrt(px*px + py*py + pz*pz);


    storeTrajectory[trackID] = false;
    if ((pdg == idElectron || pdg == idPositron) && trackID == 1) {
    	  storeTrajectory[trackID] = true;
    	  analysisManager->FillNtupleIColumn(0, 0, eventID);
      	analysisManager->FillNtupleIColumn(0, 1, trackID);
      	analysisManager->FillNtupleDColumn(0, 2, InitialMomentum);
      	analysisManager->AddNtupleRow(0);
      	eventManager->KeepTheCurrentEvent();
        G4cout << "KR: trajectory stored: z position (in cm): " << z << ", pdg = " << pdg << G4endl;
    } 
    // G4cout << "KR: z position (in cm): " << z << ", pdg = " << pdg << G4endl;

  }
  for (size_t i=0; i<hc->GetSize(); ++i) {

  	B2TrackerHit* th = (B2TrackerHit*) hc->GetHit(i);

  	if (/*storeTrajectory[th->GetTrackID()] &&*/ th->GetTrackID()==1) {

      analysisManager->FillNtupleIColumn(1, 0, eventID);
      analysisManager->FillNtupleIColumn(1, 1, th->GetTrackID());
      analysisManager->FillNtupleIColumn(1, 2, th->GetChamberNb());
      analysisManager->FillNtupleDColumn(1, 3, th->GetPos().x()/m);
      analysisManager->FillNtupleDColumn(1, 4, th->GetPos().y()/m);
      analysisManager->FillNtupleDColumn(1, 5, th->GetPos().z()/m);
      // analysisManager->FillNtupleDColumn(0, 6, InitialMomentum[th->GetTrackID()]);
      analysisManager->FillNtupleDColumn(1, 6, th->GetMomentum());

      analysisManager->AddNtupleRow(1);

  	}

    
    
    
  }

  

}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
