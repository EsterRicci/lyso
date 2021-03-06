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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

//#include "Run.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Hit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
:G4UserEventAction()
{ } 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  //fEdep1 = fEdep2 = fWeight1 = fWeight2 = 0.;
  //fTime0 = -1*s;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
 G4VHitsCollection* hc_LYSO = event->GetHCofThisEvent()->GetHC(0);
G4VHitsCollection* hc_NaI = event->GetHCofThisEvent()->GetHC(1);

G4double Edep_LYSO=0;
G4double Edep_NaI[2];

for(int iDet=0;iDet<2;++iDet)
  Edep_NaI[iDet]=0;

for(int iH=0;iH<hc_LYSO->GetSize();++iH){
    HitClass* hit=static_cast<HitClass*>(hc_LYSO->GetHit(iH));
    Edep_LYSO+=hit->GetEdep(); //adding the energies of the steps inside each detector, identified with chamber number
}

for(int iH=0;iH<hc_NaI->GetSize();++iH){
    HitClass* hit=static_cast<HitClass*>(hc_NaI->GetHit(iH));
    Edep_NaI[hit->GetReplicaNb()]+=hit->GetEdep(); //adding the energies of the steps inside each detector, identified with chamber number
  }

  analysisManager->FillNtupleDColumn(0,Edep_LYSO);
  analysisManager->FillNtupleDColumn(1,Edep_NaI[0]);
  //analysisManager->FillNtupleDColumn(2,Edep_NaI[1]);

  analysisManager->AddNtupleRow(); 

  //G4cout<<"Edep LYSO "<<Edep_LYSO<<" Edep NaI 1 "<<Edep_NaI[0]<<" Edep NaI 2 "<<Edep_NaI[1]<<G4endl;
      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


