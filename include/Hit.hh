#ifndef Hit_h
#define Hit_h 1

#include "G4ThreeVector.hh"
#include "globals.hh"

struct Hit
{
    G4ThreeVector pos;   // deposition position
    G4double edep;       // raw deposited energy
    G4double time;       // global time

    G4String particle;   // particle species
    G4String process;   // NEW
    G4int trackID;       // unique track
    G4int parentID;      // ancestry
    G4int stepNumber;    // step index
};

#endif
