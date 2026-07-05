#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

class ActionInitialization : public G4VUserActionInitialization
{
public:
    ActionInitialization(DetectorConstruction* detector);
    virtual ~ActionInitialization();

    virtual void Build() const;

private:
    DetectorConstruction* fDetector;
};

#endif
