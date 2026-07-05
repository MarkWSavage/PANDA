#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(
    DetectorConstruction* detector)
 : fDetector(detector)
{}

ActionInitialization::~ActionInitialization()
{}

void ActionInitialization::Build() const
{
    // Primary particle source
    SetUserAction(new PrimaryGeneratorAction());

    // Event container
    auto eventAction = new EventAction();
    SetUserAction(eventAction);

    // Step-level hit recorder
    SetUserAction(
        new SteppingAction(fDetector, eventAction)
    );
}
