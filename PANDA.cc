#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
 #include "QGSP_BIC_HP.hh"
//#include "QGSP_INCLXX_HP.hh"
#include "G4GenericBiasingPhysics.hh"
#include "SEEBiasingOperator.hh"

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

#include <vector>

int main(int argc, char** argv)
{
    // Construct the default run manager
    auto* runManager =
        G4RunManagerFactory::CreateRunManager(
            G4RunManagerType::SerialOnly
        );

    // Detector geometry
    auto* detector = new DetectorConstruction();
    runManager->SetUserInitialization(detector);

    // Physics list
    auto* physicsList = new QGSP_BIC_HP();
        //new QGSP_INCLXX_HP()

    // Wrap physics for biasing. Following the official Geant4 biasing
    // examples (GB01/GB07) exactly: explicitly name the ONE process to
    // wrap -- proton inelastic scattering -- rather than wrapping
    // every physics process for the proton. This guarantees the
    // ionization process (responsible for essentially all "normal",
    // non-biased energy deposition) is never touched by the biasing
    // framework at all.
    //
    // This does NOT bias anything by itself -- whether biasing
    // actually does anything, and by how much, is entirely controlled
    // by SEEBiasingOperator, created and attached to the relevant
    // logical volumes in DetectorConstruction::ConstructSDandField(),
    // using the cross-section factor set via
    // /sim/biasCrossSectionFactor in run.mac (default 1.0 = no bias).
    auto* biasingPhysics = new G4GenericBiasingPhysics();
    std::vector<G4String> protonProcessesToBias;
    protonProcessesToBias.push_back("protonInelastic");
    biasingPhysics->PhysicsBias("proton", protonProcessesToBias);
    physicsList->RegisterPhysics(biasingPhysics);

    runManager->SetUserInitialization(physicsList);

    // User actions
    runManager->SetUserInitialization(
        new ActionInitialization(detector)
    );

    // Visualization
    auto* visManager = new G4VisExecutive();
    visManager->Initialize();

    // UI manager
    auto* UImanager =
        G4UImanager::GetUIpointer();

    if(argc == 1)
    {
        // Interactive mode
        auto* ui = new G4UIExecutive(argc, argv);
        UImanager->ApplyCommand("/control/execute run.mac");
        ui->SessionStart();
        delete ui;
    }
    else
    {
        // Batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }

    // Print ground-truth biasing interaction counts now, immediately
    // after the run finishes. Deliberately NOT relying on
    // ~SEEBiasingOperator() for this: that destructor may never run
    // (DetectorConstruction::ConstructSDandField() creates each
    // thread's operator with `new` and never deletes it, a common
    // Geant4 biasing example pattern) -- see SEEBiasingOperator.hh for
    // the full explanation. This call site is guaranteed to execute.
    SEEBiasingOperator::PrintTotals();

    // Cleanup
    delete visManager;
    delete runManager;

    return 0;
}
