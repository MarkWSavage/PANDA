# PANDA

**PANDA** — *Protons And Neutron charge Deposition in mAterials* — is a Geant4 Monte Carlo simulation for calculating charge deposition from energetic particles (protons, neutrons, ions) in semiconductor structures, for Single Event Effect (SEE) analysis.

The core physical quantity computed is **collected charge** per event; upset probability, cross-section spectra, and other SEE metrics are derived from it in post-processing (`PANDA_Analyze.py`), keeping particle transport and interpretation cleanly separated.

Validated against CREME-MC. Includes optional cross-section biasing (`SEEBiasingOperator`) to efficiently sample the rare nuclear-recoil tail, with correct raw-charge/event-weight separation so biased and unbiased runs produce statistically consistent spectra.

## Build

```
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

## Run

```
./PANDA run.mac
python3 PANDA_Analyze.py
```

See `Documentation/PANDA_MASTER_DESIGN` for the full design philosophy.
