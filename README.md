# Degradation Mapper

Degradation Mapper is a reactive polymer MLIP project built on `OPoly26` and `FAIRChem/UMA` to model local polymer failure chemistry and to test whether those learned local signals can be reused for formulation screening.

## Project Aim

The project has two linked aims:

1. Train a leakage-safe reactive polymer model that improves on a zero-shot molecular/polymer foundation model for degradation-like structures.
2. Reuse the trained model as a mechanistic screening layer for formulation questions such as API-polymer compatibility, solvent competition, polymer cohesion, ion sensitivity, and local reactive risk.

## What Was Done

### Reactive polymer modeling

- Extracted a reactive `OPoly26` slice with `666,220` structures.
- Built leakage-safe AFIR-family splits:
  - train pool: `655,969`
  - validation: `5,151`
  - frozen test: `5,100`
- Benchmarked:
  - zero-shot `UMA-s-1p1`
  - reactive pilot `v1`
  - hard-slice top-up `v2`
  - `100k` scaling run

### Pathway showcases

- Built a peptoid AFIR time-lapse demo.
- Built a fluoropolymer AFIR time-lapse demo.

### Formulation-oriented extension

- Built mapper-derived polymer descriptor scans for:
  - `PEG/PEO`
  - `Gelatin`
  - `Pectin`
  - `HPC`
  - `Xanthan`
- Merged descriptors into `SSE_ML_V3`.
- Ran a first leave-one-family-out SSE printability pilot.
- Built a mechanistic preformulation screen and ran a relaxed GPU demo for neutral and ionic candidate clusters.

## Main Results

### Reactive model performance on the frozen 5,100-structure test set

| model | train rows | epochs | energy MAE (eV) | force MAE (eV/A) | energy gain vs zero-shot |
| --- | ---: | ---: | ---: | ---: | ---: |
| UMA-s-1p1 zero-shot | 0 | 0 | 0.5991 | 0.02484 | 0.0% |
| Reactive pilot v1 | 40,308 | 3 | 0.3077 | 0.02023 | 48.6% |
| Reactive hard-slice v2 | 18,654 | 2 | 0.2949 | 0.02036 | 50.8% |
| Reactive scale 100k | 100,375 | 3 | 0.2889 | 0.01999 | 51.8% |

The finished `100k` checkpoint is the strongest completed model so far.

### Hard reactive slices improved

- `ion_inserted` energy improved from `1.1255 eV` to `0.4895 eV` (`56.5%`)
- neutral `charge = 0` energy improved from `0.8753 eV` to `0.4513 eV` (`48.4%`)
- `226-250` atom systems improved from `1.4250 eV` to `0.7573 eV` (`46.9%`)

### AFIR showcase summaries

- Peptoid pathway:
  - `11` frames
  - `121` atoms
  - hotspot bond `C-C`
  - bond stretch `1.538 A -> 4.844 A`
- Fluoropolymer pathway:
  - `19` frames
  - `132` atoms
  - hotspot bond `C-C`
  - bond stretch `1.427 A -> 5.223 A`

### Formulation descriptor pilot

The descriptor engine produced interpretable bond-fragility descriptors for five formulation-relevant polymer families:

| polymer | hotspot bond | barrier proxy (eV) | stretch tolerance (A) | peak force (eV/A) |
| --- | --- | ---: | ---: | ---: |
| PEG/PEO | O-C | 4.3301 | 1.7278 | 3.1942 |
| Gelatin | C-N | 4.0654 | 3.3429 | 3.3869 |
| Pectin | O-C | 4.5937 | 3.4284 | 3.1147 |
| HPC | O-C | 4.4491 | 1.7356 | 2.6085 |
| Xanthan | O-C | 4.4279 | 3.4381 | 2.9208 |

### SSE leave-one-family-out pilot

This first formulation ML pilot was useful as a negative result:

- `Baseline A` and `Baseline B` both reached mean balanced accuracy `0.8448`
- the descriptor-augmented model dropped to `0.5982` on the core mixed-family subset

So the current descriptor set is chemically meaningful, but it does **not** yet add predictive value for SSE LOFO classification.

### Relaxed mechanistic screening demo

For the relaxed neutral paracetamol demo, the ranking was:

1. `PVP + water`
2. `HPC + ethanol/water`
3. `PEG + water`
4. `Pectin + glycerol/water`

The top neutral candidate stayed `PVP + water` after relaxation, but the middle ordering changed, which shows that relaxation is necessary before interpreting the scores.

## Current Findings

- Reactive fine-tuning on `OPoly26` clearly works.
- Leakage-safe splits were necessary and were enforced at the AFIR-family level.
- Scaling from pilot to `100k` improved the best finished model further.
- The model is strongest on the reactive chemistry it was fine-tuned for.
- Polymer descriptor generation is now real and reproducible.
- Descriptor usefulness for downstream formulation ML is not proven yet.
- Relaxed local compatibility screening is more trustworthy than static proxy ranking.

## Current Limits

- The current formulation screen is still a local-cluster proxy, not a free-energy workflow.
- Ionic screening is separated but not yet cross-comparable to neutral screening because monoatomic ion references are still approximate.
- `HPMC` and `Carbopol` proxies are not yet in the mechanistic screen.
- No experimental validation is in this repository yet.

## Repository Layout

- [`docs/project_summary.md`](docs/project_summary.md): fuller project narrative
- [`results/`](results): compact CSV tables for the main outcomes
- [`figures/`](figures): selected figures for scaling, hard slices, and pathway demos

## Next Steps

1. Add `HPMC` and `Carbopol` proxies to the mechanistic screen.
2. Replace the demo matrix with the real target API and excipient systems.
3. Run DFT spot checks on top, middle, and weak relaxed candidates.
4. Use the validated screen to prioritize experimental formulation work.
