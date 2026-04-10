# Project Summary

## Scientific Objective

This project started from a narrow question:


It then expanded into a second question:

> Can those learned local signals be reused as chemistry-aware descriptors and mechanistic screens for formulation materials?

The project is therefore split into two layers:

1. `Degradation Mapper`: a reactive polymer model and AFIR-pathway analysis workflow.
2. `Preformulation screen`: a formulation-oriented reuse layer built on top of the trained model.

## Dataset and Modeling Setup

### Reactive polymer training data

- Source dataset: `OPoly26`
- Reactive slice extracted: `666,220` structures
- Leakage-safe train pool: `655,969`
- Validation: `5,151`
- Frozen test: `5,100`

The split was built to avoid AFIR-family leakage between train and test.

### Model line

- Zero-shot baseline: `UMA-s-1p1`
- Fine-tuned checkpoints:
  - reactive pilot `v1`
  - hard-slice top-up `v2`
  - reactive scale `100k`

### Formulation reuse layer

- Descriptor families completed:
  - `PEG/PEO`
  - `Gelatin`
  - `Pectin`
  - `HPC`
  - `Xanthan`
- SSE formulation table merged from:
  - `SSE_ML_V3.xlsx`

## Main Experimental Findings

### 1. Reactive fine-tuning is a large win

The strongest finished checkpoint was the `100k` reactive model:

- energy MAE: `0.2889 eV`
- force MAE: `0.01999 eV/A`

Compared with zero-shot:

- energy improved by `51.8%`
- force improved by `19.5%`

This is the clearest core result in the project.

### 2. Gains are strongest on hard reactive chemistry

The subgroup audit showed strong improvements on slices that matter for degradation-like chemistry:

- `ion_inserted`
- neutral charge
- larger atom-count systems

This supports the claim that the model is learning reactive polymer behavior rather than only easy near-equilibrium structures.

### 3. Pathway-level demos are viable

Two AFIR time-lapse demos were built:

- peptoid degradation pathway
- fluoropolymer degradation pathway

These are not just static structures. They show ordered frame sequences, hotspot bonds, bond stretching, and energy progression along the path.

### 4. Descriptor generation works

A standardized bond-stretch workflow generated polymer-specific descriptors:

- barrier proxy
- hotspot bond type
- stretch tolerance
- peak force

The descriptor profiles are distinct across the five polymer families in the repository.

### 5. Descriptor utility for SSE prediction is not established yet

The first leave-one-family-out SSE printability pilot produced a negative result:

- `Baseline A` balanced accuracy mean: `0.8448`
- `Baseline B` balanced accuracy mean: `0.8448`
- `Augmented` balanced accuracy mean: `0.5982`

So the current descriptor layer is chemically interpretable but not yet predictive enough for that downstream task.

### 6. Relaxation materially changes formulation ranking

The relaxed GPU screen was an important correction to the original static proxy screen.

Neutral ranking after relaxation:

1. `paracetamol + pvp + water`
2. `paracetamol + hpc + ethanol/water`
3. `paracetamol + peg + water`
4. `paracetamol + pectin + glycerol/water`

Because the order changed after relaxation, the project now treats the relaxed screen as the minimum trustworthy version.

## Practical Meaning

What this project can support now:

- ranking local polymer degradation risk
- generating chemically interpretable polymer descriptors
- screening local API-polymer-solvent clusters
- prioritizing a small number of DFT or wet-lab follow-ups

What it does not support yet:

- direct prediction of final formulation success
- direct prediction of printability
- direct comparison of ionic and neutral cases in one score table
- broad claims about bulk rheology or long-term stability

## Most Important Remaining Work

1. Add real excipient chemistry:
   - `HPMC`
   - `Carbopol`
2. Replace the paracetamol demo matrix with the actual target API systems.
3. Run DFT spot checks on top, middle, and weak relaxed candidates.
4. Validate whether the mechanistic ranking agrees with experiment.

## Files To Start From

- [`../results/reactive_model_performance.csv`](../results/reactive_model_performance.csv)
- [`../results/formulation_descriptor_family_summary.csv`](../results/formulation_descriptor_family_summary.csv)
- [`../results/mechanistic_screen_relaxed_neutral.csv`](../results/mechanistic_screen_relaxed_neutral.csv)
- [`../results/mechanistic_screen_relaxed_ionic.csv`](../results/mechanistic_screen_relaxed_ionic.csv)
- [`../results/sse_lofo_summary.csv`](../results/sse_lofo_summary.csv)
