# RRP-RD
RRP-RD is a Ridge Detector (RD) based on Relevant Ridge Portions (RRP) corresponding to [1]. `RRP_alg/RRP_RD.m` is the matlab implementation of this ridge detector.

[1] N. Laurent and S. Meignen, "A Novel Ridge Detector for Nonstationary Multicomponent Signals: Development and Application to Robust Mode Retrieval," in IEEE Transactions on Signal Processing, vol. 69, pp. 3325-3336, 2021, doi: 10.1109/TSP.2021.3085113.

## Installation

Clone this repository using one of the two options:
- option 1: `git clone --recurse-submodules git@github.com:Nils-Laurent/RRP-RD.git`
- option 2:
  1. `git clone git@github.com:Nils-Laurent/RRP-RD.git`
  1.  and `git submodule init`
  1.  and `git submodule update`

## Usage
A minimal example is provided in `minimal_example.m`

## Project structure
- `RRP_alg/` contains the implementation of RRP-RD
- `test/` contains necessary functions to run tests
- `genfig_R1_*.m` matlab files generate the same figures as in the paper
- `minimal_example.m` provides a minimal example
