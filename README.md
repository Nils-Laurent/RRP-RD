# RRP-RD
RRP-RD is a Ridge Detector (RD) based on Relevant Ridge Portions (RRP) corresponding to [1]. `RRP_alg/RRP_RD.m` is the matlab implementation of this ridge detector.

[1] N. Laurent and S. Meignen, "A Novel Ridge Detector for Nonstationary Multicomponent Signals: Development and Application to Robust Mode Retrieval," in IEEE Transactions on Signal Processing, vol. 69, pp. 3325-3336, 2021, doi: 10.1109/TSP.2021.3085113.

## Installation

Clone this repository using option a or b:
 a. `git clone --recurse-submodules git@github.com:Nils-Laurent/RRP-RD.git`
 b. `git clone git@github.com:Nils-Laurent/RRP-RD.git`
    and `git submodule init`
    and `git submodule update`

## Project structure
- `RRP_alg` contains the implementation of RRP-RD
- `test` contains necessary functions to run tests
- genfig matlab files generate the same figures as in the paper
