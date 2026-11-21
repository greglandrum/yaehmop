# YAeHMOP: Yet Another extended Hueckel Molecular Orbital Package

This software was original written by Greg Landrum during his time in the Hoffmann Group at Cornell University. Modifications to support COHP analysis and f-orbitals were provided by Wingfield Glassey.

The license has been updated to the standard BSD 2-clause license for the github version of the code.


## Spin polarization and spin-orbit coupling

`bind` can do spin-polarized runs with one parameter table per spin channel
(`--parm-up`, `--parm-down`, optionally `--common-fermi`) and on-site spin-orbit
coupling (`--soc`), also combined. See [docs/SPIN_SOC.md](docs/SPIN_SOC.md); the
tests are in `tests/spin_soc/`. SOC needs the build with LAPACK
(`-DUSE_BLAS_LAPACK=ON`).
