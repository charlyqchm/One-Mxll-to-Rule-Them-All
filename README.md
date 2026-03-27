# One-Mxll-to-Rule-Them-All (OMxRTA)

OMxRTA propagates macroscopic Maxwell's equations using the Finite-Difference Time-Domain (FDTD) method, coupled to the evolution of quantum-mechanical systems embedded in the simulation grid under the local dipole approximation. This implementation is based on the work of Maxim Sukharev:

* Sukharev, Maxim, and Abraham Nitzan, "Numerical studies of the interaction of an atomic sample with the electromagnetic field in two dimensions". _Phys. Rev. A._, 84(4), 043802 (2011).
* Sukharev, Maxim, "Efficient parallel strategy for molecular plasmonics–A numerical tool for integrating Maxwell-Schrödinger equations in three dimensions". _J. Comput. Phys._, 477, 111920 (2023).
* Sukharev, Maxim, Joseph Subotnik, and Abraham Nitzan, "Dissociation slowdown by collective optical response under strong coupling conditions". _J. Chem. Phys._, 158(8), 084104 (2023).
* Sukharev, Maxim, Joseph E. Subotnik, and Abraham Nitzan, "Unveiling the dance of molecules: Rovibrational dynamics of molecules under intense illumination at complex plasmonic interfaces". _JCTC_, 21(5), 2165-2178 (2025).

The code is designed to let you implement quantum models from scratch, or to interface with external open-source libraries for representing the quantum subsystems.
Currently, the only supported external backend is the open-source code [DFTB+](https://github.com/dftbplus/dftbplus).

The current implementation can run on a single node (with optional OpenMP parallelism) and also has the option to run in parallel using MPI.

## Requirements

- [DFTB+](https://github.com/dftbplus/dftbplus)
- LAPACK
- GCC / GNU Fortran (gfortran) >= 12.5.0

Optional:
- Open MPI >= 4.1.0 (for MPI parallel runs)

## Compilation

OMxRTA currently links against [DFTB+](https://github.com/dftbplus/dftbplus). Build and install DFTB+ first (run the following commands in the DFTB+ source directory):

```bash
FC=gfortran CC=gcc cmake -DINSTANCE_SAFE_BUILD=True -B _build .
cd _build/
cmake --build . -- -j 4
cmake --install .
```
Once DFTB+ has been installed, set the `DFTB_INSTALL` environment variable to the installation prefix. For example:

```bash
export DFTB_INSTALL=./dftbplus/_build/_install
```
Then build OMxRTA from the cloned repository directory:
```bash
make
```
To build the MPI-enabled version:
```bash
make mpi
```

## Running the code

Compilation produces the executable `OMxRTA.e`. Run it from the directory where you want to execute the simulation.
That directory must contain an input file named `inp`; the available input variables are documented in `input_mod.f90`.
Documentation for additional inputs (classical media and quantum systems) will be added in future updates.

## Output

The variables `mxll_n_detector` and `mxll_dt_det_print` of the `inp` file control how many outputs will be printed and how often. These outputs correspond to the components of the electric and magnetic fields at a point, along a line, on a surface or in an entire volume, depending on the dimensionality of the simulation box. These details have to be included in a file named `detectors.in`, and the number of lines must be equal to `mxll_n_detector`. Each detector will be printed in one or many files in a directory called `detector_xxxxxxx`. The number of the directory is equal to the corresponding line of the detector declared in the file `detectors.in`. For more details about the content of this file, check the subroutine `init_outputs` in `outputs_mod.f90`, and `init_detector` in `detector_mod.f90`.
More detailed information will be provided in the future documentation.

## Publications

1. [Bustamante, Carlos M., Franco P. Bonafé, Maxim Sukharev, Michael Ruggenthaler, Abraham Nitzan, and Angel Rubio, "Molecular polariton dynamics in realistic cavities". _JCTC_ 21(19), 9823-9831 (2025).](https://pubs.acs.org/doi/full/10.1021/acs.jctc.5c01318)
2. [Sidler, Dominik, Carlos M. Bustamante, Franco P. Bonafé, Michael Ruggenthaler, Maxim Sukharev, and Angel Rubio, "Density-functional tight binding meets Maxwell: unraveling the mysteries of (strong) light–matter coupling efficiently". _Nanophotonics_, 14(27), 4941-4955 (2025).](https://www.degruyterbrill.com/document/doi/10.1515/nanoph-2025-0453/html)
3. [Bustamante, Carlos M., Franco P. Bonafé, Richard Richardson, Michael Ruggenthaler, Wenxiang Ying, Abraham Nitzan, Maxim Sukharev, and Angel Rubio, "Collective Rabi-driven vibrational activation in molecular polaritons". _arXiv_, preprint arXiv:2601.16299 (2026).](https://arxiv.org/abs/2601.16299)
