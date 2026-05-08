FC=gfortran
# FFLAGS=-O3 -cpp -fopenmp
FFLAGS=-O0 -g -cpp -fopenmp -fbacktrace -fcheck=all -Wall -Wextra -ffpe-trap=invalid,zero,overflow
USE_MPI ?= 0
MPI_FC ?= mpif90

ifeq ($(USE_MPI),1)
FC := $(MPI_FC)
FFLAGS += -DUSE_MPI
endif

SRC:=constants_mod.f90
SRC+=input_mod.f90
SRC+=rs_vec_base_mod.f90
SRC+=rs_vec_dimensions_mod.f90
SRC+=design_base_mod.f90
SRC+=design_1D_mod.f90
SRC+=design_2D_mod.f90
SRC+=design_3D_mod.f90
SRC+=parallel_subs_mod.f90
SRC+=medium_mod.f90
SRC+=operator_mod.f90
SRC+=allocator_multidim_mod.f90
SRC+=output_mod.f90
SRC+=target_source_mod.f90
SRC+=rs_operations_subs_mod.f90
SRC+=bicgstab_mod.f90
SRC+=optimization_problem_mod.f90

MAIN=ch_inv_design.f90
OBJ=$(SRC:.f90=.o) $(MAIN:.f90=.o)
EXC=OMxTRA_inv_design.e

%.o: %.f90 #wildcard rule, creation of *.o depends on *.f90
	$(FC) $(FFLAGS) -o $@ -c $<

$(EXC): $(OBJ)
	$(FC) -o $@ $^ ${FFLAGS}

.PHONY: clean
clean:
	rm -f *.o *.mod $(EXC)

.PHONY: mpi
mpi:
	$(MAKE) USE_MPI=1