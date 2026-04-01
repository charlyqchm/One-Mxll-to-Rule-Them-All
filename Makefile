FC ?= gfortran
USE_MPI ?= 0
MPI_FC ?= mpif90
FFLAGS=-O3 -cpp -fopenmp 
DBG_FFLAGS=-O0 -g -cpp -fopenmp -fbacktrace -fcheck=all -Wall -Wextra -ffpe-trap=invalid,zero,overflow

ifeq ($(USE_MPI),1)
FC := $(MPI_FC)
FFLAGS += -DUSE_MPI
endif
LIBS=-ldftbplus -lblas -llapack
INCLUDE_DIR=$(DFTB_INSTALL)/include/dftbplus/modfiles
LIB_DIR=$(DFTB_INSTALL)/lib

SRC:=constants_mod.f90
SRC+=input_mod.f90
SRC+=detector_mod.f90
SRC+=source_mod.f90
SRC+=classical_medium_mod.f90
SRC+=mxll_base_mod.f90
SRC+=mxll_1D_mod.f90
SRC+=mxll_2D_mod.f90
SRC+=mxll_3D_mod.f90
SRC+=q_sys_base_mod.f90
SRC+=q_sys_dftb_mod.f90
SRC+=factory_mod.f90
SRC+=q_group_mod.f90
SRC+=parallel_subs_mod.f90
SRC+=write_fields_subs_mod.f90
SRC+=output_mod.f90
SRC+=interactions_mod.f90
MAIN=mxim_mxll.f90
MOD=${SRC:.f90=.mod}
OBJ=${SRC:.f90=.o}
OBJ+=mxim_mxll.o
EXC=OMxRTA.e

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< -I$(INCLUDE_DIR)


$(EXC): $(OBJ)
			$(FC) -o $@ $^ ${FFLAGS} -L$(LIB_DIR) $(LIBS) 

.PHONY: mpi
mpi:
	$(MAKE) USE_MPI=1


.PHONY: mpi_debug
mpi_debug:
	$(MAKE) USE_MPI=1 FFLAGS='$(DBG_FFLAGS)'


.PHONY: clean
clean:
	      rm $(OBJ) *.mod $(EXC)