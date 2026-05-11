FC ?= gfortran
OPT_FFLAGS ?= -O3 -cpp -fopenmp
DBG_FFLAGS ?= -O0 -g -cpp -fopenmp -fbacktrace -fcheck=all -Wall -Wextra -ffpe-trap=invalid,zero,overflow
FFLAGS ?= $(OPT_FFLAGS)
USE_MPI ?= 0
MPI_FC ?= mpif90
DEBUG ?= 0

ifeq ($(DEBUG),1)
FFLAGS := $(DBG_FFLAGS)
endif

ifeq ($(USE_MPI),1)
FC := $(MPI_FC)
FFLAGS += -DUSE_MPI
endif

BUILD ?= ../../build
PROJECT ?= design
OBJ_DIR ?= $(BUILD)/$(PROJECT)/obj
MOD_DIR ?= $(BUILD)/$(PROJECT)/mod
BIN_DIR ?= $(BUILD)/bin

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
OBJ=$(addprefix $(OBJ_DIR)/,$(SRC:.f90=.o))
OBJ+=$(OBJ_DIR)/$(MAIN:.f90=.o)
EXC=$(BIN_DIR)/OMxRTA_inv_design.e

.PHONY: all clean mpi debug mpi_debug

all: $(EXC)

$(OBJ_DIR) $(MOD_DIR) $(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: %.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -o $@ -c $< -I$(MOD_DIR) -J$(MOD_DIR)


$(EXC): $(OBJ) | $(BIN_DIR)
	$(FC) -o $@ $^ ${FFLAGS}

clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR) $(EXC)

mpi:
	$(MAKE) USE_MPI=1

debug:
	$(MAKE) DEBUG=1

mpi_debug:
	$(MAKE) USE_MPI=1 DEBUG=1