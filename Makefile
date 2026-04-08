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

SRC_DIR ?= src
BUILD ?= build
BUILD_DIR ?= $(BUILD)
OBJ_DIR ?= $(BUILD_DIR)/obj
MOD_DIR ?= $(BUILD_DIR)/mod
BIN_DIR ?= $(BUILD_DIR)/bin

SRC_FILES:=constants_mod.f90
SRC_FILES+=input_mod.f90
SRC_FILES+=detector_mod.f90
SRC_FILES+=source_mod.f90
SRC_FILES+=classical_medium_mod.f90
SRC_FILES+=mxll_base_mod.f90
SRC_FILES+=mxll_1D_mod.f90
SRC_FILES+=mxll_2D_mod.f90
SRC_FILES+=mxll_3D_mod.f90
SRC_FILES+=q_sys_base_mod.f90
SRC_FILES+=q_sys_dftb_mod.f90
SRC_FILES+=factory_mod.f90
SRC_FILES+=q_group_mod.f90
SRC_FILES+=parallel_subs_mod.f90
SRC_FILES+=write_fields_subs_mod.f90
SRC_FILES+=output_mod.f90
SRC_FILES+=interactions_mod.f90
MAIN_FILE:=mxim_mxll.f90

SRC := $(addprefix $(SRC_DIR)/,$(SRC_FILES))
MAIN := $(SRC_DIR)/$(MAIN_FILE)

OBJ := $(addprefix $(OBJ_DIR)/,$(SRC_FILES:.f90=.o))
OBJ += $(OBJ_DIR)/$(MAIN_FILE:.f90=.o)

EXE_NAME ?= OMxRTA.e
EXC := $(BIN_DIR)/$(EXE_NAME)

all: $(EXC)

$(OBJ_DIR) $(MOD_DIR) $(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -o $@ -c $< -I$(INCLUDE_DIR) -I$(MOD_DIR) -J$(MOD_DIR)


$(EXC): $(OBJ) | $(BIN_DIR)
	$(FC) -o $@ $^ $(FFLAGS) -L$(LIB_DIR) $(LIBS)

.PHONY: mpi
mpi:
	$(MAKE) USE_MPI=1


.PHONY: mpi_debug
mpi_debug:
	$(MAKE) USE_MPI=1 FFLAGS='$(DBG_FFLAGS)'


.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)