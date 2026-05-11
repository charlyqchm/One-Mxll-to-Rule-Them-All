FC ?= gfortran
MPI_FC ?= mpif90
USE_MPI ?= 0
DEBUG ?= 0

# Build selectors: by default, compile both independent programs.
FDTD ?= 1
DESIGN ?= 1

BUILD ?= build
BIN_DIR ?= $(BUILD)/bin

FDTD_DIR := src/FDTD
DESIGN_DIR := src/freq_domain_inverse_design

SELECTED_TARGETS :=
ifeq ($(FDTD),1)
SELECTED_TARGETS += fdtd
endif
ifeq ($(DESIGN),1)
SELECTED_TARGETS += design
endif

ifeq ($(strip $(SELECTED_TARGETS)),)
$(error At least one build target must be enabled: use FDTD=1 and/or DESIGN=1)
endif

.PHONY: all fdtd design mpi debug mpi_debug clean clean_fdtd clean_design

all: $(SELECTED_TARGETS)

fdtd:
	$(MAKE) -C $(FDTD_DIR) \
		FC='$(FC)' MPI_FC='$(MPI_FC)' \
		USE_MPI='$(USE_MPI)' DEBUG='$(DEBUG)' \
		BUILD='../../$(BUILD)' BIN_DIR='../../$(BIN_DIR)'

design:
	$(MAKE) -C $(DESIGN_DIR) \
		FC='$(FC)' MPI_FC='$(MPI_FC)' \
		USE_MPI='$(USE_MPI)' DEBUG='$(DEBUG)' \
		BUILD='../../$(BUILD)' BIN_DIR='../../$(BIN_DIR)'

mpi:
	$(MAKE) USE_MPI=1 FDTD='$(FDTD)' DESIGN='$(DESIGN)'

debug:
	$(MAKE) DEBUG=1 FDTD='$(FDTD)' DESIGN='$(DESIGN)'

mpi_debug:
	$(MAKE) USE_MPI=1 DEBUG=1 FDTD='$(FDTD)' DESIGN='$(DESIGN)'

clean: clean_fdtd clean_design
	rm -rf $(BUILD)

clean_fdtd:
	$(MAKE) -C $(FDTD_DIR) clean BUILD='../../$(BUILD)' BIN_DIR='../../$(BIN_DIR)'

clean_design:
	$(MAKE) -C $(DESIGN_DIR) clean BUILD='../../$(BUILD)' BIN_DIR='../../$(BIN_DIR)'