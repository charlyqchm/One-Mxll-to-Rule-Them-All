FC=gfortran
# FFLAGS=-g -fopenmp -ffpe-trap=invalid,zero,overflow -fbacktrace
FFLAGS=-O3 -cpp -fopenmp

SRC:=constants_mod.f90
SRC+=input_mod.f90
SRC+=rs_vec_base_mod.f90
SRC+=rs_vec_dimensions_mod.f90
SRC+=allocator_multidim_mod.f90
SRC+=medium_mod.f90
SRC+=target_source_mod.f90
SRC+=parallel_subs_mod.f90
SRC+=rs_operations_subs_mod.f90
SRC+=operators_mod.f90
SRC+=design_mod.f90
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
