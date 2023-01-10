############################
# COMPILER AND DEBUG MODE
############################
#
#
#
ifdef OPT
	DEBUG  := false
else
	DEBUG  := true
endif
#
#
#
# COMP := ifort OR gfortran OR idris
# Choose automatically the compiler
hostname=$(shell hostname)
#
ifeq ($(hostname),turing1)
	COMP  := idris
else
	ifeq ($(hostname),turing2)
	   COMP  := idris
	else
	   COMP  := ifort
	endif
endif
#
COMP  := gfortran
#ulimit -s unlimited
#
$(info Hostname is $(hostname))
$(info Compiler is $(COMP))
#
#
#
############################
# SUBDIRECTORIES : USE OF VPATH
############################
#
#
#
SUB_DIR    := . chemkin eglib tranfit thd numerical_recipes cosiness postcosiness io
SRC_DIR    := $(addprefix src/,$(SUB_DIR))
VPATH      := $(SRC_DIR)
BUILD_DIR  := build
MODS_DIR   := mods
HDF5       := /home/rb/Desktop/codes/devel-petsc/linux-mpich-gcc
HDF5_OPTC  := -I${HDF5}/include
HDF5_OPTL  := -L${HDF5}/lib  -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz
OPTC       := ${HDF5_OPTC}
OPTL       := ${HDF5_OPTL}

#BUILD_DIR  := /workgpfs/rech/coh/rcoh940/run_2048/
#
#
#
############################
# BINARY FILES
############################
#
#
#
COSINESS_BIN      := creams
POSTCOSINESS_BIN  := postcreams
#
#
#
############################
# COMPILATION FLAGS
############################
#
#
#
ifeq ($(DEBUG),true)
#
	ifeq ($(COMP),ifort)
	    F77       := ifort
	    F90       := mpif90
	    F77FLAGS  := -cpp -traceback -O3
	    F90FLAGS  := -cpp -traceback -g -warn alignments -warn interfaces -warn unused -check all -CB -ftrapuv -mp1
#	    F90FLAGS  := -cpp -traceback -g -warn alignments -warn interfaces -warn unused -check noarg_temp_created -CB -ftrapuv -mp1  # Hide temporary array creation messages
	else
	    ifeq ($(COMP),gfortran)
	        F77       := /home/rb/Desktop/codes/devel-petsc/linux-mpich-gcc/bin/mpif77
	        F90       := /home/rb/Desktop/codes/devel-petsc/linux-mpich-gcc/bin/mpif90
	        F77FLAGS  := -cpp -Wintrinsics-std -O3
	        F90FLAGS  := -cpp -fbacktrace -g -fbounds-check -ffpe-trap=invalid,zero,overflow -mieee-fp -Wall
	    else
	        F77       := mpixlf77_r
	        F90       := mpixlf90_r
	        F77FLAGS  := -qfixed -O3 -qstrict -qsuffix=f=f90:cpp=f
	        F90FLAGS  := -qfree=f90 -qstrict -qsuffix=f=f90:cpp=f90
	    endif
	endif
#
else
#
	ifeq ($(COMP),ifort)
	    F77       := ifort
	    F90       := mpif90
	    F77FLAGS  := -cpp -O3
	    F90FLAGS  := -cpp -O3
	else
	    ifeq ($(COMP),gfortran)
	        F77       := /home/rb/Desktop/codes/devel-petsc/linux-mpich-gcc/bin/mpif77
	        F90       := /home/rb/Desktop/codes/devel-petsc/linux-mpich-gcc/bin/mpif90
	        F77FLAGS  := -cpp -O3
	        F90FLAGS  := -cpp -O3 -msse3 -march=native -mtune=native -funroll-loops -ftree-vectorize
	    else
	        F77       := mpixlf77_r
	        F90       := mpixlf90_r
	        F77FLAGS  := -qfixed -O3 -qstrict -qsuffix=f=f90:cpp=f
	        F90FLAGS  := -qfree=f90 -O3 -qstrict -qsuffix=f=f90:cpp=f90
	    endif
	endif
#
endif
#
#
#
############################
# SOURCES
############################
#
#
#
COSINESS_SRC := nrtype.f90  nrutil.f90  module_nr_subroutines.f90  random.f90  \
	          module_parameters.f90  \
	          module_parallel.f90  \
	          module_adim.f90  \
	          param_thd.f90  type_thd.f90  common_file.f90  file_tools.f90  thd_tools_perfectgas.f90  \
	          module_deriv.f90  \
	          module_input.f90  \
	          module_thermodynamics.f90  \
	          module_Rankine_Hugoniot.f90  \
	          module_weno.f90  \
	          module_tools.f90  \
	          module_iohdf5.f90 \
	          module_SGS_models.f90 \
	          module_eg_lib.f90  \
	          module_viscflux.f90  \
	          module_BCs.f90  \
	          module_ICs.f90  \
	          module_reaction.f90  \
	          module_solver.f90  \
	          ckinterp36.f  tranfit.f  \
	          ckinterp39.f  xerror.f  tran.f  \
	          initchemkin.f  vode.f  psr_chemkin.f  psr_simp_chemkin.f  reaction_v.f  \
	          egfrmc.f  EGSlib.f  EGini.f  \
	          cosiness.f90
#
#
#
POSTCOSINESS_SRC := nrtype.f90  nrutil.f90  module_nr_subroutines.f90  random.f90  \
                  module_parameters.f90  \
	          module_parallel.f90  \
	          module_adim.f90  \
	          param_thd.f90  type_thd.f90  common_file.f90  file_tools.f90  thd_tools_perfectgas.f90  \
	          module_deriv_post.f90  \
	          module_input.f90  \
	          module_variables.f90  \
	          module_thermodynamics.f90  \
	          module_Rankine_Hugoniot.f90  \
	          module_weno.f90  \
	          module_tools.f90  \
	          module_SGS_models.f90 \
	          module_BCs.f90  \
	          module_reaction.f90  \
	          module_eg_lib.f90  \
	          module_tools_post.f90  \
	          module_budget_tke.f90  module_budget_vorticity.f90  module_budget_scal_var.f90  module_budget_Reynolds.f90 \
	          module_autocorrelations.f90  \
	          module_inst_plot.f90  \
	          module_stat_plot.f90  \
	          module_statistics.f90  \
	          module_paraview.f90  \
	          ckinterp36.f  tranfit.f  \
	          ckinterp39.f  xerror.f  tran.f  \
	          initchemkin.f  vode.f  psr_chemkin.f  psr_simp_chemkin.f  psr_simp_chemkin_qi.f  reaction_v.f  \
	          egfrmc.f  EGSlib.f  EGini.f  \
	          postcosiness.f90
#
#
#
############################
# F77 AND F90 OBJECTS
############################
#
#
#
COSINESS_OBJ     :=  $(subst .f90,.o,$(COSINESS_SRC)) # F90 OBJECTS
COSINESS_OBJ     :=  $(subst .f,.o,$(COSINESS_OBJ))   # F90 & F77 OBJECTS
#
#
#
POSTCOSINESS_OBJ :=  $(subst .f90,.o,$(POSTCOSINESS_SRC)) # F90 OBJECTS
POSTCOSINESS_OBJ :=  $(subst .f,.o,$(POSTCOSINESS_OBJ))   # F90 & F77 OBJECTS
#
#
#
############################
# TARGETS
############################
#
#
#
$(COSINESS_BIN): $(COSINESS_OBJ)
	       $(F90) $(COSINESS_OBJ) ${OPTL} -o  $(BUILD_DIR)/$(COSINESS_BIN)
	       mv *.o *.mod $(MODS_DIR)
#
#
#
$(POSTCOSINESS_BIN): $(POSTCOSINESS_OBJ)
	           $(F90)  $(POSTCOSINESS_OBJ) ${OPTL} -o  $(BUILD_DIR)/$(POSTCOSINESS_BIN)
	           mv *.o *.mod $(MODS_DIR)

#
#
#
all:
	make $(COSINESS_BIN)
	make $(POSTCOSINESS_BIN)
#
#
#
clean:
	rm -f *.o *~ *.mod *mod.f90
#
#
#
clean_mod:
	rm -f *.mod
#
#
#
clean_all:
	$(foreach bdir,$(SRC_DIR), rm -f $(bdir)/*~)  \
	rm -f *.o *~ *.mod *mod.f90 $(BUILD_DIR)/*~  \
	$(BUILD_DIR)/*.bin $(BUILD_DIR)/cosiness \
	$(BUILD_DIR)/postcosiness
#
#
#
############################
# SUFFIXES
############################
#
#
#
%.o: %.f
	${F77} -c ${F77FLAGS} $<
%.o: %.f90
	${F90} -c ${F90FLAGS} ${OPTC} $<
