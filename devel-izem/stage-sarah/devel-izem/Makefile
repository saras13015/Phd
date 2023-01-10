############################
# COMPILER AND DEBUG MODE
############################
#
ifdef OPT
	DEBUG  := false
else
	DEBUG  := true
endif
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
#
$(info Hostname is $(hostname))
$(info Compiler is $(COMP))
#
############################
# SUBDIRECTORIES : USE OF VPATH
############################
#
SUB_DIR    := .  chemkin  eglib  tranfit  thd  numerical_recipes
SRC_DIR    := $(addprefix src/,$(SUB_DIR))
VPATH      := $(SRC_DIR)
BUILD_DIR  := bin
#
############################
# BINARY FILES
############################
#
IZEM_BIN      := izem.exe
POSTIZEM_BIN  := postizem.exe
#
############################
# COMPILATION FLAGS
############################
#
HDF5_HOME := /home/rb/Desktop/codes/devel-petsc/linux-mpich-gcc
HDF5_INC  := -I$(HDF5_HOME)/include -DHDF5
LDFLAGS   := -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lhdf5_hl -lz
F77       := /home/rb/Desktop/codes/devel-petsc/linux-mpich-gcc/bin/mpif77
F90       := /home/rb/Desktop/codes/devel-petsc/linux-mpich-gcc/bin/mpif90
# F77FLAGS  := -cpp -O3 -msse3 -march=native -mtune=native -funroll-loops -ftree-vectorize $(HDF5_INC) #-Wno-argument-mismatch
# F90FLAGS  := -cpp -O3 -msse3 -march=native -mtune=native -funroll-loops -ftree-vectorize $(HDF5_INC) #-Wno-argument-mismatch

F77FLAGS  := -cpp -fbacktrace -ffixed-form -O3 $(HDF5_INC) #-Wno-argument-mismatch
F90FLAGS  := -cpp -fbacktrace -g -fbounds-check -ffpe-trap=invalid,zero,overflow -mieee-fp -Wall $(HDF5_INC) #-Wno-argument-mismatch


# ifeq ($(DEBUG),true)
# #
# 	ifeq ($(COMP),ifort)
# 	    F77       := ifort
# 	    F90       := mpif90
# 	    F77FLAGS  := -cpp -traceback -O3
# 	    F90FLAGS  := -cpp -traceback -g -warn alignments -warn interfaces  -check all -CB -ftrapuv -mp1
# 	else
# 	    ifeq ($(COMP),gfortran)
# 	        F77       := mpif77 
# 	        F90       := mpif90
# 	        F77FLAGS  := -cpp -fbacktrace -ffixed-form -O3
# 	        F90FLAGS  := -cpp -fbacktrace -g -fbounds-check -ffpe-trap=invalid,zero,overflow -mieee-fp -Wall
# 	    else
# 	        F77       := mpixlf77_r
# 	        F90       := mpixlf90_r
# 	        F77FLAGS  := -qfixed -O3 -qstrict -qsuffix=f=f90:cpp=f
# 	        F90FLAGS  := -qfree=f90 -qstrict -qsuffix=f=f90:cpp=f90
# 	    endif
# 	endif
# #
# else
# #
# 	ifeq ($(COMP),ifort)
# 	    F77       := ifort
# 	    F90       := mpif90
# 	    F77FLAGS  := -cpp -O3
# 	    F90FLAGS  := -cpp -O3
# 	else
# 	    ifeq ($(COMP),gfortran)
# 	        F77       := ifort
# 	        F90       := mpif90
# 	        F77FLAGS  := -cpp -ffixed-form -O3
# 	        F90FLAGS  := -cpp -O3 -msse3 -march=native -mtune=native -funroll-loops -ftree-vectorize
# 	    else
# 	        F77       := mpixlf77_r
# 	        F90       := mpixlf90_r
# 	        F77FLAGS  := -qfixed -O3 -qstrict -qsuffix=f=f90:cpp=f
# 	        F90FLAGS  := -qfree=f90 -O3 -qstrict -qsuffix=f=f90:cpp=f90
# 	    endif
# 	endif
# #
# endif
#
#
#
############################
# SOURCES
############################
#
#
#
IZEM_SRC :=nrtype.f90  \
			nrutil.f90  \
			module_nr_subroutines.f90  \
			random.f90  \
			module_parameters.f90  \
			module_parallel.f90  \
			module_hdf5.f90   \
			module_adim.f90  \
			module_parsing.f90  \
			module_input.f90  \
			param_thd.f90  \
			type_thd.f90 \
			common_file.f90 \
			file_tools.f90  \
			thd_tools_perfectgas.f90  \
			module_thermodynamics.f90  \
			module_Rankine_Hugoniot.f90  \
			module_weno.f90  \
			module_deriv.f90  \
			module_tools.f90  \
			module_SGS_models.f90 \
			module_eg_lib.f90  \
			module_viscflux.f90  \
			module_BCs.f90  \
			module_ICs.f90  \
			module_reaction.f90  \
			module_solver.f90  \
			ckinterp36.f \
			tranfit.f  \
			ckinterp39.f \
			xerror.f \
			tran.f  \
			initchemkin.f  \
			vode.f  \
			psr_chemkin.f \
			psr_simp_chemkin.f \
			reaction_v.f  \
			egfrmc.f  \
			EGSlib.f  \
			EGini.f \
			module_io.f90 \
			izem.f90
#
POSTIZEM_SRC :=nrtype.f90  \
				nrutil.f90  \
				module_nr_subroutines.f90  \
				random.f90  \
				module_parameters.f90  \
				module_parallel.f90  \
				module_hdf5.f90  \
				module_adim.f90  \
				module_parsing.f90           \
				module_input.f90  \
				param_thd.f90  type_thd.f90  common_file.f90  file_tools.f90  thd_tools_perfectgas.f90  \
				module_variables.f90  \
				module_thermodynamics.f90  \
				module_Rankine_Hugoniot.f90  \
				module_weno.f90  \
				module_deriv_post.f90  \
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
				module_io_post.f90 \
				postizem.f90

#
############################
# F77 AND F90 OBJECTS
############################
#
#
#
IZEM_OBJ     :=  $(subst .f90,.o,$(IZEM_SRC)) # F90 OBJECTS
IZEM_OBJ     :=  $(subst .f,.o,$(IZEM_OBJ))   # F90 & F77 OBJECTS
#
#
#
POSTIZEM_OBJ :=  $(subst .f90,.o,$(POSTIZEM_SRC)) # F90 OBJECTS
POSTIZEM_OBJ :=  $(subst .f,.o,$(POSTIZEM_OBJ))   # F90 & F77 OBJECTS
#
#
#
############################
# TARGETS
############################
#
#
#
$(IZEM_BIN): $(IZEM_OBJ)
	       $(F90)  $(IZEM_OBJ) $(LDFLAGS) -o  $(BUILD_DIR)/$(IZEM_BIN)
#
#
#
$(POSTIZEM_BIN): $(POSTIZEM_OBJ)
	           $(F90)  $(POSTIZEM_OBJ) $(LDFLAGS) -o  $(BUILD_DIR)/$(POSTIZEM_BIN)
#
#
#
all:
	make $(IZEM_BIN)
	make $(POSTIZEM_BIN)
#
#
#
clean:
	rm -f *.o *~ *.mod *mod.f90
#
#
veryclean:
	$(foreach bdir,$(SRC_DIR), rm -f $(bdir)/*~)  \
	rm -f *.o *~ *.mod *mod.f90 $(BUILD_DIR)/*~  \
	$(BUILD_DIR)/*.bin $(BUILD_DIR)/*.exe \
	$(BUILD_DIR)/*.out
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
	${F90} -c ${F90FLAGS} $<

