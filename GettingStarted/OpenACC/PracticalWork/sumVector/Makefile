COMP=gcc

ifeq ($(COMP),gcc)
FC=gfortran
#FFLAGS=-g -O3 -fopenacc -fopt-info-note-omp
FFLAGS=-g -O3
RUN=
endif

EXE=sum

all: $(EXE)

$(EXE): sum.F90
	$(FC) $(FFLAGS) $^ -o $@ 

clean:
	rm -f $(EXE) *.o *~

