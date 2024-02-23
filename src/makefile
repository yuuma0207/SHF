.PHONY: all move_mods clean
TARGET  = main
MODULES = constants_and_parameters.f90 global_variables.f90 math_derivation.f90 spacial_function.f90 harmonic_occ.f90 math_integrate.f90 write_to_file.f90 solve_poisson_problem.f90 determine_meanfield_potential.f90 determine_observables.f90
OBJECTS = constants_and_parameters.o global_variables.o math_derivation.o spacial_function.o harmonic_occ.o math_integrate.o write_to_file.o solve_poisson_problem.o determine_meanfield_potential.o determine_observables.o main.o
F90 	= gfortran
FFLAGS  = -Wall -fcheck=all -fall-intrinsics -fimplicit-none -fbacktrace
OPENMP  = -fopenmp
LDFLAGS = -llapack -lblas
SRC_DIR = ./../src
MOD_DIR = ./../mod
OBJ_DIR = ./../obj
ZUNDA_DIR = ./../zunda

all: ${TARGET} move_mods move_objs

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o:
	${F90} -g -c $< ${FFLAGS} ${OPENMP}

${TARGET} : ${OBJECTS}
	${F90} -g -o $@ ${OBJECTS} ${OPENMP} ${LDFLAGS}

${OBJECTS} : ${MODULES}


move_mods:
	mkdir -p $(MOD_DIR)
	mv *.mod $(MOD_DIR)
	
move_objs:
	mkdir -p $(OBJ_DIR)
	mv *.o $(OBJ_DIR)
	rm $(OBJ_DIR)/main.o

clean:
	$(RM) $(OBJECTS) $(TARGET) *.mod
