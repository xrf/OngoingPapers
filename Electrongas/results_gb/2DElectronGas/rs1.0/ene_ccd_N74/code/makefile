
#
# Compiler options:
#         'srl' for serial job
#         'prl' for parallel MPI job
#
CMP_OPTION = 'prl'

# Flag options:
#         'normal' for normal run
#         'mkl_abel' for normal run with MKL on Abel
#         'prof' for profiling
#         'idbc' for debugging with idbc
#         'valgrind' for memory checks with Valgrind
#         'check' to check the code
#
#         Start the multiprocessing daemon: mpd --ncpus=4 &
#         End mpd: mpdallexit
#
#         Run parallel:
#         mpirun -np 4 ./ccd.exe input.dat < /dev/null > outlist &
#
#         Run valgrind: 
#         valgrind --leak-check=full --show-reachable=yes --log-file="valgrind_out" ./ccd.exe input.dat
#
#         Run gprof:
#         ./ccd.exe input.dat
#         gprof ./ccd.exe gmon.out > profile
#
FLG_OPTION = 'mkl_abel'

# Link line options:
#         'mkl_laptop' when using MKL, static linking, 
#                      Fortran 95 interface and parallel
#                      routines on my laptop
#         'st_blas' for standard BLAS and LAPACK 
#         'goto' for GotoBLAS on Titan
#         'mkl_abel' for MKL on Abel
#
LN_OPTION = 'mkl_abel'


ifeq ($(CMP_OPTION), 'srl')
CMP = ifort 
else 
CMP = mpif90
endif

MODPATH = ./modules

ifeq ($(FLG_OPTION), 'normal')
FLG = -O2 -fp-model strict -heap-arrays 1024 -openmp -module ${MODPATH} 
else ifeq ($(FLG_OPTION), 'mkl_abel')
FLG = -O2 -fp-model strict -heap-arrays 1024 -openmp -module ${MODPATH} \
	-mkl=parallel
else ifeq ($(FLG_OPTION), 'prof')
FLG = -O0 -pg -module ${MODPATH} 
else ifeq ($(FLG_OPTION), 'idbc')
FLG = -heap_arrays -debug -O0 -module ${MODPATH} 
else ifeq ($(FLG_OPTION), 'valgrind')
FLG = -g -O0 -module ${MODPATH} 
else ifeq ($(FLG_OPTION), 'check')
FLG = -g -check all -trace -O3 -module ${MODPATH} 
endif

ifeq ($(LN_OPTION), 'mkl_laptop')
LN = -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include \
	-I$(MKLROOT)/include/intel64/lp64 -lmkl_blas95_lp64 \
	-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
	$(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 \
	-lpthread
else ifeq ($(LN_OPTION), 'st_blas')
LN = -llapack -lblas
else ifeq ($(LN_OPTION), 'goto')
LN = -L/site/VERSIONS/Goto2/ -limf -lm -llapack -lgoto
else ifeq ($LN_OPTION), 'mkl_abel')
LN = 
endif

#     Object files
#
OBJ = Constants.o Mpi.o MemoryCounter.o Library.o Library2.o \
	Matrix.o IntMatrix.o RealMatrix.o MatList.o \
	IntMatList.o Quadrature.o SpBasis.o SpBasisPW2d.o \
	SpBasisPW3d.o PairList.o Block1.o BlockList.o \
	BlockListRe.o Hamiltonian.o Yukawa2dHam.o \
	Yukawa3dHam.o FockMatrix.o RefEnergy.o \
	CcdAmplitude.o CorrEnergy.o Pt2CorrEnergy.o \
	CcdCorrEnergy.o Energy.o Solver.o Tests.o Main.o

ccd.exe: ${OBJ}
	${CMP} ${FLG} -o ccd.exe ${OBJ} ${LN}  

%.o: 		%.f90
		${CMP} ${FLG} -c $<

%.o: 		%.f
		${CMP} ${FLG} -c -static $<

clean : 
	rm -f *.o ${MODPATH}/*.mod *.exe


