CMP = mpif90
MODPATH = ./modules
FLG = -O2 -openmp
LN = -llapack -lblas
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


