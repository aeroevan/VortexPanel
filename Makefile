FC=mpif90
#FC=ifort
FCOPTS=-O3 -fopenmp -mtune=native -march=native
#FCOPTS=-fast -openmp -parallel

COMMON_OBJ=matrixtools.o workingprecision.o

all: project_qr project_cg project_cg_single

clean:
	rm -f *.mod
	rm -f *.o
	rm -f project_qr project_cg project_cg_single

test_qr.o: test_qr.f90 $(COMMON_OBJ) qr.o
	$(FC) $(FCOPTS) -c test_qr.f90

test_qr: test_qr.o $(COMMON_OBJ) qr.o
	$(FC) $(FCOPTS) -o test_qr test_qr.o $(COMMON_OBJ) qr.o

project_qr.o: project_qr.f90 $(COMMON_OBJ) qr.o
	$(FC) $(FCOPTS) -c project_qr.f90

project_qr: project_qr.o $(COMMON_OBJ) qr.o
	$(FC) $(FCOPTS) -o project_qr project_qr.o $(COMMON_OBJ) qr.o

project_cg.o: project_cg.f90 $(COMMON_OBJ) cg.o
	$(FC) $(FCOPTS) -c project_cg.f90

project_cg_single.o: project_cg_single.f90 $(COMMON_OBJ) cg.o
	$(FC) $(FCOPTS) -c project_cg_single.f90

project_cg: project_cg.o $(COMMON_OBJ) cg.o
	$(FC) $(FCOPTS) -o project_cg project_cg.o $(COMMON_OBJ) cg.o

project_cg_single: project_cg_single.o $(COMMON_OBJ) cg.o
	$(FC) $(FCOPTS) -o project_cg_single project_cg_single.o $(COMMON_OBJ) cg.o

workingprecision.o: workingprecision.f90
	$(FC) $(FCOPTS) -c workingprecision.f90

matrixtools.o: matrixtools.f90 workingprecision.o
	$(FC) $(FCOPTS) -c matrixtools.f90

qr.o: qr.f90 $(COMMON_OBJ)
	$(FC) $(FCOPTS) -c qr.f90

cg.o: cg.f90 $(COMMON_OBJ)
	$(FC) $(FCOPTS) -c cg.f90
