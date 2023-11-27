SOURCE	= kuo_schemes.f90 thermo_functions2.f90

FC = gfortran
#FFLAGS = -fdefault-real-8 -O3 
#FFLAGS = -g -fdefault-real-8 -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace \
#-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -frecursive
FFLAGS = -g -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace \
-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

OBJECTS	=	$(SOURCE:.f90=.o)

.SUFFIXES:
.SUFFIXES:	.o .f90

.f90.o :
	$(FC)  $(FFLAGS) -c $*.f90

main:	$(OBJECTS)
	$(FC) -o main ${OBJECTS} 

clean:
	 \rm -f $(OBJECTS) *~

fort:
	 \rm -f  fort.* *~
