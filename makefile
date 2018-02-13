# Makefile to build lbstime.a library and driver routines 
#
# Andrew J. Pounds, Ph.D.
# Departments of Chemistry and Computer Science
# Mercer University
# Spring 2007
#

# Use gfortran to avoid g77/gcc incompatibilities on Fedora
# g77 works fine on systems where gcc/g77 were linked against the same library

# with this makefile if you type just "make" it will compile the code with the default
# dimension of 50.  If you type "make size=100" it will compile with dimension 100 or whatever
# other size you want to use

F77 = gfortran    
F90 = gfortran
CC  = gcc 
TIMINGLIB_FLAGS = -O3
OPT_FLAGS = -fcompare-elim -fipa-pure-const -fsplit-wide-types -ftree-ccp -ftree-ch -DOPT1 -DOPT2 -DOPT3 -DOPT4 -DOPT5 -DOPT6 -fcombine-stack-adjustments -fipa-reference -fshrink-wrap -fssa-backprop -ftree-coalesce-vars -ftree-copy-prop

DEF_FLAGS = -DDO_TIMING

TIMINGLIBS =  -L./ -llbstime 
CLIBS = -lm

OBJS = cputime.o walltime.o  

ifeq ($(origin size), undefined)
	PROFLAGS = -D CPPFLAGS=100
else
	PROFLAGS = -D CPPFLAGS=$(size)
endif

all: lib dusty77 dusty90 dusty 

cputime.o : cputime.cc   
	$(CC) $(TIMINGLIB_FLAGS) -c cputime.cc  

walltime.o : walltime.cc   
	$(CC) $(TIMINGLIB_FLAGS) -c walltime.cc  

dusty77.o : dusty77.f   
	$(F77) -cpp $(PROFLAGS) $(DEF_FLAGS) $(OPT_FLAGS) -c dusty77.f   

# Don't forget the -lstdc++
dusty77 : dusty77.o  $(OBJS) 
	$(F77) -o dusty77 dusty77.o  $(TIMINGLIBS) -lstdc++  

dusty90.o : dusty90.f90  
	$(F90) -cpp $(PROFLAGS) $(DEF_FLAGS) $(OPT_FLAGS) -c dusty90.f90   

dusty90 : dusty90.o  $(OBJS) 
	$(F90) -o dusty90 dusty90.o  $(TIMINGLIBS) -lstdc++  

dusty.o : dusty.c   
	$(CC) $(PROFLAGS) $(DEF_FLAGS) $(OPT_FLAGS) -c dusty.c   

# Don't forget the -lstdc++
dusty : dusty.o  $(OBJS) 
	$(CC) -o dusty dusty.o  $(TIMINGLIBS) $(CLIBS) -lstdc++   

# Default Targets for Cleaning up the Environment
clean :
	rm *.o
	rm *.a

pristine :
	rm *.o
	rm *.a
	touch *.c *.f *.f90 

ctags :
	ctags  *.c *.f *.f90

# Target for making the library

lib: $(OBJS) 
	ar -rc liblbstime.a $(OBJS) 
	ranlib liblbstime.a

run :
	dusty77 | tee testdata77.dat
	dusty77 | tee -a testdata77.dat
	dusty90 | tee testdata90.dat
	dusty90 | tee -a testdata90.dat
	dusty | tee testdata.dat
	dusty | tee -a testdata.dat

average :
	echo "mydata = read.table(\"testdata77.dat\", header=FALSE)" >mean77.R
	echo "mean(mydata[[\"V2\"]])" >>mean77.R
	R -q --vanilla -f mean77.R
	rm mean77.R
	echo "mydata = read.table(\"testdata90.dat\", header=FALSE)" >mean90.R
	echo "mean(mydata[[\"V2\"]])" >>mean90.R
	R -q --vanilla -f mean90.R
	rm mean90.R
	echo "mydata = read.table(\"testdata.dat\", header=FALSE)" >mean.R
	echo "mean(mydata[[\"V2\"]])" >>mean.R
	R -q --vanilla -f mean.R
	rm mean.R

